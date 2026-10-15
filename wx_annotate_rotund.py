#!/usr/bin/python2.5

from math import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *

try: from pydb import debugger
except ImportError:
    import pdb
    def debugger(): pdb.set_trace()



import wx
from wxBasicGLCanvas import BasicGLCanvas
from ui_actions import *
from wx_util import Timers
from mesh_gui import TriMesh, Camera, set_headlight, draw_mesh_faces_flat, draw_mesh_faces_smooth
from helpers import Struct, solve_for_fixed_points, smooth_open_stroke, closest_point_index, contour_contains

from vecutil import *
from OpenGL.GL import *
from glhelpers import *

kDebugging = False
from vecutil import kEps
#from annotations_rotund import kDebugging, kEps

kUseTexture = True
kOnlyShowAxis0 = True
## Screen-space (pixels) picking distance threshold
kPickScreenThreshold = 20
kAutoPlaceInDepth = True
kUseOblique = True

draw_mesh_faces = draw_mesh_faces_flat

from annotations_rotund import RotundRegion, Connection, mesh_llconstraints_from_constraints2, mesh_llconstraints_from_constraints2_oblique
import fibermesh


class RotundAnnotations( object ):
    def __init__( self ):
        self.Reset()
    
    def Reset( self, keep_texture = False ):
        self.rotund_regions = []
        self.connections = []
        self.oblique = [0.,0.]
        if not keep_texture:
            self.image_path = None
            self.gltex = None
    
    def RotundAnnotationsFromFile( filename ):
        result = RotundAnnotations()
        result.Load( filename )
        return result
    RotundAnnotationsFromFile = staticmethod( RotundAnnotationsFromFile )
    
    def SetImage( self, image_path ):
        from GLTexture import GLTexture
        self.image_path = image_path
        try:
            self.gltex = GLTexture( image_path )
        except IOError:
            print('Failed to load texture:', image_path)
            self.image_path = None
            return
        
        ## Set the rectangle's size according to the texture's aspect ratio.
        self.tex_size = arrayf( ( self.gltex.width, self.gltex.height ) )
        self.tex_size *= 2./self.tex_size.max()
        
        ## Center the image.
        self.tex_pos = zeros(2) - .5*self.tex_size
    
    FFORMAT_V1 = 'wx_annotate_rotund file format v1'
    FFORMAT_V2 = 'wx_annotate_rotund file format v2'
    
    def Save( self, filename, overwrite = False ):
        import os
        if not overwrite and os.path.exists( filename ):
            print('File exists, not overwriting:', filename)
            return
        
        if overwrite and os.path.isdir( filename ):
            print('Directory with same name as desired file already exists, can\'t overwrite:', filename)
            return
        
        ## If we made it here, we're fine.
        if overwrite and os.path.exists( filename ):
            print('Overwriting file:', filename)
        
        ## This is our file format:
        fformat = ( self.FFORMAT_V2, self.rotund_regions, self.connections, self.image_path, self.oblique )
        
        import pickle
        f = open( filename, 'wb' )
        pickle.dump( fformat, f )
    
    def Load( self, filename ):
        import os
        if not os.path.isfile( filename ):
            print('Can\'t load from non-existent file:', filename)
            return
        
        self.Reset()
        import pickle
        f = open( filename, 'rb' )
        ## Specify encoding = 'latin1' to load Python 2 pickles with Python 3.
        ## https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3/28218598#28218598
        fformat = pickle.load( f, encoding = 'latin1' )
        
        if fformat[0] not in (self.FFORMAT_V1,self.FFORMAT_V2):
            print('Wrong file format, cannot load.')
            return
        
        self.rotund_regions = fformat[1]
        ## Fixup for some changes to the RotundRegion data structure
        for region in self.rotund_regions:
            if not hasattr( region, '_mesh' ): region._mesh = region.__dict__['mesh']
            if not hasattr( region, 'last_iterate_diff' ): region.last_iterate_diff = RotundRegion([],[],[]).last_iterate_diff
            if not hasattr( region, 'iterate_threshold' ): region.iterate_threshold = RotundRegion([],[],[]).iterate_threshold
        
        self.connections = fformat[2]
        
        self.SetImage( fformat[3] )
        
        if fformat[0] == self.FFORMAT_V1:
            self.oblique = [0., 0.]
        else:
            self.oblique = fformat[4]
    
    def DeleteRegion( self, region ):
        self.rotund_regions.remove( region )
        
        ## Also delete any connections which reference this region.
        connections_delete_me = []
        for connection in self.connections:
            if region in connection.regions:
                connections_delete_me.append( connection )
        for con in connections_delete_me: self.connections.remove( con )
    
    def PlaceInDepth( self ):
        from annotations_rotund import place_in_depth_according_to_connection
        
        frozen_regions = set()
        
        ## A map from a region to connections touching it.
        region2connections = {}
        for region in self.rotund_regions: region2connections.setdefault( region, [] )
        for connection in self.connections:
            for region in connection.regions:
                region2connections[ region ].append( connection )
        
        def dfs( region ):
            assert region in frozen_regions
            
            for connection in region2connections[ region ]:
                assert len( connection.regions ) == 2
                if connection.regions[0] in frozen_regions and connection.regions[1] in frozen_regions: continue
                place_in_depth_according_to_connection( connection, frozen_regions, self.oblique )
                dfs( [ child for child in connection.regions if child != region ][0] )
        
        for region in self.rotund_regions:
            if region not in frozen_regions:
                frozen_regions.add( region )
                dfs( region )

class RotundAnnotateCanvasBase( BasicGLCanvas ):
    
    EYE_DISTANCE = 100.
    
    def __init__( self, parent, annotations, pos = wx.DefaultPosition, size = wx.DefaultSize, shared = None ):
        BasicGLCanvas.__init__( self, parent, pos, size, debugging = True, shared = shared )
        
        self.camera = Camera()
        ## Move the default camera EYE_DISTANCE away.
        self.camera.eye[2] = self.EYE_DISTANCE
        self.camera.near_far[1] = 1000.
        
        self.annotations = annotations
        self.show_annotations = True
        
        if kUseTexture:
            import os
            #self.annotations.SetImage( '/mixed media/pictures/ColossalKatamari.png' )
            #self.annotations.SetImage( os.path.join( 'scenes', 'cactus.png' ) )
            self.annotations.SetImage( os.path.join( 'scenes', 'artwork', 'crow bakery - rabbit.png' ) )
        
        self.debugging = False
    
    def GetBestSize( self ):
        return (500,500)
    
    def prepare2D( self ):
        glDisable( GL_DEPTH_TEST )
        glDisable( GL_LIGHTING )
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        ## Center the square inside the viewport.
        glOrtho(*( list( self.world_bounds ) + [ -100, 100 ] ))
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
    
    def prepare3D( self ):
        glEnable( GL_DEPTH_TEST )
        glEnable( GL_LIGHTING )
        
        ## Keep the up vector pointed in positive y, so we don't get twisted.
        self.camera.up[:] = (0,1,0)
        ## Center the square inside the viewport.
        if kUseOblique:
            self.camera.proj = Camera.Oblique( *(
                list( self.world_bounds ) + list( self.annotations.oblique ) + [self.EYE_DISTANCE]
                ) )
        else:
            self.camera.proj = Camera.Ortho( *self.world_bounds )
        self.camera.Apply()
        ## TODO: Something with the light.  Right now, we're getting OpenGL's default light.
        
        set_headlight( self.camera )
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE )
        glShadeModel( GL_SMOOTH )


class ModeProxy( EventHandler ):
    def __init__( self, canvas ):
        self.canvas = canvas
    
    ## This magic allows non-overloaded member names to magically pass through to self.canvas.
    def __getattr__( self, name ):
        if not name.startswith('_'): return getattr( self.canvas, name )
        else: raise AttributeError("I am not proxying attribute '%s'" % (name,))

class ObliqueMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def doMouseDown( self, evt ):
        return self.doMouseDrag( evt )
    
    def doMouseDrag( self, evt ):
        
        world_pos = self.event2world( self.GetMousePosition() )
        self.annotations.oblique[0] = -world_pos[0]
        self.annotations.oblique[1] = -world_pos[1]
        
        print('oblique:', self.annotations.oblique)
        
        self.GetTopLevelParent().Refresh()
        return True
    
    def doDraw( self ):
        def draw_cube():
            self.prepare2D()
            
            glColor3f( 0,0,0 )
            glLineWidth( 8 )
            
            glMatrixMode( GL_MODELVIEW )
            glPushMatrix()
            glTranslatef( -1, -1, 0 )
            
            a,b = self.annotations.oblique
            
            ## Front face of cube
            draw_lineloop( [
                [0,0,0],
                [1,0,0],
                [1,1,0],
                [0,1,0]
                ] )
            ## Back face of cube
            draw_lineloop( [
                [0 - a,0 - b,-1],
                [1 - a,0 - b,-1],
                [1 - a,1 - b,-1],
                [0 - a,1 - b,-1]
                ] )
            ## Sides
            draw_lines( [
                [0,0,0],
                [0 - a,0 - b,-1],
                [1,0,0],
                [1 - a,0 - b,-1],
                [1,1,0],
                [1 - a,1 - b,-1],
                [0,1,0],
                [0 - a,1 - b,-1]
                ] )
            
            glPopMatrix()
        
        self.DrawScene( callback = draw_cube )

class AnnotateViewRotateMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.rotator = RotateHandler( self.camera, self.Refresh, self.event2world, debugging = self.debugging )
    
    def doMouseDown( self, evt ):
        return self.rotator.doMouseDown( evt )
    
    def doMouseDrag( self, evt ):
        self.rotator.doMouseDrag( evt )
        return True
    
    def doMouseUp( self, evt ):
        self.camera.center = arrayf( ( 0,0,0 ) )
        self.camera.eye = arrayf( ( 0,0,self.EYE_DISTANCE ) )
        self.camera.up = arrayf( ( 0,1,0 ) )
    
    def doDraw( self ):
        self.DrawScene()

class AnnotationCreateMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.stroke = None
        self.stroke_handler = None
    
    def doMouseDown( self, evt ):
        ## Forward events to the stroke drawing class.
        self.stroke = []
        if self.stroke_handler: unforward_events( self.stroke_handler )
        self.stroke_handler = StrokeHandler(
            self.stroke,
            lambda pt: self.Refresh(),
            self.doStrokeHandlerFinished,
            self.event2world,
            debugging = self.debugging
            )
        forward_events( self, self.stroke_handler )
        ## Now that we've forwarded events, call doMouseDown() again.  This time it will forward.
        return self.doMouseDown( evt )
    
    def doStrokeHandlerFinished( self, stroke ):
        self.stroke = None
        unforward_events( self.stroke_handler )
        self.stroke_handler = None
        
        self.NewRegionBoundaryStrokeFinished( stroke )
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
        ## in wxPython.
        if keycode in (wx.WXK_DELETE,8):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            try:
                region = self.PickRegion( screen_pos )
                print('deleting a region')
                self.annotations.DeleteRegion( region )
                
                if kAutoPlaceInDepth: self.annotations.PlaceInDepth()
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        else:
            evt.Skip()
        
        return self
    
    def doDraw( self ):
        def draw_stroke():
            if self.stroke is None: return
            
            glColor4ub( *( self.Color[ 'silhouette' ] + (200,) ) )
            glLineWidth( 8 )
            draw_linestrip( self.stroke )
        
        self.DrawScene( callback = draw_stroke )

class CrossSectionDragMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.stroke = None
        self.stroke_handler = None
    
    def doMouseDown( self, evt ):
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return self
        
        screen_pos = self.event2screen( self.GetMousePosition() )
        region = self.PickRegion( screen_pos )
        if region is None: return False
        
        ## Check if we pick an existing cross section.
        csi = self.PickCrossSection( screen_pos, region )
        if csi is not -1:
            self.stroke = [ region.cross_sections[ csi ][0] ]
            if self.stroke_handler: unforward_events( self.stroke_handler )
            self.stroke_handler = StrokeHandler(
                self.stroke,
                #lambda pt: self.Refresh(),
                lambda pt: self.doStrokeHandlerMoved( region, csi, self.stroke ),
                lambda stroke: self.doStrokeHandlerFinished( region, csi, stroke ),
                self.event2world,
                debugging = self.debugging
                )
            forward_events( self, self.stroke_handler )
            ## Now that we've forwarded events, call doMouseDown() again.  This time it will forward.
            mdresult = self.doMouseDown( evt )
            
            ## doMouseDown() clears the stroke, but I really want it to keep the first point in the
            ## cross section.
            del self.stroke[:]
            self.stroke.append( region.cross_sections[ csi ][0] )
            return mdresult
        
        
        ## Forward events to the stroke drawing class.
        self.stroke = []
        if self.stroke_handler: unforward_events( self.stroke_handler )
        region.cross_sections.append( None )
        self.stroke_handler = StrokeHandler(
            self.stroke,
            #lambda pt: self.Refresh(),
            lambda pt: self.doStrokeHandlerMoved( region, -1, self.stroke ),
            lambda stroke: self.doStrokeHandlerFinished( region, -1, stroke ),
            self.event2world,
            debugging = self.debugging
            )
        forward_events( self, self.stroke_handler )
        ## Now that we've forwarded events, call doMouseDown() again.  This time it will forward.
        return self.doMouseDown( evt )
    
    def doStrokeHandlerMoved( self, region, csi, stroke ):
        if len( stroke ) > 1 and stroke[0] != stroke[-1]:
            region.cross_sections[csi] = ( stroke[0], stroke[-1] )
            
            if kUseOblique:
                mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2_oblique(
                    region, self.annotations.oblique
                    )
            else:
                mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2(
                    region
                    )
            expect( len( xyz_cross_sections ) > 0 )
            
            #region.mesh = mesh
            region.llconstraints = llconstraints
            region.cross_sections3d = xyz_cross_sections
            
            self.GetTopLevelParent().Refresh()
    
    def doStrokeHandlerFinished( self, region, csi, stroke ):
        self.stroke = None
        unforward_events( self.stroke_handler )
        self.stroke_handler = None
        
        if len( stroke ) > 1 and stroke[0] != stroke[-1]:
            region.cross_sections[csi] = ( stroke[0], stroke[-1] )
            
            if kUseOblique:
                mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2_oblique(
                    region, self.annotations.oblique
                    )
            else:
                mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2(
                    region
                    )
            expect( len( xyz_cross_sections ) > 0 )
            
            #region.mesh = mesh
            region.llconstraints = llconstraints
            region.cross_sections3d = xyz_cross_sections
            
            self.GetTopLevelParent().Refresh()
        else:
            pass
            #del region.cross_sections[csi]
    
    def doContextMenu( self, evt ):
        evt_pos = self.GetMousePosition()
        #print 'context menu:', self.event2screen( evt_pos )
        
        screen_pos = self.event2screen( evt_pos )
        region = self.PickRegion( screen_pos )
        if region is None:
            evt.Skip()
            return self
        
        ## Check if we pick an existing cross section.
        csi = self.PickCrossSection( screen_pos, region )
        if csi == -1:
            evt.Skip()
            return self
        
        
        ## Show a contextual menu for the cross section
        menu = wx.Menu()
        menu_ids = []
        
        def delete( evt ):
            del region.cross_sections[csi]
            del region.cross_sections3d[csi]
            
            mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2(
                region
                )
            
            #region.mesh = mesh
            region.llconstraints = llconstraints
            region.cross_sections3d = xyz_cross_sections
            
            self.GetTopLevelParent().Refresh()
        
        item = menu.Append( wx.ID_ANY, "Delete" )
        self.Bind( wx.EVT_MENU, delete, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        self.PopupMenu( menu )
        menu.Destroy()
        for id in menu_ids: self.Unbind( wx.EVT_MENU, id = id )
        
        return self
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
        ## in wxPython.
        if keycode in (wx.WXK_DELETE,8):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            try:
                region = self.PickRegion( screen_pos )
                print('deleting a region')
                self.annotations.DeleteRegion( region )
                
                if kAutoPlaceInDepth: self.annotations.PlaceInDepth()
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        else:
            evt.Skip()
        
        return self
    
    def doDraw( self ):
        def draw_stroke():
            if self.stroke is None or len( self.stroke ) < 2: return
            
            glColor4ub( *( self.Color[ 'cross_section' ] + (200,) ) )
            glLineWidth( 8 )
            draw_linestrip( [ self.stroke[0], self.stroke[-1] ] )
        
        self.DrawScene( callback = draw_stroke )

class ConnectionCurveDrawMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.stroke = None
        self.stroke_handler = None
    
    def doMouseDown( self, evt ):
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return self
        
        screen_pos = self.event2screen( self.GetMousePosition() )
        ## If we don't pick exactly two regions under the mouse, don't draw a connection.
        regions = self.PickRegions( screen_pos )
        if len( regions ) != 2: return False
        ## If there's already a connection between the two regions, overdraw it.
        ## TODO Q: Should this drag the connection curve instead?
        previous_connection = [ con for con in self.annotations.connections if set( con.regions ) == set( regions ) ]
        assert len( previous_connection ) <= 1
        previous_connection = None if len( previous_connection ) == 0 else previous_connection[0]
        
        ## TODO: Check if we would create a cycle in the connection graph; if so, bail.
        
        ## Forward events to the stroke drawing class.
        self.stroke = []
        if self.stroke_handler: unforward_events( self.stroke_handler )
        self.stroke_handler = StrokeHandler(
            self.stroke,
            lambda pt: self.Refresh(),
            lambda stroke: self.doStrokeHandlerFinished( regions, previous_connection, stroke ),
            self.event2world,
            debugging = self.debugging
            )
        forward_events( self, self.stroke_handler )
        ## Now that we've forwarded events, call doMouseDown() again.  This time it will forward.
        return self.doMouseDown( evt )
    
    def doStrokeHandlerFinished( self, regions, previous_connection, stroke ):
        assert len( regions ) == 2
        
        self.stroke = None
        unforward_events( self.stroke_handler )
        self.stroke_handler = None
        
        if len( stroke ) > 1 and tuple( stroke[0] ) != tuple( stroke[-1] ):
            
            stroke = smooth_open_stroke( stroke, resample_N = 10 )
            if previous_connection:
                previous_connection.ReplaceConnectionCurve( stroke )
            else:
                connection = Connection( stroke, regions[0], True, regions[1], True )
                self.annotations.connections.append( connection )
            
            if kAutoPlaceInDepth: self.annotations.PlaceInDepth()
            self.GetTopLevelParent().Refresh()
    
    def doContextMenu( self, evt ):
        evt_pos = self.GetMousePosition()
        #print 'context menu:', self.event2screen( evt_pos )
        
        screen_pos = self.event2screen( evt_pos )
        regions = self.PickRegions( screen_pos )
        #print 'context menu picked regions:', regions
        ## If we don't pick exactly two regions, bail.
        if len( regions ) != 2:
            evt.Skip()
            return self
        
        connection = None
        for con in self.annotations.connections:
            if set( con.regions ) == set( regions ):
                connection = con
                break
        #print 'context menu picked connection:', connection
        ## If the regions don't have a connection between them, bail.
        if connection is None:
            evt.Skip()
            return self
        
        
        ## Show a contextual menu for the closest two overlapping annotations.
        menu = wx.Menu()
        menu_ids = []
        
        def disconnect( evt ):
            self.annotations.connections.remove( connection )
            
            if kAutoPlaceInDepth: self.annotations.PlaceInDepth()
            self.GetTopLevelParent().Refresh()
        
        item = menu.Append( wx.ID_ANY, "Disconnect" )
        self.Bind( wx.EVT_MENU, disconnect, id = item.GetId() )
        menu_ids.append( item.GetId() )
        menu.AppendSeparator()
        
        def connect_with_signature( signature ):
            connection.regions_attaches_front = signature
            
            if kAutoPlaceInDepth: self.annotations.PlaceInDepth()
            self.GetTopLevelParent().Refresh()
        
        def connect_front_front( evt ):
            print('connect_front_front')
            connect_with_signature( (True,True) )
        def connect_back_back( evt ):
            print('connect_back_back')
            connect_with_signature( (False,False) )
        
        item = menu.Append( wx.ID_ANY, "Connect, Front" )
        self.Bind( wx.EVT_MENU, connect_front_front, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        item = menu.Append( wx.ID_ANY, "Connect, Front/Back" )
        self.Bind( wx.EVT_MENU, lambda evt: connect_with_signature( (True,False) ), id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        item = menu.Append( wx.ID_ANY, "Connect, Back/Front" )
        self.Bind( wx.EVT_MENU, lambda evt: connect_with_signature( (False,True) ), id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        item = menu.Append( wx.ID_ANY, "Connect, Back" )
        self.Bind( wx.EVT_MENU, connect_back_back, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        
        self.PopupMenu( menu )
        menu.Destroy()
        for id in menu_ids: self.Unbind( wx.EVT_MENU, id = id )
        
        return self
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
        ## in wxPython.
        if keycode in (wx.WXK_DELETE,8):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            try:
                region = self.PickRegion( screen_pos )
                print('deleting a region')
                self.annotations.DeleteRegion( region )
                
                if kAutoPlaceInDepth: self.annotations.PlaceInDepth()
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        else:
            evt.Skip()
        
        return self
    
    def doDraw( self ):
        def draw_stroke():
            if self.stroke is None or len( self.stroke ) < 2: return
            
            glColor4ub( *( self.Color[ 'sharp_cut' ] + (200,) ) )
            glLineWidth( 8 )
            draw_linestrip( self.stroke )
        
        self.DrawScene( callback = draw_stroke )

class ConnectionCurveDragMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.connection = None
    
    def doMouseDown( self, evt ):
        screen_pos = self.event2screen( self.GetMousePosition() )
        self.connection = self.PickConnection( screen_pos )
        if self.connection is None: return False
    
    def doMouseDrag( self, evt ):
        if self.connection is None: return False
        
        raise NotImplementedError
        ## 1 shoot a ray into both meshes.
        ## 2 move one region in depth so that the intersection point has the same depth on both regions.
        ## 3 update connection_curve_precise (and make sure place_in_depth_according_to_connection()
        ##   uses connection_curve_precise if it exists).
        
        self.GetTopLevelParent().Refresh()
    
    def doMouseUp( self, evt ):
        if self.connection is None: return False
        
        raise NotImplementedError
        ## 1 shoot a ray into both meshes.
        ## 2 move one region in depth so that the intersection point has the same depth on both regions.
        ## 3 update connection_curve_precise (and make sure place_in_depth_according_to_connection()
        ##   uses connection_curve_precise if it exists).
        
        self.GetTopLevelParent().Refresh()
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
        ## in wxPython.
        if keycode in (wx.WXK_DELETE,8):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            try:
                region = self.PickRegion( screen_pos )
                print('deleting a region')
                self.annotations.DeleteRegion( region )
                
                if kAutoPlaceInDepth: self.annotations.PlaceInDepth()
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        else:
            evt.Skip()
        
        return self
    
    def doDraw( self ):
        self.DrawScene()

class RotundAnnotateCanvas( RotundAnnotateCanvasBase ):
    
    def __init__( self, *args, **kwargs ):
        RotundAnnotateCanvasBase.__init__( self, *args, **kwargs )
        
        self.SetMode( AnnotationCreateMode( self ) )
    
    def SetMode( self, mode ):
        ## Oops: I can't use forward_events because it's broken and only works to one level.
        ## My modes need to be able to forward events, too.
        #if self.mode is not None: unforward_events( self.mode )
        self.mode = mode
        #forward_events( self, self.mode )
    
    kAnnotationCreateMode = 1
    kCrossSectionDragMode = 2
    kConnectionCurveDrawMode = 3
    kObliqueMode = 4
    kAnnotateViewRotateMode = 5
    
    @staticmethod
    def GetAllModeIDs():
        return [
            RotundAnnotateCanvas.kAnnotationCreateMode,
            RotundAnnotateCanvas.kCrossSectionDragMode,
            RotundAnnotateCanvas.kConnectionCurveDrawMode,
            RotundAnnotateCanvas.kObliqueMode,
            RotundAnnotateCanvas.kAnnotateViewRotateMode
            ]
    
    @staticmethod
    def GetLongNameForModeID( ID ):
        return {
            RotundAnnotateCanvas.kAnnotationCreateMode: 'Create shapes',
            RotundAnnotateCanvas.kCrossSectionDragMode: 'Drag cross sections',
            RotundAnnotateCanvas.kConnectionCurveDrawMode: 'Draw connection curves',
            RotundAnnotateCanvas.kObliqueMode: 'Specify the oblique camera',
            RotundAnnotateCanvas.kAnnotateViewRotateMode: 'Rotate the camera, temporarily'
            }[ ID ]
    
    @staticmethod
    def GetShortNameForModeID( ID ):
        return {
            RotundAnnotateCanvas.kAnnotationCreateMode: 'Create',
            RotundAnnotateCanvas.kCrossSectionDragMode: 'Tilt',
            RotundAnnotateCanvas.kConnectionCurveDrawMode: 'Connect',
            RotundAnnotateCanvas.kObliqueMode: 'Oblique',
            RotundAnnotateCanvas.kAnnotateViewRotateMode: 'Rotate'
            }[ ID ]
    
    def SetModeByID( self, ID ):
        mode_type = {
            RotundAnnotateCanvas.kAnnotationCreateMode: AnnotationCreateMode,
            RotundAnnotateCanvas.kCrossSectionDragMode: CrossSectionDragMode,
            RotundAnnotateCanvas.kConnectionCurveDrawMode: ConnectionCurveDrawMode,
            RotundAnnotateCanvas.kObliqueMode: ObliqueMode,
            RotundAnnotateCanvas.kAnnotateViewRotateMode: AnnotateViewRotateMode
            }[ ID ]
        print('Switching to mode:', self.GetLongNameForModeID( ID ))
        self.SetMode( mode_type( self ) )
    
    def PickCrossSection( self, screen_pos, region ):
        world_pos = self.screen2world( screen_pos )
        
        assert len( region.cross_sections ) == len( region.cross_sections3d ) or debugger()
        
        for csi in range(len( region.cross_sections3d )):
            if contour_contains( [ pt[:2] for pt in region.cross_sections3d[ csi ] ], world_pos ):
                return csi
        
        return -1
    
    def PickRegions( self, screen_pos ):
        from surface_intersections import ray_mesh_intersections
        
        world_pos = self.screen2world( screen_pos )
        ## Shoot the ray in the negative screen-space z direction.
        ## (Invert the (oblique) projection matrix to get the actual ray.)
        ray = ( append( world_pos, 0 ), ( self.annotations.oblique[0], self.annotations.oblique[1], -1) )
        
        t2region = []
        
        for region in self.annotations.rotund_regions:
            if region.mesh is None: continue
            
            isects = ray_mesh_intersections( ray, region.mesh )
            if len( isects ) == 0: continue
            
            minisect = min( isects )
            t2region.append( ( minisect[0], region ) )
        
        t2region.sort()
        result = [ region for t, region in t2region ]
        #print 'PickRegions( %s ):' % (screen_pos,), t2region
        return result
    
    def PickRegion( self, screen_pos ):
        regions = self.PickRegions( screen_pos )
        if len( regions ) == 0: return None
        else: return regions[0]
    
    def NewRegionBoundaryStrokeFinished( self, stroke ):
        
        ## Reset if not enough points in the stroke.
        if len( stroke ) < 3: return
        
        import fibermesh
        
        region_boundary = fibermesh.smooth_curve( stroke, smooth = 11 )
        #print 'self.data.silhouette:', self.data.silhouette
        free_vertices = set([
            vi
            for vi, pos in enumerate( region_boundary )
            if closest_point_index( pos, stroke ) in (0,len(stroke)-1)
            ])
        ## Make all edges without an endpoint in free_vertices silhouette_edges.
        #silhouette_edges = range(len( region_boundary ))
        silhouette_edges = [
            ei for ei in range(len( region_boundary ))
            if ei not in free_vertices and ((ei+1) % len(region_boundary)) not in free_vertices
            ]
        
        fibermesh.save_stroke( region_boundary, 'last_stroke.py' )
        
        #region = RotundRegion( region_boundary, range(len( region_boundary )), [] )
        region = RotundRegion( region_boundary, silhouette_edges, [] )
        ## TODO: For debugging, append a no-op cross section constraint, see if anything happens.
        #region.cross_sections.append( ( anchor, anchor ) )
        
        if kUseOblique: #self.annotations.oblique != [0., 0.]:
            mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2_oblique(
                region, self.annotations.oblique
                )
        else:
            mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2(
                region
                )
        expect( len( xyz_cross_sections ) == 0 )
        
        region.mesh = mesh
        region.llconstraints = llconstraints
        region.cross_sections3d = xyz_cross_sections
        
        self.annotations.rotund_regions.append( region )
        
        ## TODO: place in depth
        
        ## Refresh all
        self.GetTopLevelParent().Refresh()
    
    Color = {
        'rotund_region': ( 200, 200, 255 ),
        'silhouette': ( 234, 145, 248 ),
        'sharp_cut': ( 174//2, 237//2, 111//2 ),
        'cross_section': ( 213, 99, 191 ),
        'spine': ( 74, 74, 74 )
        }
    
    def DrawRotundRegion( self, region, what_to_draw = None ):
        if what_to_draw is None: what_to_draw = list(self.Color.keys())
        
        self.prepare3D()
        glDisable( GL_LIGHTING )
        
        if 'rotund_region' in what_to_draw:
            ## Draw the interior
            ## TODO Q: Should I use draw it 2d, using gluTessPolygon (see my layers2d app for reference)
            if region.mesh is not None:
                glColor4ub( *( self.Color['rotund_region'] + (50,) ) )
                ## draw_mesh_faces_smooth() is much faster than draw_mesh_faces_flat()
                draw_mesh_faces_smooth( region.mesh, True )
        
        if 'silhouette' in what_to_draw:
            if region.mesh is not None:
                glLineWidth( 8 )
                glColor4ub( *( self.Color['silhouette'] + (200,) ) )
                
                ## TODO: Use weakref.WeakKeyDictionary()
                if not hasattr( region.mesh, '_silhouette_edges' ) or region.mesh._silhouette_edges[0] != region.mesh.lifetime_counter:
                    from surface_intersections import mesh_silhouette_edges
                    mesh_skew = TriMesh()
                    mesh_skew.faces = array( region.mesh.faces )
                    skew = identity(3)
                    skew[0,2] = self.annotations.oblique[0]
                    skew[1,2] = self.annotations.oblique[1]
                    mesh_skew.vs = dot( skew, asarray( region.mesh.vs ).T ).T
                    mesh_skew.topology_changed()
                    region.mesh._silhouette_edges = region.mesh.lifetime_counter, mesh_silhouette_edges( mesh_skew )
                
                pts = []
                for edge in region.mesh._silhouette_edges[1]:
                    pts.append( region.mesh.vs[ edge[0] ] )
                    pts.append( region.mesh.vs[ edge[1] ] )
                draw_lines( pts )
        
        ## Draw the boundary stroke.
        assert region.region_boundary is not None
        if region.region_boundary is not None:
            self.prepare2D()
            
            if 'silhouette' in what_to_draw:
                ## Draw silhouette edges
                glLineWidth( 8 )
                glColor4ub( *( self.Color['silhouette'] + (50,) ) )
                #draw_lineloop( region.region_boundary )
                draw_lines( asarray( [
                    ( region.region_boundary[ ei ], region.region_boundary[ ( ei + 1 ) % len( region.region_boundary ) ] )
                    for ei in region.silhouette_edges
                    ] ).reshape( ( 2*len( region.silhouette_edges ), 2 ) ) )
            
            if 'sharp_cut' in what_to_draw:
                ## Draw sharp cut edges
                glColor4ub( *( self.Color['sharp_cut'] + (50,) ) )
                draw_lines( asarray( [
                    ( region.region_boundary[ ei ], region.region_boundary[ ( ei + 1 ) % len( region.region_boundary ) ] )
                    for ei in region.sharp_cut_edges
                    ] ).reshape( ( 2*len( region.sharp_cut_edges ), 2 ) ) )
        
        
        if 'cross_section' in what_to_draw:
            ## Draw cross sections.
            if region.cross_sections3d:
                self.prepare3D()
                glDisable( GL_LIGHTING )
                
                glDisable( GL_DEPTH_TEST )
                ## This doesn't work so well at making only the front of the cross section visible,
                ## because the cross section is independent of the geometry it's intersecting, so
                ## we might need quite a bit of an offset.
                ## Instead, we'll disable the depth test and set the alpha to the z of each point.
                #glMatrixMode( GL_MODELVIEW )
                #glPushMatrix()
                #glTranslated( *( -.05 * append( self.annotations.oblique, -1 ) ) )
                
                glLineWidth( 3 )
                #glColor4ub( *( self.Color['cross_section'] + (200,) ) )
                #for cs3d in region.cross_sections3d: draw_lineloop( cs3d )
                
                for cs3d in region.cross_sections3d:
                    colors = empty( ( len(cs3d), 4 ), dtype = ubyte )
                    colors[:,0:3] = self.Color['cross_section']
                    
                    cs3d = asarray( cs3d )
                    minz = cs3d[:,2].min()
                    maxz = cs3d[:,2].max()
                    if maxz == minz: colors[:,3] = 200
                    else: colors[:,3] = 200 * ( (cs3d[:,2] - minz) / (maxz-minz) )
                    
                    glEnableClientState( GL_COLOR_ARRAY )
                    glColorPointer( len(colors[0]), GL_UNSIGNED_BYTE, 0, asarray( colors ) )
                    
                    ## OK, I know this is implemented by with glVertexPointer, but it was
                    ## until a moment ago implemented using GL_LINE_STRIP and tacking the first
                    ## point onto the end.  Better to just put the four lines of code here,
                    ## to ensure consistency.
                    #draw_lineloop( cs3d )
                    glEnableClientState( GL_VERTEX_ARRAY )
                    glVertexPointer( len(cs3d[0]), GL_DOUBLE, 0, asarray( cs3d ) )
                    glDrawArrays( GL_LINE_LOOP, 0, len(cs3d) )
                    glDisableClientState( GL_VERTEX_ARRAY )
                    
                    glDisableClientState( GL_COLOR_ARRAY )
                
                #glPopMatrix()
        
        if 'spine' in what_to_draw:
            ## Draw spine.
            if region.spine is not None:
                self.prepare2D()
                
                glLineWidth( 4 )
                glColor4ub( *( self.Color['spine'] + (200,) ) )
                draw_linestrip( region.spine )
                
                glEnable( GL_DEPTH_TEST )
    
    def DrawConnection( self, connection, alpha = 200 ):
        glLineWidth( 8 )
        glColor4ub( *( self.Color['sharp_cut'] + (alpha,) ) )
        if connection.connection_curve_precise is not None:
            self.prepare3D()
            glDisable( GL_LIGHTING )
            glDisable( GL_DEPTH_TEST )
            draw_linestrip( connection.connection_curve_precise )
        else:
            self.prepare2D()
            draw_linestrip( connection.connection_curve )
    
    def DrawScene( self, callback = None ):
        glClearColor( .9,.9,.9,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        ## Draw the gltex if there is one.
        if self.annotations.gltex is not None:
            self.prepare2D()
            
            glColor4f( 1,1,1, .5 )
            draw_gltexture( self.annotations.gltex, self.annotations.tex_pos, self.annotations.tex_size )
        
        if self.show_annotations:
            ## First draw the rotund regions, which sets the zbuffer
            for region in self.annotations.rotund_regions:
                self.DrawRotundRegion( region, what_to_draw = 'rotund_region' )
            
            ## Now draw the decorating lines
            decs = set( self.Color.keys() )
            decs.remove( 'rotund_region' )
            for region in self.annotations.rotund_regions:
                self.DrawRotundRegion( region, what_to_draw = decs )
            
            ## Draw the connections either lightly or darkly, depending on whether we're in
            ## connecting mode
            connection_alpha = 200 if type( self.mode ) is ConnectionCurveDrawMode else 50
            for connection in self.annotations.connections:
                self.DrawConnection( connection, alpha = connection_alpha )
        
        if callback is not None: callback()
        
        self.SwapBuffers()
    
    def doDraw( self ):
        return self.mode.doDraw()
    
    def doMouseDown( self, evt ):
        return self.mode.doMouseDown( evt )
    
    def doMouseUp( self, evt ):
        return self.mode.doMouseUp( evt )
    
    def doMouseMove( self, evt ):
        return self.mode.doMouseMove( evt )
    
    def doMouseDrag( self, evt ):
        return self.mode.doMouseDrag( evt )
    
    def doKeyDown( self, evt ):
        return self.mode.doKeyDown( evt )
    
    def doKeyUp( self, evt ):
        return self.mode.doKeyUp( evt )
    
    def doContextMenu( self, evt ):
        return self.mode.doContextMenu( evt )
    
    def _AddPerfectRegion( self, region_boundary, silhouette_edges = None, sharp_cut_edges = None, spine = None ):
        ## NOTE: I can't make [] the default parameter, because this always uses the same instance,
        ##       so the instance is shared by all objects instantiated with the default parameter,
        ##       and when I append to that instance, it is seen by all objects sharing it.
        if silhouette_edges is None: silhouette_edges = []
        if sharp_cut_edges is None: sharp_cut_edges = []
        
        import fibermesh
        
        fibermesh.save_stroke( region_boundary, 'last_stroke.py' )
        
        region = RotundRegion( region_boundary, silhouette_edges, sharp_cut_edges, spine )
        if kUseOblique: #self.annotations.oblique != [0.,0.]:
            mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2_oblique(
                region, self.annotations.oblique
                )
        else:
            mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2(
                region
                )
        expect( len( xyz_cross_sections ) == 0 )
        
        region.mesh = mesh
        region.llconstraints = llconstraints
        region.cross_sections3d = xyz_cross_sections
        
        self.annotations.rotund_regions.append( region )
        self.GetTopLevelParent().Refresh()
    
    def AddPerfectCircleRegion( self ):
        
        region_boundary = [ ( .5 * cos( theta ), .5 * sin( theta ) ) for theta in linspace( 0, 2*pi, 100 ) ]
        
        silhouette_edges = list(range(len( region_boundary )))
        self._AddPerfectRegion( region_boundary, silhouette_edges )
    
    def AddOvalRegion( self ):
        
        region_boundary = [ ( cos( theta ), .5 * sin( theta ) ) for theta in linspace( 0, 2*pi, 100 ) ]
        ## Scale by .9
        region_boundary = ( .9 * asarrayf( region_boundary ) ).tolist()
        
        silhouette_edges = list(range(len( region_boundary )))
        self._AddPerfectRegion( region_boundary, silhouette_edges )
    
    def AddRoundedCylinderRegion( self ):
        
        region_boundary = (
            [ ( .5 * cos( theta ) - .5, .5 * sin( theta ) ) for theta in linspace( pi/2, 3*pi/2, 100 ) ] +
            list(zip( linspace( -.5, .5, 10 ), -.5 * ones(10) ))[1:-1] +
            [ ( .5 * cos( theta ) + .5, .5 * sin( theta ) ) for theta in linspace( -pi/2, pi/2, 100 ) ] +
            list(zip( linspace( .5, -.5, 10 ), .5 * ones(10) ))[1:-1]
            )
        ## Scale by .9
        region_boundary = ( .9 * asarrayf( region_boundary ) ).tolist()
        
        silhouette_edges = list(range(len( region_boundary )))
        self._AddPerfectRegion( region_boundary, silhouette_edges )

class ResultCanvas( BasicGLCanvas ):
    
    EYE_DISTANCE = 100.
    
    def __init__( self, parent, annotations, pos = wx.DefaultPosition, size = wx.DefaultSize, shared = None ):
        ## GLUT = True because this uses draw_text().
        BasicGLCanvas.__init__( self, parent, pos, size, debugging = True, shared = shared, GLUT = True )
        
        self.camera = Camera()
        ## Center the square inside the viewport.
        self.camera.proj = Camera.Ortho( *self.world_bounds )
        self.camera.eye[2] = self.EYE_DISTANCE
        self.camera.near_far[1] = 1000.
        
        self.annotations = annotations
        self.image_depth = -.5
        self.show_vertex_numbers = False
        self.show_meshes = True
        self.show_background_image = True
        
        debugging = False
        self.rotator = RotateHandler( self.camera, self.Refresh, self.event2world, debugging = debugging )
        forward_events( self, self.rotator )
    
    def GetBestSize( self ):
        return (500,500)
    
    def doDraw( self ):
        glClearColor( .5,.5,.5,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        self.prepare3D()
        
        #self.drawTexture()
        ## Clear the depth buffer.
        ## Actually, this looks weird.
        # glClear( GL_DEPTH_BUFFER_BIT )
        
        for region in self.annotations.rotund_regions:
            if region.region_boundary is not None:
                glDisable( GL_LIGHTING )
                
                ## Draw silhouette edges.
                glLineWidth( 3 )
                glColor4ub( *( RotundAnnotateCanvas.Color['silhouette'] + (200,) ) )
                #draw_lineloop( region.region_boundary )
                draw_lines( asarray( [
                    ( region.region_boundary[ ei ], region.region_boundary[ ( ei + 1 ) % len( region.region_boundary ) ] )
                    for ei in region.silhouette_edges
                    ] ).reshape( ( 2*len( region.silhouette_edges ), 2 ) ) )
                
                ## Draw sharp cut edges
                glColor4ub( *( RotundAnnotateCanvas.Color['sharp_cut'] + (200,) ) )
                draw_lines( asarray( [
                    ( region.region_boundary[ ei ], region.region_boundary[ ( ei + 1 ) % len( region.region_boundary ) ] )
                    for ei in region.sharp_cut_edges
                    ] ).reshape( ( 2*len( region.sharp_cut_edges ), 2 ) ) )
                
            
            if region.mesh is not None:
                if self.show_vertex_numbers:
                    ## Draw vertex number labels
                    glColor4f( 0,0,0, 1 )
                    for vi in range(len( region.mesh.vs )):
                        draw_text( region.mesh.vs[vi] + .01 * region.mesh.vertex_normals[vi], str(vi) )
                
                glEnable( GL_LIGHTING )
                
                glEnable( GL_COLOR_MATERIAL )
                glColor3ub( *RotundAnnotateCanvas.Color['rotund_region'] )
                if self.show_meshes:
                    draw_mesh_faces( region.mesh, False )
            
            ## Draw cross sections.
            if region.cross_sections3d:
                glDisable( GL_LIGHTING )
                glLineWidth( 3 )
                glColor4ub( *( RotundAnnotateCanvas.Color['cross_section'] + (200,) ) )
                for cs3d in region.cross_sections3d: draw_lineloop( cs3d )
            
            ## Draw spine.
            if region.spine3d is not None:
                glDisable( GL_LIGHTING )
                glDisable( GL_DEPTH_TEST )
                
                glLineWidth( 4 )
                glColor4ub( *( RotundAnnotateCanvas.Color['spine'] + (200,) ) )
                draw_linestrip( region.spine3d )
                
                glEnable( GL_DEPTH_TEST )
            
            if kDebugging:
                if hasattr( region, '_decs' ):
                    glDisable( GL_LIGHTING )
                    for dec in region._decs:
                        #print 'dec:', dec
                        if dec[0] == 'line':
                            glLineWidth( 4 )
                            [ None, None, None, glColor3f, glColor4f ][ len( dec[1] ) ]( *dec[1] )
                            draw_linestrip( dec[2:] )
                        elif dec[0] == 'point':
                            glPointSize( 8 )
                            [ None, None, None, glColor3f, glColor4f ][ len( dec[1] ) ]( *dec[1] )
                            glBegin( GL_POINTS )
                            for pt in dec[2:]:
                                [ None, None, glVertex2d, glVertex3d, glVertex4d ][ len( pt ) ]( *pt )
                            glEnd()
        
        ## Draw the gltex if there is one.
        if self.annotations.gltex is not None:
            self.drawTexture(
                self.annotations.gltex,
                self.annotations.tex_pos,
                self.annotations.tex_size
                )
        
        self.SwapBuffers()
    
    def drawTexture( self, gltex, tex_pos, tex_size ):
        if self.image_depth is None: return
        if self.show_background_image is False: return
        
        glDisable( GL_LIGHTING )
        
        glMatrixMode( GL_MODELVIEW )
        glPushMatrix()
        glTranslate( 0,0,self.image_depth )
        
        glColor4f( 1,1,1, .5 )
        draw_gltexture( gltex, tex_pos, tex_size )
        
        glPopMatrix()
    
    def prepare2D( self ):
        glDisable( GL_DEPTH_TEST )
        glDisable( GL_LIGHTING )
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(-1, 1, -1, 1, -100, 100)
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
    
    def prepare3D( self ):
        glEnable( GL_DEPTH_TEST )
        glEnable( GL_LIGHTING )
        
        ## Keep the up vector pointed in positive y, so we don't get twisted.
        self.camera.up[:] = (0,1,0)
        ## Center the square inside the viewport.
        if kUseOblique:
            self.camera.proj = Camera.Oblique( *(
                list( self.world_bounds ) + list( self.annotations.oblique ) + [self.EYE_DISTANCE]
                ) )
        else:
            self.camera.proj = Camera.Ortho( *self.world_bounds )
        self.camera.Apply()
        ## TODO: Something with the light.  Right now, we're getting OpenGL's default light.
        
        set_headlight( self.camera )
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE )
        glShadeModel( GL_SMOOTH )
    
    def reset_rotation( self ):
        #self.camera.center = arrayf( ( 0,0,0 ) )
        self.camera.eye = arrayf( ( 0,0,self.EYE_DISTANCE ) )
        self.camera.up = arrayf( ( 0,1,0 ) )
    
    def event2world( self, pos ):
        return BasicGLCanvas.event2world( self, pos )


class AnnotateFrame( wx.Frame ):
    
    RESET_BUTTON = 10
    ITERATE_BUTTON = 20
    ITERATE_10_BUTTON = 25
    ITERATE_INF_BUTTON = 30
    
    def __init__( self, parent, annotations ):
        wx.Frame.__init__( self, parent, wx.ID_ANY, 'Annotation Modeler' )
        
        self.annotations = annotations
        
        ## Add a toolbar for the modes.  Do this before adding the canvases.
        TBFLAGS = ( wx.TB_HORIZONTAL
            | wx.NO_BORDER
            | wx.TB_FLAT
            #| wx.TB_TEXT
            #| wx.TB_HORZ_LAYOUT
            )
        tb = self.CreateToolBar( TBFLAGS )
        tbsize = ( 64, 64 )
        tb.SetToolBitmapSize( tbsize )
        
        panel = wx.Panel( self, wx.ID_ANY, size=(1000,530) )
        self.sizer = wx.BoxSizer( wx.VERTICAL )
        
        ## Add the side-by-side gl views
        glviews = wx.BoxSizer( wx.HORIZONTAL )
        ## Create the AnnotateCanvas.
        canvas = RotundAnnotateCanvas( panel, self.annotations )
        canvas.Show()
        self.annotate_canvas = canvas
        canvas.SetSize( canvas.GetBestSize() )
        glviews.Add( canvas, 1, wx.EXPAND )
        ## Create the ResultCanvas.
        canvas = ResultCanvas( panel, self.annotations, shared = self.annotate_canvas.GetContext() )
        self.result_canvas = canvas
        canvas.SetSize( canvas.GetBestSize() )
        glviews.Add( canvas, 1, wx.EXPAND )
        self.sizer.Add( glviews, 1, wx.EXPAND )
        
        ## Add the buttons
        buttons = wx.BoxSizer( wx.HORIZONTAL )
        ## Reset button
        btn = wx.Button( panel, self.RESET_BUTTON, "Reset" )
        self.Bind( wx.EVT_BUTTON, self.OnClick, btn )
        btn.SetDefault()
        btn.SetSize( btn.GetBestSize() )
        buttons.Add( btn, 0, wx.ALIGN_CENTER )
        
        ## Iterate button
        btn = wx.Button( panel, self.ITERATE_BUTTON, "Iterate N" )
        self.Bind( wx.EVT_BUTTON, self.OnClick, btn )
        btn.SetDefault()
        btn.SetSize( btn.GetBestSize() )
        buttons.Add( btn, 0, wx.ALIGN_CENTER )
        ## Iterate button
        btn = wx.Button( panel, self.ITERATE_INF_BUTTON, "Iterate inf" )
        self.Bind( wx.EVT_BUTTON, self.OnClick, btn )
        btn.SetDefault()
        btn.SetSize( btn.GetBestSize() )
        buttons.Add( btn, 0, wx.ALIGN_CENTER )
        
        ## Add the buttons sizer to the parent sizer
        self.sizer.Add( buttons, 0, wx.ALIGN_BOTTOM )
        
        panel.SetSizer( self.sizer )
        self.sizer.Fit( panel )
        self.Fit()
        
        import os
        self.last_wd = os.getcwd()
        self.save_filename = None
        self.iterate_n = 1
        self.use_cpp = True
        self.use_cotangent = False
        self.use_linear = False
        #self.use_oblique = True
        
        menubar = wx.MenuBar()
        
        ##
        
        filem = wx.Menu()
        
        item = filem.Append( wx.ID_ANY, "Open Scene...\tctrl-o" )
        self.Bind( wx.EVT_MENU, self.OnFileMenuOpenScene, item )
        
        item = filem.Append( wx.ID_ANY, "Save Scene...\tctrl-s" )
        self.Bind( wx.EVT_MENU, self.OnFileMenuSaveScene, item )
        
        item = filem.Append( wx.ID_ANY, "Save Scene As...\tctrl-shift-s" )
        self.Bind( wx.EVT_MENU, self.OnFileMenuSaveSceneAs, item )
        
        filem.AppendSeparator()
        
        item = filem.Append( wx.ID_ANY, "Load Background PNG...\tctrl-l" )
        self.Bind( wx.EVT_MENU, self.OnFileMenuLoadBackground, item )
        
        filem.AppendSeparator()
        
        item = filem.Append( wx.ID_ANY, "Save Annotation Screenshot..." ) #\tctrl-3" )
        self.Bind( wx.EVT_MENU, lambda evt: self.SaveScreenShot( self.annotate_canvas, 'annotate canvas' ), item )
        
        item = filem.Append( wx.ID_ANY, "Save Result Screenshot..." ) #\tctrl-4" )
        self.Bind( wx.EVT_MENU, lambda evt: self.SaveScreenShot( self.result_canvas, 'result canvas' ), item )
        
        menubar.Append( filem, "&File" )
        
        ##
        
        modem = wx.Menu()
        mode_ids = RotundAnnotateCanvas.GetAllModeIDs()
        for mode_id in mode_ids:
            ## TODO: make these into a check item group and have only one checked.
            item = modem.Append( mode_id, RotundAnnotateCanvas.GetLongNameForModeID( mode_id ) + "\tctrl-%s" % (mode_id,) )
            self.Bind( wx.EVT_MENU, lambda evt: self.annotate_canvas.SetModeByID( evt.Id ), item )
        menubar.Append( modem, "&Mode" )
        
        ##
        
        viewm = wx.Menu()
        
        item = viewm.Append( wx.ID_ANY, "Reset 3D Result View\tctrl-r" )
        self.Bind( wx.EVT_MENU, self.OnReset3DResultView, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Result View Image Behind\tctrl-b" )
        if self.result_canvas.image_depth < 0: item.Check()
        self.Bind( wx.EVT_MENU, self.OnToggleResultViewImageDepth, item )
        
        viewm.AppendSeparator()
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show Annotations\tctrl-a" )
        if self.annotate_canvas.show_annotations: item.Check()
        self.Bind( wx.EVT_MENU, self.OnShowAnnotations, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show Background Image" )
        if self.result_canvas.show_background_image: item.Check()
        self.Bind( wx.EVT_MENU, self.OnShowBackgroundImage, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show Vertex Numbers\tctrl-n" )
        if self.result_canvas.show_vertex_numbers: item.Check()
        self.Bind( wx.EVT_MENU, self.OnShowVertexNumbers, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Smooth Shading\tctrl-h" )
        if draw_mesh_faces == draw_mesh_faces_smooth: item.Check()
        self.Bind( wx.EVT_MENU, self.OnSmoothShading, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show 3D Meshes" )
        if self.result_canvas.show_meshes: item.Check()
        self.Bind( wx.EVT_MENU, self.OnShowMeshes, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Use Oblique Camera" )
        #if self.use_oblique: item.Check()
        if kUseOblique: item.Check()
        self.Bind( wx.EVT_MENU, self.OnUseOblique, item )
        
        item = viewm.Append( wx.ID_ANY, "Refresh screen" )
        self.Bind( wx.EVT_MENU, lambda evt: self.GetTopLevelParent().Refresh(), item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show debugging annotations" )
        if kDebugging: item.Check()
        self.Bind( wx.EVT_MENU, self.OnDebugging, item )
        
        menubar.Append( viewm, "&View" )
        
        ##
        
        testm = wx.Menu()
        
        item = testm.Append( self.ITERATE_BUTTON, "Iterate N\tctrl-i" )
        self.Bind( wx.EVT_MENU, self.OnClick, item )
        
        item = testm.Append( self.ITERATE_10_BUTTON, "Iterate 10\tctrl-alt-i" )
        self.Bind( wx.EVT_MENU, self.OnClick, item )
        
        item = testm.Append( self.ITERATE_INF_BUTTON, "Iterate inf\tctrl-shift-i" )
        self.Bind( wx.EVT_MENU, self.OnClick, item )
        
        item = testm.Append( wx.ID_ANY, "Set iteration N..." )
        self.Bind( wx.EVT_MENU, self.OnSetIterationN, item )
        
        testm.AppendSeparator()
        
        item = testm.AppendCheckItem( wx.ID_ANY, "Use C++ iteration" )
        if self.use_cpp: item.Check()
        self.Bind( wx.EVT_MENU, self.OnUseCPP, item )
        
        item = testm.AppendCheckItem( wx.ID_ANY, "Use linear (bending) energy" )
        if self.use_linear: item.Check()
        self.Bind( wx.EVT_MENU, self.OnUseLinear, item )
        
        item = testm.AppendCheckItem( wx.ID_ANY, "Use cotangent weights" )
        if self.use_cotangent: item.Check()
        self.Bind( wx.EVT_MENU, self.OnUseCotangent, item )
        
        testm.AppendSeparator()
        
        item = testm.Append( wx.ID_ANY, "Regenerate constraints\tctrl-g" )
        self.Bind( wx.EVT_MENU, self.OnRegenerateConstraints, item )
        
        item = testm.Append( wx.ID_ANY, "Swap silhouette constraint type" )
        self.Bind( wx.EVT_MENU, self.OnSwapSilhouetteConstraints, item )
        
        item = testm.Append( wx.ID_ANY, "Place in depth\tctrl-d" )
        self.Bind( wx.EVT_MENU, self.OnPlaceInDepth, item )
        
        testm.AppendSeparator()
        
        item = testm.Append( wx.ID_ANY, "Draw a perfect circle" )
        self.Bind( wx.EVT_MENU, lambda evt: self.annotate_canvas.AddPerfectCircleRegion(), item )
        
        item = testm.Append( wx.ID_ANY, "Draw a perfect ellipse" )
        self.Bind( wx.EVT_MENU, lambda evt: self.annotate_canvas.AddOvalRegion(), item )
        
        item = testm.Append( wx.ID_ANY, "Draw a rounded rect" )
        self.Bind( wx.EVT_MENU, lambda evt: self.annotate_canvas.AddRoundedCylinderRegion(), item )
        
        menubar.Append( testm, "&Test" )
        
        ##
        
        self.SetMenuBar( menubar )
        
        
        ## Add modes to the toolbar
        ## TODO: make these into a check item group and have only one checked.
        mode_ids = RotundAnnotateCanvas.GetAllModeIDs()
        for mode_id in mode_ids:
            short_name = RotundAnnotateCanvas.GetShortNameForModeID( mode_id )
            long_name = RotundAnnotateCanvas.GetLongNameForModeID( mode_id )
            
            img = wx.Image( os.path.join( ".", "modes", short_name.lower() + ".png" ), wx.BITMAP_TYPE_ANY )
            img = img.Scale( tbsize[0], tbsize[1], quality = wx.IMAGE_QUALITY_HIGH )
            bitmap = img.ConvertToBitmap()
            tb.AddRadioLabelTool( mode_id, short_name, bitmap, shortHelp=short_name, longHelp=long_name )
            self.Bind( wx.EVT_TOOL, lambda evt: self.annotate_canvas.SetModeByID( evt.Id ), id=mode_id )
        
        tb.Realize()
        
        ## Do I need to call this again?
        self.Fit()
    
    def SaveScreenShot( self, canvas, suggested_name ):
        dialog = wx.FileDialog(
            self, "Save screenshot...", self.last_wd,
            suggested_name + '.png',
            "PNG|*.png|JPEG|*.jpg",
            wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
            )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            
            self.save_filename = str( dialog.GetPath() )
            canvas.SaveToFile( self.save_filename )
    
    def OnToggleResultViewImageDepth( self, evt ):
        self.result_canvas.image_depth *= -1.
        self.result_canvas.Refresh()
    
    #def OnToggleTestNewOverlap( self, evt ):
    #    global kTestNewOverlap
    #    kTestNewOverlap = not kTestNewOverlap
    
    def OnReset3DResultView( self, evt ):
        self.result_canvas.reset_rotation()
        self.result_canvas.Refresh()
    
    def OnShowAnnotations( self, evt ):
        self.annotate_canvas.show_annotations = not self.annotate_canvas.show_annotations
        self.annotate_canvas.Refresh()
    
    def OnShowBackgroundImage( self, evt ):
        self.result_canvas.show_background_image = not self.result_canvas.show_background_image
        self.result_canvas.Refresh()
    
    def OnShowVertexNumbers( self, evt ):
        self.result_canvas.show_vertex_numbers = not self.result_canvas.show_vertex_numbers
        self.result_canvas.Refresh()
    
    def OnSmoothShading( self, evt ):
        global draw_mesh_faces
        if draw_mesh_faces == draw_mesh_faces_smooth: draw_mesh_faces = draw_mesh_faces_flat
        else: draw_mesh_faces = draw_mesh_faces_smooth
        self.GetTopLevelParent().Refresh()
    
    def OnUseCPP( self, evt ):
        self.use_cpp = not self.use_cpp
    
    def OnUseLinear( self, evt ):
        self.use_linear = not self.use_linear
    
    def OnUseCotangent( self, evt ):
        self.use_cotangent = not self.use_cotangent
    
    def OnUseOblique( self, evt ):
        #self.use_oblique = not self.use_oblique
        global kUseOblique
        kUseOblique = not kUseOblique
        self.GetTopLevelParent().Refresh()
    
    def OnDebugging( self, evt ):
        global kDebugging
        kDebugging = not kDebugging
        self.GetTopLevelParent().Refresh()
    
    def OnShowMeshes( self, evt ):
        self.result_canvas.show_meshes = not self.result_canvas.show_meshes
        self.GetTopLevelParent().Refresh()
    
    def OnFileMenuOpenScene( self, evt ):
        dialog = wx.FileDialog( self, "Open scene file...", self.last_wd, '', "Annotation Scene|*.annscene", style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            
            save_filename = str( dialog.GetPath() )
            try:
                self.annotations.Load( save_filename )
            except:
                raise
            else:
                self.save_filename = save_filename
            finally:
                self.GetTopLevelParent().Refresh()
    
    def OnFileMenuLoadBackground( self, evt ):
        dialog = wx.FileDialog( self, "Load background PNG...", self.last_wd, '', "PNG|*.png", style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            
            self.annotations.SetImage( str( dialog.GetPath() ) )
            self.GetTopLevelParent().Refresh()
    
    def OnFileMenuSaveScene( self, evt ):
        if self.save_filename is not None:
            self.annotations.Save( self.save_filename, overwrite = True )
            return
        
        suggested_name = 'annotated scene'
        import os
        if self.annotations.image_path: suggested_name = os.path.split( os.path.splitext( self.annotations.image_path )[0] )[1]
        
        dialog = wx.FileDialog(
            self, "Save scene file...", self.last_wd,
            suggested_name + '.annscene',
            "Annotation Scene|*.annscene",
            wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
            )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            
            self.save_filename = str( dialog.GetPath() )
            self.annotations.Save( self.save_filename, overwrite = True )
    
    def OnFileMenuSaveSceneAs( self, evt ):
        self.save_filename = None
        self.OnFileMenuSaveScene( evt )
    
    def OnSetIterationN( self, evt ):
        dlg = wx.NumberEntryDialog( self, "", "Number of iterations:", "N",
            self.iterate_n,
            1,
            10000
            )
        if dlg.ShowModal() == wx.ID_OK:
            try:
                val = int( dlg.GetValue() )
                self.iterate_n = val
            except TypeError:
                print('Value was not a number!')
        dlg.Destroy()
    
    def OnClick( self, event ):
        if event.GetId() == self.RESET_BUTTON:
            print("reset!")
            self.annotations.Reset( keep_texture = True )
            ## Refresh all
            self.GetTopLevelParent().Refresh()
        
        elif event.GetId() == self.ITERATE_BUTTON:
            print("iterate!")
            self.IterateOnce()
        
        elif event.GetId() == self.ITERATE_10_BUTTON:
            print("iterate!")
            self.IterateOnce( 10 )
        
        elif event.GetId() == self.ITERATE_INF_BUTTON:
            if hasattr( self, 'iterate_inf_timer_id' ):
                Timers.CancelTimer( self.iterate_inf_timer_id )
                del self.iterate_inf_timer_id
            else:
                self.iterate_inf_timer_id = Timers.RequestTimer( self.IterateOnce, 10, once = False )
    
    def IterateOnce( self, N = None ):
        if N is None: N = self.iterate_n
        
        for region in self.annotations.rotund_regions:
            region.iterate( N, fast = self.use_cpp, use_cotangent = self.use_cotangent, linear = self.use_linear )
        
        ## Place all in depth
        self.annotations.PlaceInDepth()
        ## Refresh all
        self.GetTopLevelParent().Refresh()
    
    def OnRegenerateConstraints( self, event ):
        oblique = self.annotations.oblique if kUseOblique else None
        
        for region in self.annotations.rotund_regions: region.update_llconstraints( oblique )
    
    def OnSwapSilhouetteConstraints( self, event ):
        import annotations_rotund
        if annotations_rotund.kUseGenSilhouetteConstraint == annotations_rotund.gen_silhouette_constraint:
            print('Will use gen_silhouette_constraint2()')
            annotations_rotund.kUseGenSilhouetteConstraint = annotations_rotund.gen_silhouette_constraint2
        else:
            print('Will use gen_silhouette_constraint()')
            annotations_rotund.kUseGenSilhouetteConstraint = annotations_rotund.gen_silhouette_constraint
    
    def OnPlaceInDepth( self, event ):
        self.annotations.PlaceInDepth()
        self.GetTopLevelParent().Refresh()

def parse_options( argv ):
    global kDebugging
    
    try:
        argv.remove( '-d' )
        kDebugging = True
    except: pass
    try:
        argv.remove( '--debug' )
        kDebugging = True
    except: pass

def main():
    import sys
    argv = list( sys.argv[1:] )
    parse_options( argv )
    
    import wx
    #app = wx.App() ## <-- on OS X tracebacks and stdout/stderr appear in a dialog
    app = wx.PySimpleApp()
    
    ## load the PyCrust shell
    ## NOTE: This generates the following errors and is very slow switching between the app window
    ##       and the PyCrust window.
    ## Sun Jun 15 20:19:36 yotronagon.local Python[32825] <Error>: CGContextRestoreGState: invalid context
    if kDebugging:
        from wxBasicGLCanvas import load_PyCrust
        crust = load_PyCrust( app )
    
    annotations = RotundAnnotations()
    
    frame = AnnotateFrame( None, annotations )
    frame.Show(True)
    
    for fname in argv:
        if fname.lower().endswith( '.annscene' ):
            annotations.Load( fname )
            frame.Refresh()
        elif fname.lower().endswith( '.png' ):
            annotations.SetImage( fname )
            frame.Refresh()
    
    app.annotateframe = frame
    app.MainLoop()
    app.Destroy()

if __name__ == '__main__': main()
