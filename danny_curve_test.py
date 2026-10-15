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
from mesh_gui import TriMesh, Camera, set_headlight, draw_mesh_faces_flat, draw_mesh_faces_smooth, draw_mesh_edges
from helpers import Struct, solve_for_fixed_points, smooth_open_stroke, closest_point_index, contour_contains

from vecutil import *
from OpenGL.GL import *
from glhelpers import *

kDebugging = False
from vecutil import kEps
#from annotations_rotund import kDebugging, kEps

kUseTexture = False
kOnlyShowAxis0 = True
## Screen-space (pixels) picking distance threshold
kPickScreenThreshold = 20
kAutoPlaceInDepth = True
kUseOblique = True

#draw_mesh_faces = draw_mesh_faces_flat
draw_mesh_faces = draw_mesh_faces_smooth

from danny_curve_shape import CurveBoundaryRegion, mesh_llconstraints_from_closed_curve_3D, curve_3D_mesh_llconstraints_for_N_by_M_grid
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
    
    FFORMAT_V1 = 'curve_test file format v1'
    
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
        fformat = ( self.FFORMAT_V1, self.rotund_regions, self.connections, self.image_path, self.oblique )
        
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
        
        if fformat[0] not in (self.FFORMAT_V1,):
            print('Wrong file format, cannot load.')
            return
        
        self.rotund_regions = fformat[1]
        ## Fixup for some changes to the CurveBoundaryRegion data structure
        for region in self.rotund_regions:
            if not hasattr( region, '_mesh' ): region._mesh = region.__dict__['mesh']
            if not hasattr( region, 'last_iterate_diff' ): region.last_iterate_diff = CurveBoundaryRegion([]).last_iterate_diff
            if not hasattr( region, 'iterate_threshold' ): region.iterate_threshold = CurveBoundaryRegion([]).iterate_threshold
        
        self.connections = fformat[2]
        
        self.SetImage( fformat[3] )
        
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
            #RotundAnnotateCanvas.kCrossSectionDragMode,
            #RotundAnnotateCanvas.kConnectionCurveDrawMode,
            #RotundAnnotateCanvas.kObliqueMode,
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
            ## See wx_annotate_rotund.py
            # from wx_annotate_rotund import ...
            #RotundAnnotateCanvas.kCrossSectionDragMode: CrossSectionDragMode,
            #RotundAnnotateCanvas.kConnectionCurveDrawMode: ConnectionCurveDrawMode,
            #RotundAnnotateCanvas.kObliqueMode: ObliqueMode,
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
        fibermesh.save_stroke( region_boundary, 'last_stroke.py' )
        
        curve_3D = zeros( ( len( region_boundary ), 3 ), dtype = float )
        curve_3D[:,:2] = region_boundary
        
        region = CurveBoundaryRegion( curve_3D )
        mesh, llconstraints = mesh_llconstraints_from_closed_curve_3D(
            curve_3D
            )
        
        region.mesh = mesh
        region.llconstraints = llconstraints
        
        self.annotations.rotund_regions.append( region )
        
        ## TODO: place in depth
        
        ## Refresh all
        self.GetTopLevelParent().Refresh()
    
    Color = {
        'interior': ( 200, 200, 255 ),
        'silhouette': ( 234, 145, 248 ),
        'sharp_cut': ( 174//2, 237//2, 111//2 ),
        'cross_section': ( 213, 99, 191 ),
        'spine': ( 74, 74, 74 )
        }
    
    def DrawCurveBoundaryRegion( self, region, what_to_draw = None ):
        if what_to_draw is None: what_to_draw = list(self.Color.keys())
        
        self.prepare3D()
        glDisable( GL_LIGHTING )
        
        if 'interior' in what_to_draw:
            ## Draw the interior
            ## TODO Q: Should I use draw it 2d, using gluTessPolygon (see my layers2d app for reference)
            if region.mesh is not None:
                glColor4ub( *( self.Color['interior'] + (50,) ) )
                ## draw_mesh_faces_smooth() is much faster than draw_mesh_faces_flat()
                draw_mesh_faces_smooth( region.mesh, True )
        
        ## Draw the boundary stroke.
        assert region.region_boundary is not None
        if region.region_boundary is not None:
            self.prepare2D()
            
            if 'silhouette' in what_to_draw:
                ## Draw silhouette edges
                glLineWidth( 8 )
                glColor4ub( *( self.Color['silhouette'] + (50,) ) )
                draw_lineloop( region.region_boundary )
    
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
                self.DrawCurveBoundaryRegion( region, what_to_draw = 'interior' )
            
            ## Now draw the decorating lines
            decs = set( self.Color.keys() )
            decs.remove( 'interior' )
            for region in self.annotations.rotund_regions:
                self.DrawCurveBoundaryRegion( region, what_to_draw = decs )
            
            ## Draw the connections either lightly or darkly, depending on whether we're in
            ## connecting mode
            #connection_alpha = 200 if type( self.mode ) is ConnectionCurveDrawMode else 50
            connection_alpha = 50
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
    
    def _AddPerfectRegion( self, curve_3D ):
        import fibermesh
        
        fibermesh.save_stroke( curve_3D[:,:2], 'last_stroke.py' )
        
        region = CurveBoundaryRegion( curve_3D )
        mesh, llconstraints = mesh_llconstraints_from_closed_curve_3D(
            curve_3D
            )
        
        region.mesh = mesh
        region.llconstraints = llconstraints
        
        self.annotations.rotund_regions.append( region )
        self.GetTopLevelParent().Refresh()
    
    def AddPerfectCircleRegion( self ):
        
        region_boundary = [ ( .5 * cos( theta ), .5 * sin( theta ) ) for theta in linspace( 0, 2*pi, 100 ) ]
        
        curve_3D = zeros( ( len( region_boundary ), 3 ), dtype = float )
        curve_3D[:,:2] = region_boundary
        self._AddPerfectRegion( curve_3D )
    
    def AddOvalRegion( self ):
        
        region_boundary = [ ( cos( theta ), .5 * sin( theta ) ) for theta in linspace( 0, 2*pi, 100 ) ]
        ## Scale by .9
        region_boundary = ( .9 * asarrayf( region_boundary ) ).tolist()
        
        curve_3D = zeros( ( len( region_boundary ), 3 ), dtype = float )
        curve_3D[:,:2] = region_boundary
        self._AddPerfectRegion( curve_3D )
    
    def AddRoundedCylinderRegion( self ):
        
        region_boundary = (
            [ ( .5 * cos( theta ) - .5, .5 * sin( theta ) ) for theta in linspace( pi/2, 3*pi/2, 100 ) ] +
            list(zip( linspace( -.5, .5, 10 ), -.5 * ones(10) ))[1:-1] +
            [ ( .5 * cos( theta ) + .5, .5 * sin( theta ) ) for theta in linspace( -pi/2, pi/2, 100 ) ] +
            list(zip( linspace( .5, -.5, 10 ), .5 * ones(10) ))[1:-1]
            )
        ## Scale by .9
        region_boundary = ( .9 * asarrayf( region_boundary ) ).tolist()
        
        curve_3D = zeros( ( len( region_boundary ), 3 ), dtype = float )
        curve_3D[:,:2] = region_boundary
        self._AddPerfectRegion( curve_3D )
    
    def AddSaddleRegion( self ):
        
        region_boundary = [ ( .5 * cos( theta ), .5 * sin( theta ) ) for theta in linspace( 0, 2*pi, 100 ) ]
        region_zs = [ 1*( x**2 - y**2 ) for x, y in region_boundary ]
        
        curve_3D = zeros( ( len( region_boundary ), 3 ), dtype = float )
        curve_3D[:,:2] = region_boundary
        curve_3D[:,2] = region_zs
        self._AddPerfectRegion( curve_3D )
    
    def AddGridRegion( self, N, M ):
        
        curve_3D, mesh, llconstraints = curve_3D_mesh_llconstraints_for_N_by_M_grid(
            N, M
            )
        
        import fibermesh
        fibermesh.save_stroke( curve_3D[:,:2], 'last_stroke.py' )
        
        region = CurveBoundaryRegion( curve_3D )
        region.mesh = mesh
        region.llconstraints = llconstraints
        
        self.annotations.rotund_regions.append( region )
        self.GetTopLevelParent().Refresh()

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
        self.show_wireframe = True
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
                draw_lineloop( region.region_boundary )
            
            if region.mesh is not None:
                if self.show_vertex_numbers:
                    ## Draw vertex number labels
                    glColor4f( 0,0,0, 1 )
                    for vi in range(len( region.mesh.vs )):
                        draw_text( region.mesh.vs[vi] + .01 * region.mesh.vertex_normals[vi], str(vi) )
                
                glEnable( GL_COLOR_MATERIAL )
                glColor3ub( *RotundAnnotateCanvas.Color['interior'] )
                if self.show_meshes:
                    ## Draw the faces.
                    glEnable( GL_LIGHTING )
                    draw_mesh_faces( region.mesh, True )
                    glDisable( GL_LIGHTING )
                if self.show_wireframe:
                    ## Draw the edges.
                    glColor3ub( 0, 0, 0 )
                    glLineWidth( 1. )
                    draw_mesh_edges( region.mesh )
            
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
        wx.Frame.__init__( self, parent, wx.ID_ANY, 'Test Environment' )
        
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
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show 3D Mesh Wireframe" )
        if self.result_canvas.show_wireframe: item.Check()
        self.Bind( wx.EVT_MENU, self.OnShowWireframe, item )
        
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
        
        item = testm.Append( wx.ID_ANY, "Draw a perfect saddle" )
        self.Bind( wx.EVT_MENU, lambda evt: self.annotate_canvas.AddSaddleRegion(), item )
        
        item = testm.Append( wx.ID_ANY, "Add a 3x3 grid" )
        self.Bind( wx.EVT_MENU, lambda evt: self.annotate_canvas.AddGridRegion( 3, 3 ), item )
        
        item = testm.Append( wx.ID_ANY, "Add a 5x5 grid" )
        self.Bind( wx.EVT_MENU, lambda evt: self.annotate_canvas.AddGridRegion( 5, 5 ), item )
        
        item = testm.Append( wx.ID_ANY, "Add a 7x7 grid" )
        self.Bind( wx.EVT_MENU, lambda evt: self.annotate_canvas.AddGridRegion( 7, 7 ), item )
        
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
    
    def OnShowWireframe( self, evt ):
        self.result_canvas.show_wireframe = not self.result_canvas.show_wireframe
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
