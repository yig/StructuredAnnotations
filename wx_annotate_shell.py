#!/usr/bin/env python3

from math import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *

import os



import wx
from wxBasicGLCanvas import BasicGLCanvas
from ui_actions import *
from wx_util import Timers, SaveWindowPositionAndSize, LoadWindowPositionAndSize, SaveDirectoryPath, LoadDirectoryPath
from mesh_gui import TriMesh, Camera, set_headlight, draw_mesh_faces_flat, draw_mesh_faces_smooth
from helpers import Struct, solve_for_fixed_points, smooth_open_stroke, smooth_closed_stroke, closest_point_index, contour_contains

from vecutil import *
from OpenGL.GL import *
from glhelpers import *
from GLTexture import GLTexture

## for session saving filenames we use timestamps
from datetime import datetime
import time

kImageFileExtensions = ( '.png', '.jpg', '.tiff', '.tif', '.gif', '.bmp' )

kDebugging = False
kPaperFigureMode = False
from vecutil import kEps
#from annotations_rotund import kDebugging, kEps

kUseTexture = True
kUseOblique = True

kSaveSessions = True

draw_mesh_faces = draw_mesh_faces_flat

def DeepCopyAnns( anns ):
    ##print '+++ begin deepcopy'
    from copy import deepcopy
    result = deepcopy( anns )
    ##print '--- end deepcopy'
    return result

class list_type_queryable( list ):
    def of_type( self, the_type ):
        return [ el for el in self if isinstance( el, the_type ) ]

class Annotations( object ):
    def __init__( self ):
        self.Reset()
    
    def Reset( self, keep_texture = False ):
        self.anns = list_type_queryable()
        self.anns_beautified = self.anns
        self.oblique = [0.,0.]
        if not keep_texture:
            self.image_path = None
            ## Shared OpenGL contexts just aren't working anymore, so I need to make
            ## a per-context dictionary.
            #self.gltex = None
            self._gltex = {}
    
    def Beautify( self ):
        pass
    
    def MoveToFront( self, ann ):
        assert ann in self.anns
        self.anns.remove( ann )
        self.anns.insert( 0, ann )
    
    def AnnotationsFromFile( filename ):
        result = Annotations()
        result.Load( filename )
        return result
    AnnotationsFromFile = staticmethod( AnnotationsFromFile )
    
    def SetImage( self, image_path ):
        #from GLTexture import GLTexture
        from PIL import Image
        try:
            ## We can't create the GLTexture yet, because we
            ## may not have an OpenGL context yet.
            #self.gltex = GLTexture( image_path )
            self._gltex = {}
            self.img = Image.open( image_path )
            self.image_path = image_path
        except IOError:
            print('Failed to load texture:', image_path)
            return
        
        ## Set the rectangle's size according to the texture's aspect ratio.
        self.tex_size = arrayf( self.img.size )
        self.tex_size *= 2./self.tex_size.max()
        
        ## Center the image.
        self.tex_pos = zeros(2) - .5*self.tex_size
    
    def gltex( self, current_context ):
        if self.image_path is None: return None
        
        if current_context not in self._gltex:
            self._gltex[ current_context ] = GLTexture( self.image_path )
        
        return self._gltex[ current_context ]
    
    
    FFORMAT_V1 = 'wx_annotate_shell file format v1'
    FFORMAT_V2 = 'wx_annotate_shell file format v2'
    
    def SaveInBackground( self, *args, **kwargs ):
        from threading import Thread
        copy = Annotations()
        ## NOTE: If ever the fields Save() uses change, I need to update these lines.
        copy.oblique = list( self.oblique )
        copy.image_path = self.image_path
        copy.anns = DeepCopyAnns( self.anns )
        thread = Thread( target = lambda: copy.Save( *args, **kwargs ) )
        thread.start()
    
    def Save( self, filename, overwrite = False ):
        if not overwrite and os.path.exists( filename ):
            print('File exists, not overwriting:', filename)
            return
        
        if overwrite and os.path.isdir( filename ):
            print('Directory with same name as desired file already exists, can\'t overwrite:', filename)
            return
        
        ## If we made it here, we're fine.
        if overwrite and os.path.exists( filename ):
            print('Saving file (overwriting):', filename)
        else:
            print('Saving file:', filename)
        
        ## This is our file format:
        fformat = ( self.FFORMAT_V2, self.anns, self.image_path, self.oblique )
        ## NOTE: If ever these fields change, I need to update SaveInBackground().
        
        import pickle
        f = open( filename, 'wb' )
        pickle.dump( fformat, f, pickle.HIGHEST_PROTOCOL )
        print('...Saved')
    
    def Load( self, filename ):
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
        
        self.anns = fformat[1]
        
        self.SetImage( fformat[2] )
        
        self.oblique = fformat[3]
        
        self.Beautify()


class AnnotateCanvasBase( BasicGLCanvas ):
    
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
            #self.annotations.SetImage( '/mixed media/pictures/ColossalKatamari.png' )
            #self.annotations.SetImage( os.path.join( 'scenes', 'cactus.png' ) )
            #self.annotations.SetImage( os.path.join( 'scenes', 'artwork', 'crow bakery - rabbit.png' ) )
            #self.annotations.SetImage( os.path.join( 'scenes', 'artwork', 'bugs bunny', 'LT-Valentine-Bugs-Bunny-2.jpg' ) )
            self.annotations.SetImage( os.path.join( 'scenes', 'artwork', 'xy grid2.png' ) )
        
        self.debugging = False
    
    def GetBestSize( self ):
        return (300,300)
    
    def prepare2D( self ):
        glDisable( GL_DEPTH_TEST )
        glDisable( GL_LIGHTING )
        glDisable( GL_COLOR_MATERIAL )
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        ## Center the square inside the viewport.
        glOrtho(*( list( self.world_bounds ) + [ -100, 100 ] ))
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
    
    def prepare3D( self ):
        glEnable( GL_DEPTH_TEST )
        glEnable( GL_LIGHTING )
        glEnable( GL_COLOR_MATERIAL )
        
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


class AnnotateCanvas( AnnotateCanvasBase ):
    
    def __init__( self, *args, **kwargs ):
        AnnotateCanvasBase.__init__( self, *args, **kwargs )
        
        self.ModeID_to_ModeFactory = {}
        self.ModeID_to_LongName = {}
        self.ModeID_to_ShortName = {}
        #self.mode = None
        self.mode = ModeProxy( self )
        #self.SetMode( AnnotationCreateMode( self ) )
    
    def SetMode( self, mode ):
        assert type( mode ) in list(self.ModeID_to_ModeFactory.values())
        
        ## Don't set the mode if we're already in it.
        ## TODO Q: Should I turn this on?  I've never had it on before.  It seems useful
        ##         when debugging to reset all state.
        if type( mode ) is type( self.mode ): return
        
        ## Oops: I can't use forward_events because it's broken and only works to one level.
        ## My modes need to be able to forward events, too.
        #if self.mode is not None: unforward_events( self.mode )
        self.mode = mode
        #forward_events( self, self.mode )
        
        ## Modes have different draw functions, so refresh now that we have a new one.
        self.Refresh()
        
        ## Some modes (or event handlers used by modes) may disallow undo/redo.
        ## Always re-allow it when switching modes.
        self.EditMenuUndoRedoAllow()
    
    def EditMenuUndoRedoAllow( self ):
        '''
        This method will be overridden by the frame, who will set a different method to this
        property.
        '''
        pass
    
    def RegisterMode( self, mode_proxy_factory, long_name, short_name ):
        '''
        Returns the mode id
        '''
        
        mode_id = 1 if len( self.ModeID_to_LongName ) == 0 else max( self.ModeID_to_LongName.keys() ) + 1
        
        assert mode_id not in self.ModeID_to_ModeFactory
        assert mode_id not in self.ModeID_to_LongName
        assert mode_id not in self.ModeID_to_ShortName
        
        self.ModeID_to_ModeFactory[ mode_id ] = mode_proxy_factory
        self.ModeID_to_LongName[ mode_id ] = long_name
        self.ModeID_to_ShortName[ mode_id ] = short_name
        
        return mode_id
    
    def GetAllModeIDs( self ):
        return list(self.ModeID_to_ModeFactory.keys())
    
    def GetLongNameForModeID( self, ID ):
        return self.ModeID_to_LongName[ ID ]
    
    def GetShortNameForModeID( self, ID ):
        return self.ModeID_to_ShortName[ ID ]
    
    def GetModeIDForShortName( self, short_name ):
        for mode_id, name in self.ModeID_to_ShortName.items():
            if name == short_name: return mode_id
        raise RuntimeError("No mode with given short name:").with_traceback(short_name)
    
    def SetModeByID( self, ID ):
        mode_type = self.ModeID_to_ModeFactory[ ID ]
        if kDebugging or ID != self.GetCurrentModeID():
            print('Switching to mode:', self.GetLongNameForModeID( ID ))
        self.SetMode( mode_type( self ) )
    
    def GetCurrentModeID( self ):
        for id, mode_factory in self.ModeID_to_ModeFactory.items():
            if type( self.mode ) == mode_factory:
                return id
        
        return None
    
    def PickAnnotationsByType( self, screen_pos, pick_type ):
        ## TODO Q: Do I want this method or can I call this directly on the anns[foo] object?
        ##         Do I need screen2world()?
        pass
    
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
    
    def DrawScene( self, **kwargs ):
        glClearColor( .9,.9,.9,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        ## Draw the gltex if there is one.
        if self.annotations.gltex( self.context ) is not None:
            self.prepare2D()
            
            glColor4f( 1,1,1, .5 )
            draw_gltexture( self.annotations.gltex( self.context ), self.annotations.tex_pos, self.annotations.tex_size )
        
        if 'callback' in kwargs: kwargs[callback]()
        
        self.SwapBuffers()

class ResultCanvas( BasicGLCanvas ):
    
    EYE_DISTANCE = 100.
    
    def __init__( self, parent, annotations, pos = wx.DefaultPosition, size = wx.DefaultSize, shared = None ):
        BasicGLCanvas.__init__( self, parent, pos, size, debugging = True, shared = shared )
        
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
        self.up_is_y = True
        
        debugging = False
        self.rotator = RotateHandler( self.camera, self.Refresh, self.event2world, debugging = debugging )
        forward_events( self, self.rotator )
    
    def GetBestSize( self ):
        return (300,300)
    
    def doDraw( self ):
        glClearColor( .5,.5,.5,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        self.prepare3D()
        
        #self.drawTexture()
        ## Clear the depth buffer.
        ## Actually, this looks weird.
        # glClear( GL_DEPTH_BUFFER_BIT )
        
        #for anns in self.annotations.anns.values(): anns.DrawAnnotations( self.prepare2D, self.prepare3D, None )
        #for anns in self.annotations.anns.values(): anns.DrawResults( self.prepare2D, self.prepare3D, self )
        
        ## Draw the gltex if there is one.
        if self.annotations.gltex( self.context ) is not None:
            self.drawTexture(
                self.annotations.gltex( self.context ),
                self.annotations.tex_pos,
                self.annotations.tex_size
                )
        
        self.SwapBuffers()
    
    def drawTexture( self, gltex, tex_pos, tex_size ):
        if self.image_depth is None: return
        if self.show_background_image is False: return
        
        cull_face = glIsEnabled( GL_CULL_FACE )
        if cull_face: glDisable( GL_CULL_FACE )
        
        glDisable( GL_LIGHTING )
        
        glMatrixMode( GL_MODELVIEW )
        glPushMatrix()
        glTranslate( 0,0,self.image_depth )
        
        glColor4f( 1,1,1, .5 )
        draw_gltexture( gltex, tex_pos, tex_size )
        
        glPopMatrix()
        
        if cull_face: glEnable( GL_CULL_FACE )
    
    def prepare2D( self ):
        glDisable( GL_DEPTH_TEST )
        glDisable( GL_LIGHTING )
        glDisable( GL_COLOR_MATERIAL )
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        #glOrtho(-1, 1, -1, 1, -100, 100)
        glOrtho( *( list( self.world_bounds ) + [ -100, 100 ] ) )
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
    
    def prepare3D( self ):
        glEnable( GL_DEPTH_TEST )
        glEnable( GL_LIGHTING )
        glEnable( GL_COLOR_MATERIAL )
        
        ## Keep the up vector pointed in positive y, so we don't get twisted.
        if self.up_is_y: self.camera.up[:] = (0,1,0)
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
    
    def __init__( self, parent, annotations, annotate_canvas_factory = AnnotateCanvas, result_canvas_factory = ResultCanvas ):
        wx.Frame.__init__( self, parent, wx.ID_ANY, 'Annotation Modeler', style = wx.DEFAULT_FRAME_STYLE )
        # wx.ConfigBase.Get().SetAppName('wx_annotate_shell')
        wx.App.Get().SetAppName( 'com.yotamgingold.wx_annotate_shell' )
        wx.App.Get().SetAppDisplayName( 'Structured Annotations for 2D-to-3D Modeling' )
        
        self.annotations = annotations
        
        self.undo_stack = []
        self.redo_stack = []
        
        ### Add a toolbar for the modes.  Do this before adding the canvases.
        ## Vertical looks more like Illustrator, but (a) I can't get items to
        ## visually deselect and (b) when the window isn't tall enough, there is no
        ## indication that there are more tools.
        #TBFLAGS = wx.TB_VERTICAL | wx.TB_FLAT | wx.TB_TEXT
        ## Horizontal toolbar looks more OSX.
        TBFLAGS = wx.TB_HORIZONTAL | wx.TB_FLAT | wx.TB_TEXT
        tb = self.CreateToolBar( TBFLAGS )
        #tbsize = ( 64, 64 )
        #tbsize = ( 48, 48 )
        ## For the video I want smaller buttons
        tbsize = ( 32, 32 )
        tb.SetToolBitmapSize( tbsize )
        
        panel = wx.Panel( self, wx.ID_ANY, size=(700,700) )
        self.sizer = wx.BoxSizer( wx.VERTICAL )
        
        ## Add the side-by-side gl views
        glviews = wx.BoxSizer( wx.HORIZONTAL )
        
        
        ## Create the AnnotateCanvas.
        #canvas = annotate_canvas_factory( panel, self.annotations, shared = self.result_canvas.context )
        canvas = annotate_canvas_factory( panel, self.annotations )
        canvas.Show()
        self.annotate_canvas = canvas
        canvas.SetSize( canvas.GetBestSize() )
        glviews.Add( canvas, 1, wx.EXPAND )
        
        
        ## Create the ResultCanvas.
        ## Creating a wx.MiniFrame makes it not require an extra click when operating the two views.
        ## Unfortunately, it also makes the window float (turning off STAY_ON_TOP and
        ## FRAME_FLOAT_ON_PARENT has no effect).
        #results_window = wx.Frame( self, wx.ID_ANY, 'Results',
        results_window = wx.MiniFrame( self, wx.ID_ANY, 'Results',
            style = wx.DEFAULT_FRAME_STYLE & ~wx.CLOSE_BOX
            #style = wx.DEFAULT_FRAME_STYLE & ~wx.CLOSE_BOX
            
            #style = wx.DEFAULT_FRAME_STYLE & ~( wx.CLOSE_BOX | wx.STAY_ON_TOP | wx.FRAME_FLOAT_ON_PARENT )
            ## Neither STAY_ON_TOP, FRAME_FLOAT_ON_PARENT, nor FRAME_TOOL_WINDOW make clicks
            ## pass through to the window.
            #style = ( wx.DEFAULT_FRAME_STYLE | wx.STAY_ON_TOP ) & ~( wx.CLOSE_BOX )
            #style = ( wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT ) & ~( wx.CLOSE_BOX )
            ## This generates a window with no frame, and it still doesn't float on top.
            #style = wx.FRAME_TOOL_WINDOW
            )
        # canvas = result_canvas_factory( results_window, self.annotations, shared = self.annotate_canvas.context )
        canvas = result_canvas_factory( results_window, self.annotations )
        self.result_canvas = canvas
        canvas.SetSize( canvas.GetBestSize() )
        #glviews.Add( canvas, 1, wx.EXPAND )
        results_window.SetSize( (300,300) )
        results_window.Show()
        
        
        self.sizer.Add( glviews, 1, wx.EXPAND )
        
        
        ## Users switch between the annotate and rotate views a lot, so let's let the focusing click
        ## pass through. (Frames have .SetFocus() and .Raise() methods.)
        #results_window.Bind( wx.EVT_SET_FOCUS, lambda evt: evt.Skip() )
        #self.Bind( wx.EVT_SET_FOCUS, lambda evt: evt.Skip() )
        #results_window.Bind( wx.EVT_SET_FOCUS, lambda evt: evt.Skip() )
        
        ## I don't really use these buttons anymore, I use the menus.
        '''
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
        '''
        
        panel.SetSizer( self.sizer )
        ## Don't call Fit: https://wiki.wxpython.org/UsingSizers
        # self.sizer.Fit( panel )
        # self.Fit()
        
        ## Set the last working directory to the current directory if no preference file is saved.
        try:
            self.last_wd = LoadDirectoryPath( '*' )
        except ValueError:
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
        
        item = filem.Append( wx.ID_ANY, "Save Annotation Screenshot...\tctrl-[" ) #\tctrl-3" )
        self.Bind( wx.EVT_MENU, lambda evt: self.SaveScreenShot( self.annotate_canvas, 'annotate canvas' ), item )
        
        item = filem.Append( wx.ID_ANY, "Save Result Screenshot...\tctrl-]" ) #\tctrl-4" )
        self.Bind( wx.EVT_MENU, lambda evt: self.SaveScreenShot( self.result_canvas, 'result canvas' ), item )
        
        item = filem.Append( wx.ID_ANY, "Export to OBJ..." )
        self.Bind( wx.EVT_MENU, lambda evt: self.ExportToOBJ( self.annotations.anns_beautified ), item )
        
        ## This should be ignored on OSX:
        filem.AppendSeparator()
        item = filem.Append( wx.ID_EXIT, "&Exit" )
        
        menubar.Append( filem, "&File" )
        
        ##
        
        editm = wx.Menu()
        
        undo_item = editm.Append( wx.ID_UNDO, "Undo\tctrl-z" )
        redo_item = editm.Append( wx.ID_REDO, "Redo\tctrl-shift-z" )
        
        self.Bind( wx.EVT_MENU, self.OnEditMenuUndo, undo_item )
        self.Bind( wx.EVT_MENU, self.OnEditMenuRedo, redo_item )
        
        ## The undo and redo stacks are initially empty.
        undo_item.Enable( False )
        redo_item.Enable( False )
        ## Undo and redo are allowed (were there something to undo or redo).
        self.undo_redo_allowed = True
        
        self.undo_item = undo_item
        self.redo_item = redo_item
        self.annotate_canvas.SaveStateForUndo = self.SaveStateForUndo
        self.annotate_canvas.EditMenuUndoRedoAllow = self.EditMenuUndoRedoAllow
        self.annotate_canvas.EditMenuUndoRedoDisallow = self.EditMenuUndoRedoDisallow
        self.annotate_canvas.SwitchToolbarTools = self.SwitchToolbarTools
        
        editm.AppendSeparator()
        
        item = editm.Append( wx.ID_ANY, "Unlock All Shapes" )
        self.Bind( wx.EVT_MENU, self.OnUnlockAllShapes, item )
        
        editm.AppendSeparator()
        
        item = editm.Append( wx.ID_ANY, "Reset" )
        self.Bind( wx.EVT_MENU, lambda evt: self.OnMenuReset(), item )
        
        menubar.Append( editm, "&Edit" )
        
        ##
        mode_ids = self.annotate_canvas.GetAllModeIDs()
        
        modem = wx.Menu()
        for mode_id in mode_ids:
            ## TODO: make these into a check item group and have only one checked.
            accel = ""
            if mode_id < 10:
                accel = "\tctrl-%s" % (mode_id,)
            item = modem.Append( mode_id, self.annotate_canvas.GetLongNameForModeID( mode_id ) + accel )
            
            self.Bind( wx.EVT_MENU, lambda evt: self.SwitchToolbarTools( evt.Id ), item )
        menubar.Append( modem, "&Tools" )
        
        ##
        
        viewm = wx.Menu()
        
        item = viewm.Append( wx.ID_ANY, "Reset 2D Annotation View" )
        self.Bind( wx.EVT_MENU, self.OnReset2DAnnotationView, item )
        
        item = viewm.Append( wx.ID_ANY, "Reset 3D Result View\tctrl-r" )
        self.Bind( wx.EVT_MENU, self.OnReset3DResultView, item )
        
        item = viewm.Append( wx.ID_ANY, "Set 3D Result View...\tctrl-shift-d" )
        self.Bind( wx.EVT_MENU, self.OnSet3DResultView, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Auto-Spin 3D View\tctrl-shift-r" )
        self.spinner_id = None
        self.Bind( wx.EVT_MENU, self.OnAutoSpin3DResultView, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Result View Image Behind\tctrl-b" )
        if self.result_canvas.image_depth < 0: item.Check()
        self.Bind( wx.EVT_MENU, self.OnToggleResultViewImageDepth, item )
        
        viewm.AppendSeparator()
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show Annotations\tctrl-a" )
        if self.annotate_canvas.show_annotations: item.Check()
        self.Bind( wx.EVT_MENU, self.OnShowAnnotations, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Paper Figure Mode" )
        if kPaperFigureMode: item.Check()
        self.Bind( wx.EVT_MENU, self.OnPaperFigureMode, item )
        
        try:
            import annotations_sweep
            item = viewm.AppendCheckItem( wx.ID_ANY, "Pinch-Free Shapes" )
            if annotations_sweep.kNoPinch: item.Check()
            self.Bind( wx.EVT_MENU, self.OnNoPinch, item )
        except ImportError: pass
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show Background Image" )
        if self.result_canvas.show_background_image: item.Check()
        self.Bind( wx.EVT_MENU, self.OnShowBackgroundImage, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Maintain Up Vector" )
        if self.result_canvas.up_is_y: item.Check()
        self.Bind( wx.EVT_MENU, self.OnUpIsY, item )
        
        '''
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show Vertex Numbers\tctrl-n" )
        if self.result_canvas.show_vertex_numbers: item.Check()
        self.Bind( wx.EVT_MENU, self.OnShowVertexNumbers, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Smooth Shading" )
        if draw_mesh_faces == draw_mesh_faces_smooth: item.Check()
        self.Bind( wx.EVT_MENU, self.OnSmoothShading, item )
        '''
        
        if hasattr( self.result_canvas, "npr" ):
            item = viewm.AppendCheckItem( wx.ID_ANY, "Cartoon Shading" )
            if self.result_canvas.npr: item.Check()
            self.Bind( wx.EVT_MENU, self.OnCartoonShading, item )
        
        if hasattr( self.result_canvas, "silhouettes" ):
            item = viewm.AppendCheckItem( wx.ID_ANY, "Silhoutte Rendering" )
            if self.result_canvas.silhouettes: item.Check()
            self.Bind( wx.EVT_MENU, self.OnSilhouetteRendering, item )
        
        if hasattr( self.annotate_canvas, "shape_fill" ):
            item = viewm.AppendCheckItem( wx.ID_ANY, "Filled Shapes" )
            if self.annotate_canvas.shape_fill: item.Check()
            self.Bind( wx.EVT_MENU, self.OnShapeFill, item )
        
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show 3D Meshes" )
        if self.result_canvas.show_meshes: item.Check()
        self.Bind( wx.EVT_MENU, self.OnShowMeshes, item )
        
        try:
            import sweep_surface
            item = viewm.AppendCheckItem( wx.ID_ANY, "Low-Resolution 3D Meshes" )
            if sweep_surface.kLowRes: item.Check()
            self.Bind( wx.EVT_MENU, self.OnLowRes, item )
        except ImportError: pass
        
        '''
        item = viewm.AppendCheckItem( wx.ID_ANY, "Use Oblique Camera" )
        #if self.use_oblique: item.Check()
        if kUseOblique: item.Check()
        self.Bind( wx.EVT_MENU, self.OnUseOblique, item )
        '''
        
        item = viewm.Append( wx.ID_ANY, "Refresh screen" )
        self.Bind( wx.EVT_MENU, lambda evt: self.GetTopLevelParent().Refresh(), item )
        
        if kDebugging:
            item = viewm.AppendCheckItem( wx.ID_ANY, "Show debugging annotations" )
            item.Check()
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
        
        if kDebugging: menubar.Append( testm, "&Test" )
        
        ##
        
        self.SetMenuBar( menubar )
        
        
        ## Add modes to the toolbar
        ## TODO: make these into a check item group and have only one checked.
        ## UPDATE: Create a wx.LogNull() so we don't get color profile warnings.
        donterror = wx.LogNull()
        for mode_id in mode_ids:
            short_name = self.annotate_canvas.GetShortNameForModeID( mode_id )
            long_name = self.annotate_canvas.GetLongNameForModeID( mode_id )
            
            img = wx.Image( os.path.join( ".", "modes", short_name.lower() + ".png" ), wx.BITMAP_TYPE_ANY )
            img = img.Scale( tbsize[0], tbsize[1], quality = wx.IMAGE_QUALITY_HIGH )
            bitmap = img.ConvertToBitmap()
            # tb.AddRadioLabelTool( mode_id, short_name, bitmap, shortHelp=short_name, longHelp=long_name )
            tb.AddRadioTool( mode_id, short_name, bitmap, shortHelp=short_name, longHelp=long_name )
            #tb.AddTool( mode_id, bitmap, isToggle = True )
            #self.Bind( wx.EVT_TOOL, lambda evt: self.annotate_canvas.SetModeByID( evt.Id ), id=mode_id )
            self.Bind( wx.EVT_TOOL, lambda evt: self.SwitchToolbarTools( evt.Id ), id = mode_id )
        
        del donterror
        
        tb.Realize()
        self.toolbar = tb
        ## The default mode is set by the annotate canvas.
        self.SwitchToolbarTools( self.annotate_canvas.GetCurrentModeID() )
        
        ## Q: Do I need to call this again?
        ## A: No. See: https://wiki.wxpython.org/UsingSizers
        # self.Fit()
        ## The vertical toolbar isn't being taken into account... :-p
        ## TODO Q: Why is the vertical height off by 22?  I have no idea.
        self.SetSize( array( self.GetSize() ) + (tbsize[0],-22) )
        
        ## Move the results window to the right-hand-side of the annotation window.
        results_window.SetPosition( array( self.GetPosition() ) + (self.GetSize()[0],0) )
        
        
        ## Save window positions when we quit.
        def window_saver():
            #print 'window_saver()'
            
            SaveWindowPositionAndSize( self, 'annotate' )
            SaveWindowPositionAndSize( results_window, 'result' )
            ## I save this every time it is updated, so I don't need to save it here (on exit).
            #SaveDirectoryPath( self.last_wd, '*' )
            if len( self.annotations.anns ) > 0: self.annotations.Save( 'last_session.annscene', overwrite = True )
            self.SaveSession( save_in_background = False )
        
        ## We save windows upon receiving the close event, and then "Skip" the event
        ## to tell wxPython that the default close event handling should proceed.
        self.Bind( wx.EVT_CLOSE, lambda evt: ( window_saver(), evt.Skip() ) )
        self.Bind( wx.EVT_MENU, lambda evt: ( window_saver(), evt.Skip() ), id = wx.ID_EXIT )
        results_window.Bind( wx.EVT_MENU, lambda evt: ( window_saver(), evt.Skip() ), id = wx.ID_EXIT )
        
        ## TODO Q: Veto the close event so no one can close the results window?
        #results_window.Bind( wx.EVT_CLOSE, lambda evt: ( window_saver(), evt.Veto() ) )
        
        
        try:
            ann_pos, ann_size = LoadWindowPositionAndSize( 'annotate' )
            res_pos, res_size = LoadWindowPositionAndSize( 'result' )
            
            self.SetSize( ann_size )
            self.SetPosition( ann_pos )
            
            results_window.SetSize( res_size )
            results_window.SetPosition( res_pos )
        except ValueError:
            ## Calling int() on an empty string raises a ValueError.
            pass
            #raise
    
    def SwitchToolbarTools( self, tool_id ):
        if tool_id is None: return
        
        ## TODO Q: Why doesn't this affect which tool is highlighted?  Is RadioTool the wrong kind?
        #tb.SetToggle( evt.Id, True )
        ## A1: From http://aspn.activestate.com/ASPN/Mail/Message/wxpython-users/631657 :
        ##     (Still doesn't work.)
        for id in self.annotate_canvas.GetAllModeIDs(): self.toolbar.ToggleTool( id, False )
        self.toolbar.ToggleTool( tool_id, True )
        self.toolbar.Refresh()
        #self.Refresh()
        self.annotate_canvas.SetModeByID( tool_id )
    
    def Refresh( self, *args, **kwargs ):
        
        # print( "Refresh called" )
        
        ## Because the results canvas is a separate window, we must override the frame's refresh
        ## and refresh it too.
        self.result_canvas.Refresh()
        
        ## Call through to the superclass for the regular behavior.
        super().Refresh( *args, **kwargs )
        
        ## The call to wx.Frame.Refresh() doesn't seem to call self.annotate_canvas.Refresh().
        ## Let's refresh it manually.
        self.annotate_canvas.Refresh()
    
    def EditMenuUndoRedoDisallow( self ):
        self.undo_redo_allowed = False
        
        self.undo_item.Enable( False )
        self.redo_item.Enable( False )
    
    def EditMenuUndoRedoAllow( self ):
        self.undo_redo_allowed = True
        
        self.undo_item.Enable( len( self.undo_stack ) > 0 )
        self.redo_item.Enable( len( self.redo_stack ) > 0 )
    
    def OnEditMenuUndo( self, evt ):
        assert len( self.undo_stack ) > 0
        
        #import pickle
        #self.redo_stack.append( pickle.dumps( self.annotations.anns, pickle.HIGHEST_PROTOCOL ) )
        #self.annotations.anns = pickle.loads( self.undo_stack[-1] )
        
        self.redo_stack.append( DeepCopyAnns( ( self.annotations.anns, self.annotate_canvas.locked_shapes ) ) )
        self.annotations.anns, self.annotate_canvas.locked_shapes = self.undo_stack.pop()
        
        self.annotations.Beautify()
        
        if self.undo_redo_allowed:
            self.undo_item.Enable( len( self.undo_stack ) > 0 )
            self.redo_item.Enable( len( self.redo_stack ) > 0 )
        
        self.GetTopLevelParent().Refresh()
    
    def OnEditMenuRedo( self, evt ):
        assert len( self.redo_stack ) > 0
        
        #import pickle
        #self.undo_stack.append( pickle.dumps( self.annotations.anns, pickle.HIGHEST_PROTOCOL ) )
        #self.annotations.anns = pickle.loads( self.redo_stack[-1] )
        
        self.undo_stack.append( DeepCopyAnns( ( self.annotations.anns, self.annotate_canvas.locked_shapes ) ) )
        self.annotations.anns, self.annotate_canvas.locked_shapes = self.redo_stack.pop()
        
        self.annotations.Beautify()
        
        if self.undo_redo_allowed:
            self.undo_item.Enable( len( self.undo_stack ) > 0 )
            self.redo_item.Enable( len( self.redo_stack ) > 0 )
        
        self.GetTopLevelParent().Refresh()
    
    def SaveStateForUndo( self ):
        '''
        Captures the current state for restoration in case of an undo.
        '''
        
        self.SaveSession()
        
        ## Q: Should I use deepcopy or pickle?  This gets slow when objects with big TriMesh's are
        ##    around.
        ## A: cPickle is faster than a simple copy.deepcopy,
        ##    but my custom TriMesh.__deepcopy__ is fastest.
        
        #import pickle
        #self.undo_stack.append( pickle.dumps( self.annotations.anns, pickle.HIGHEST_PROTOCOL ) )
        
        self.undo_stack.append( DeepCopyAnns( ( self.annotations.anns, self.annotate_canvas.locked_shapes ) ) )
        
        ## Make sure the undo stack never grows past 10 items.
        ## I'm not storing meshes in these objects, so they should be quite tiny.
        del self.undo_stack[:-10]
        ## Saving state nukes the redo stack.
        del self.redo_stack[:]
        
        if self.undo_redo_allowed:
            self.undo_item.Enable( len( self.undo_stack ) > 0 )
            self.redo_item.Enable( len( self.redo_stack ) > 0 )
    
    def SaveSession( self, save_in_background = True ):
        '''
        Save the state to a time-stamped file in a sessions directory.
        '''
        
        if not kSaveSessions: return
        
        sessions_dir = "sessions"
        
        ## Create the sessions directory if it doesn't exist.
        if not os.path.exists( sessions_dir ): os.mkdir( sessions_dir )
        
        ## Check for errors.
        ## At this point we have already attempted to create the sessions directory,
        ## so if it doesn't exist or is not a directory, we can't continue.
        if not os.path.isdir( sessions_dir ):
            ## Something has gone wrong and we don't have a sessions directory.
            return
        
        if len( self.annotations.anns ) > 0:
            # SimpleDateFormat formatter = new SimpleDateFormat( "yyyy-MM-dd HH:mm:ss.SSS z" );
            now = datetime.today()
            basename = now.strftime("%Y-%m-%d %H:%M:%S") + ( ".%d %s" % ( now.microsecond / 1000, time.strftime( "%Z" ) ) )
            
            filename = os.path.join( sessions_dir, basename + ".annscene" )
            
            ## annotations.Save() already prints out this string.
            ## print "[Saving session to: " + filename + "]"
            
            ## Q: Is this slow?
            ## A: No, not at all.
            #print 'start saving'
            if save_in_background:
                self.annotations.SaveInBackground( filename, overwrite = True )
            else:
                self.annotations.Save( filename, overwrite = True )
            #print 'done saving'
    
    def SaveScreenShot( self, canvas, suggested_name ):
        dialog = wx.FileDialog(
            self, "Save screenshot...", self.last_wd,
            suggested_name + '.png',
            "PNG|*.png|JPEG|*.jpg",
            wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
            )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            SaveDirectoryPath( self.last_wd, '*' )
            
            canvas.SaveToFile( str( dialog.GetPath() ) )
    
    def ExportToOBJ( self, anns ):
        suggested_name = 'annotated scene'
        if self.save_filename is not None:
            suggested_name = os.path.split( os.path.splitext( self.save_filename )[0] )[1]
        elif self.annotations.image_path: suggested_name = os.path.split( os.path.splitext( self.annotations.image_path )[0] )[1]
        
        dialog = wx.FileDialog(
            self, "Export to OBJ...", self.last_wd,
            suggested_name + '.obj',
            "OBJ|*.obj",
            wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
            )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            SaveDirectoryPath( self.last_wd, '*' )
            
            from mesh_gui import TriMesh
            all_meshes = TriMesh.FromTriMeshes( [
                shape.mesh for shape in anns if hasattr( shape, 'mesh' )
                ] )
            all_meshes.write_OBJ( str( dialog.GetPath() ) )
    
    def OnToggleResultViewImageDepth( self, evt ):
        self.result_canvas.image_depth *= -1.
        self.result_canvas.Refresh()
    
    #def OnToggleTestNewOverlap( self, evt ):
    #    global kTestNewOverlap
    #    kTestNewOverlap = not kTestNewOverlap
    
    def OnUnlockAllShapes( self, evt ):
        self.SaveStateForUndo()
        self.annotate_canvas.locked_shapes.clear()
        self.annotate_canvas.Refresh()
    
    def OnReset2DAnnotationView( self, evt ):
        self.annotate_canvas.ResetZoom()
        self.annotate_canvas.Refresh()
    
    def OnReset3DResultView( self, evt ):
        self.result_canvas.reset_rotation()
        self.result_canvas.Refresh()
    
    def OnAutoSpin3DResultView( self, evt ):
        if self.spinner_id is not None:
            Timers.CancelTimer( self.spinner_id )
            self.spinner_id = None
        else:
            def spin_percent( percent ):
                ## Silently fail if we're still being called-back when the window has closed.
                try:
                    self.result_canvas.camera.Rotate2( .5*percent, 0. )
                    self.result_canvas.Refresh()
                except: pass
            
            ## period in seconds
            period = 9.
            ## frames per second
            ## NOTE: We only screen record at 15 fps, so don't refresh more than that.
            fps = 15.
            self.spinner_id = Timers.RequestTimer(
                lambda: spin_percent( 2./(fps*period) ),
                ## One degree 15 times a second.
                1000./fps,
                once = False
                )
    
    def OnSet3DResultView( self, evt ):
        dlg = wx.TextEntryDialog( self, "", "Eye Position; Up Vector:",
            '; '.join( [ ' '.join( map( repr, self.result_canvas.camera.eye ) ), ' '.join( map( repr, self.result_canvas.camera.up ) ) ] )
            )
        if dlg.ShowModal() == wx.ID_OK:
            try:
                eye, up = dlg.GetValue().split(';')
                eye = list(map( float, eye.strip().split(' ') ))
                up = list(map( float, up.strip().split(' ') ))
                
                self.result_canvas.camera.eye[:] = eye
                self.result_canvas.camera.up[:] = up
                self.result_canvas.Refresh()
                
            except ( TypeError, ValueError ):
                print('Value was not a number!')
        dlg.Destroy()
    
    def OnShowAnnotations( self, evt ):
        self.annotate_canvas.show_annotations = not self.annotate_canvas.show_annotations
        self.annotate_canvas.Refresh()
    
    def OnPaperFigureMode( self, evt ):
        global kPaperFigureMode
        kPaperFigureMode = not kPaperFigureMode
        self.annotate_canvas.Refresh()
    
    def OnShowBackgroundImage( self, evt ):
        self.result_canvas.show_background_image = not self.result_canvas.show_background_image
        self.result_canvas.Refresh()
    
    def OnUpIsY( self, evt ):
        self.result_canvas.up_is_y = not self.result_canvas.up_is_y
        self.result_canvas.Refresh()
    
    def OnLowRes( self, evt ):
        try:
            import sweep_surface
            sweep_surface.set_lowres( not sweep_surface.kLowRes )
        except ImportError: pass
    
    def OnShowVertexNumbers( self, evt ):
        self.result_canvas.show_vertex_numbers = not self.result_canvas.show_vertex_numbers
        self.result_canvas.Refresh()
    
    def OnSmoothShading( self, evt ):
        global draw_mesh_faces
        if draw_mesh_faces == draw_mesh_faces_smooth: draw_mesh_faces = draw_mesh_faces_flat
        else: draw_mesh_faces = draw_mesh_faces_smooth
        self.GetTopLevelParent().Refresh()
    
    def OnCartoonShading( self, evt ):
        self.result_canvas.npr = not self.result_canvas.npr
        self.result_canvas.Refresh()
    
    def OnSilhouetteRendering( self, evt ):
        self.result_canvas.silhouettes = not self.result_canvas.silhouettes
        self.result_canvas.Refresh()
    
    def OnShapeFill( self, evt ):
        self.annotate_canvas.shape_fill = not self.annotate_canvas.shape_fill
        self.annotate_canvas.Refresh()
    
    def OnNoPinch( self, evt ):
        import annotations_sweep
        annotations_sweep.kNoPinch = not annotations_sweep.kNoPinch
        self.annotations.Beautify()
        self.result_canvas.Refresh()
    
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
        print('OnDebugging')
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
            SaveDirectoryPath( self.last_wd, '*' )
            
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
        wildcard = "Image files|" + ';'.join( [ '*' + ext for ext in kImageFileExtensions ] )
        dialog = wx.FileDialog( self, "Load background image...", self.last_wd, '', wildcard, style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            SaveDirectoryPath( self.last_wd, '*' )
            
            self.annotations.SetImage( str( dialog.GetPath() ) )
            self.GetTopLevelParent().Refresh()
    
    def OnFileMenuSaveScene( self, evt ):
        if self.save_filename is not None:
            self.annotations.Save( self.save_filename, overwrite = True )
            return
        
        suggested_name = 'annotated scene'
        if self.annotations.image_path: suggested_name = os.path.split( os.path.splitext( self.annotations.image_path )[0] )[1]
        
        dialog = wx.FileDialog(
            self, "Save scene file...", self.last_wd,
            suggested_name + '.annscene',
            "Annotation Scene|*.annscene",
            wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
            )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            SaveDirectoryPath( self.last_wd, '*' )
            
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
    
    def OnMenuReset( self ):
        self.SaveStateForUndo()
        print("reset!")
        self.annotations.Reset( keep_texture = True )
        ## Refresh all
        self.GetTopLevelParent().Refresh()
    
    def OnClick( self, event ):
        if event.GetId() == self.RESET_BUTTON:
            self.OnMenuReset()
        
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
    global kDebugging, kPaperFigureMode
    
    try:
        argv.remove( '-d' )
        kDebugging = True
    except: pass
    try:
        argv.remove( '--debug' )
        kDebugging = True
    except: pass
    
    try:
        argv.remove( '--paper' )
        kPaperFigureMode = True
    except: pass

def main():
    import sys
    argv = list( sys.argv[1:] )
    parse_options( argv )
    
    import wx
    app = wx.App()
    
    ## load the PyCrust shell
    ## NOTE: This generates the following errors and is very slow switching between the app window
    ##       and the PyCrust window.
    ## Sun Jun 15 20:19:36 yotronagon.local Python[32825] <Error>: CGContextRestoreGState: invalid context
    if kDebugging:
        from wxBasicGLCanvas import load_PyCrust
        crust = load_PyCrust( app )
    
    annotations = Annotations()
    
    frame = AnnotateFrame( None, annotations )
    frame.Show(True)
    
    for fname in argv:
        ext = os.path.splitext( fname )[1].lower()
        if ext in ( '.annscene', ):
            annotations.Load( fname )
            frame.save_filename = fname
            frame.Refresh()
        elif ext in kImageFileExtensions:
            annotations.SetImage( fname )
            frame.Refresh()
    
    app.annotateframe = frame
    app.MainLoop()
    app.Destroy()

if __name__ == '__main__': main()
