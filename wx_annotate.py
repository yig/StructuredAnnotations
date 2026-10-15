from OpenGL.GL import *
from myarray import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.

import wx
from wxBasicGLCanvas import BasicGLCanvas
from ui_actions import StrokeHandler, RotateHandler, forward_events
from wx_util import Timers
from mesh_gui import TriMesh, Camera, set_headlight, draw_mesh_faces

from glhelpers import *

kIterations = 10

class AnnotateData( object ):
    def __init__( self ):
        self.Reset()
    
    def Reset( self ):
        self.silhouette = None
        self.cross_sections = []
        self.cross_sections_3d = []
        self.mesh = None
        self.llconstraints = None
    
    def SetSilhouette( self, stroke ):
        import fibermesh
        from annotations_rotund import mesh_llconstraints_from_constraints
        
        self.silhouette = fibermesh.smooth_curve( stroke )
        #print 'self.data.silhouette:', self.data.silhouette
        
        fibermesh.save_stroke( self.silhouette, 'last_stroke.py' )
        
        #self.mesh = fibermesh.inflate( self.data.silhouette, 1 )
        #self.mesh = fibermesh.test2()
        
        self.mesh, self.llconstraints = mesh_llconstraints_from_constraints(
            self.silhouette, self.cross_sections
            )[:2]
    
    def AddCrossSection( self, stroke ):
        import fibermesh
        from annotations_rotund import mesh_llconstraints_from_constraints
        
        self.cross_sections.append( fibermesh.smooth_open_stroke( stroke ) )
        self.mesh, self.llconstraints, self.cross_sections_3d = mesh_llconstraints_from_constraints(
            self.silhouette, self.cross_sections
            )
    
    def Iterate( self, N = 1 ):
        '''
        Returns True if something changed, False otherwise.
        (Will always return True if there is a mesh with constraints.)
        '''
        
        assert N >= 0
        if N == 0: return False
        
        if self.mesh and self.llconstraints:
            import fibermesh
            fibermesh.optimize_surface2( self.mesh, self.llconstraints, N )
            return True
        else:
            print("No mesh")
            return False

class AnnotateCanvas( BasicGLCanvas ):
    def __init__( self, parent, data, pos = wx.DefaultPosition, size = wx.DefaultSize ):
        BasicGLCanvas.__init__( self, parent, pos, size, debugging = True )
        
        self.data = data
        
        self.stroke = []
        
        debugging = False
        self.sketcher = StrokeHandler( self.stroke, self.doPointAdded, self.doStrokeComplete, self.event2world, debugging = debugging )
        forward_events( self, self.sketcher )
    
    def GetBestSize( self ):
        return (500,500)
    
    def doPointAdded( self, _ ):
        #print 'callback param:', _
        
        self.Refresh()
    
    def doStrokeComplete( self, _ ):
        #print 'callback param:', _
        
        ## Reset if not enough points in the silhouette
        if len( self.stroke ) < 3:
            del self.stroke[:]
            return
        
        import fibermesh
        
        ## Did we just complete the silhouette?
        if self.data.silhouette is None:
            self.data.SetSilhouette( self.stroke )
            self.data.Iterate( 1 )
            
        ## We just completed a cross section.
        else:
            ## TODO: If the stroke starts and ends near an existing cross section, overdraw (replace) it.
            self.data.AddCrossSection( self.stroke )
            self.data.Iterate( 1 )
        
        del self.stroke[:]
        ## Refresh all
        self.GetTopLevelParent().Refresh()
    
    def doDraw( self ):
        glClearColor( .9,.9,.9,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        #print 'self.stroke:', self.stroke
        self.prepare2D()
        
        if self.data.silhouette is not None:
            glLineWidth( 3 )
            glColor4ub( 174, 174, 0, 200 )
            draw_lineloop( self.data.silhouette )
        
        for cross_section in self.data.cross_sections:
            glLineWidth( 2 )
            glColor4ub( 0, 174, 174, 200 )
            draw_linestrip( cross_section )
        
        glLineWidth( 5 )
        glColor4ub( 0, 0, 0, 200 )
        draw_linestrip( self.stroke )
        
        self.SwapBuffers()
    
    def prepare2D( self ):
        glDisable( GL_DEPTH_TEST )
        glDisable( GL_LIGHTING )
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(-1, 1, -1, 1, -1, 1)
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()


class ResultCanvas( BasicGLCanvas ):
    def __init__( self, parent, data, pos = wx.DefaultPosition, size = wx.DefaultSize ):
        ## GLUT = True because this uses draw_text().
        BasicGLCanvas.__init__( self, parent, pos, size, debugging = True, GLUT = True )
        
        self.camera = Camera()
        self.camera.proj = Camera.Ortho()
        
        self.data = data
        
        debugging = False
        self.rotator = RotateHandler( self.camera, self.Refresh, self.event2world, debugging = debugging )
        forward_events( self, self.rotator )
    
    def GetBestSize( self ):
        return (500,500)
    
    def doDraw( self ):
        glClearColor( .5,.5,.5,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        self.prepare3D()
        
        ## Draw the annotations in 3D
        glDisable( GL_LIGHTING )
        
        if self.data.silhouette is not None:
            glLineWidth( 3 )
            glColor4ub( 174, 174, 0, 200 )
            draw_lineloop( self.data.silhouette )
        
        ## Draw the cross section curves
        glLineWidth( 2 )
        glColor4ub( 90, 174, 174, 200 )
        for cross_section in self.data.cross_sections_3d:
            draw_linestrip( cross_section )
        
        ## Draw vertex number labels
        if self.data.mesh is not None:
            glColor4f( 0,0,0, 1 )
            for vi in range(len( self.data.mesh.vs )):
                draw_text( self.data.mesh.vs[vi] + .01 * self.data.mesh.vertex_normals[vi], str(vi) )
        
        if self.data.mesh is not None:
            glEnable( GL_LIGHTING )
            
            glEnable( GL_COLOR_MATERIAL )
            glColor3ub( 200, 200, 255 )
            draw_mesh_faces( self.data.mesh, False )
        
        self.SwapBuffers()
    
    def prepare2D( self ):
        glDisable( GL_DEPTH_TEST )
        glDisable( GL_LIGHTING )
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(-1, 1, -1, 1, -1, 1)
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
    
    def prepare3D( self ):
        glEnable( GL_DEPTH_TEST )
        glEnable( GL_LIGHTING )
        
        self.center_camera_on_mesh()
        self.camera.Apply()
        ## TODO: Something with the light.  Right now, we're getting OpenGL's default light.
        
        set_headlight( self.camera )
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE )
        glShadeModel( GL_SMOOTH )
    
    def event2world( self, pos ):
        ## Use this opportunity to make sure the camera is focused on the object's COM.
        self.center_camera_on_mesh()
        return BasicGLCanvas.event2world( self, pos )
    
    def center_camera_on_mesh( self ):
        if self.data.mesh:
            vs = self.data.mesh.vs
            self.camera.center = asarray( vs ).sum( 0 ) / len( vs )


class AnnotateFrame( wx.Frame ):
    
    RESET_BUTTON = 10
    ITERATE_BUTTON = 20
    ITERATE_INF_BUTTON = 30
    
    def __init__( self, parent, data ):
        wx.Frame.__init__( self, parent, wx.ID_ANY, 'Annotation Modeler' )
        
        self.data = data
        
        panel = wx.Panel( self, wx.ID_ANY, size=(1000,530) )
        self.sizer = wx.BoxSizer( wx.VERTICAL )
        
        ## Add the side-by-side gl views
        glviews = wx.BoxSizer( wx.HORIZONTAL )
        ## Create the AnnotateCanvas.
        canvas = AnnotateCanvas( panel, data )
        self.annotate_canvas = canvas
        canvas.SetSize( canvas.GetBestSize() )
        glviews.Add( canvas, 1, wx.EXPAND )
        ## Create the ResultCanvas.
        canvas = ResultCanvas( panel, data )
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
        btn = wx.Button( panel, self.ITERATE_BUTTON, "Iterate Once" )
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
    
    def OnClick( self, event ):
        if event.GetId() == self.RESET_BUTTON:
            print("reset!")
            self.data.Reset()
            ## Refresh all
            self.GetTopLevelParent().Refresh()
        
        elif event.GetId() == self.ITERATE_BUTTON:
            print("iterate!")
            self.IterateOnce()
        
        elif event.GetId() == self.ITERATE_INF_BUTTON:
            if hasattr( self, 'iterate_inf_timer_id' ):
                Timers.CancelTimer( self.iterate_inf_timer_id )
                del self.iterate_inf_timer_id
            else:
                self.iterate_inf_timer_id = Timers.RequestTimer( self.IterateOnce, 10, once = False )
    
    def IterateOnce( self ):
        if self.data.Iterate( 1 ):
            ## Refresh all
            self.GetTopLevelParent().Refresh()

def main():
    import wx
    # app = wx.App() ## <-- on OS X tracebacks and stdout/stderr appear in a dialog
    app = wx.PySimpleApp()
    
    ## load the PyCrust shell
    ## NOTE: This generates the following errors and is very slow switching between the app window
    ##       and the PyCrust window.
    ## Sun Jun 15 20:19:36 yotronagon.local Python[32825] <Error>: CGContextRestoreGState: invalid context
    # crust = load_PyCrust( app )
    
    data = AnnotateData()
    
    frame = AnnotateFrame( None, data )
    frame.Show(True)
    
    #frame = wx.Frame( None, -1, 'annotate', size=(512,512) )
    #frame.Show(True)
    
    #frame = ExampleFrame( None, -1, 'example view', size = (512,512) )
    
    app.MainLoop()
    app.Destroy()

if __name__ == '__main__': main()
