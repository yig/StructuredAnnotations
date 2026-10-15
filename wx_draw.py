from OpenGL.GL import *
from myarray import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.

from wxBasicGLCanvas import BasicGLCanvas
from ui_actions import StrokeHandler, RotateHandler, forward_events
from mesh_gui import TriMesh, Camera, set_headlight, draw_mesh_faces

class FiberCanvas( BasicGLCanvas ):
    def __init__( self, parent ):
        BasicGLCanvas.__init__( self, parent, debugging = True )
        
        self.stroke = []
        
        self.camera = Camera()
        self.camera.proj = Camera.Ortho()
        self.mesh = None
        
        debugging = False
        self.sketcher = StrokeHandler( self.stroke, self.doPointAdded, self.doStrokeComplete, self.event2world, debugging = debugging )
        self.rotator = RotateHandler( self.camera, self.Refresh, self.event2world, debugging = debugging )
        forward_events( self, self.sketcher )
    
    def doPointAdded( self, _ ):
        #print 'callback param:', _
        
        self.Refresh()
    
    def doStrokeComplete( self, _ ):
        #print 'callback param:', _
        
        import fibermesh
        self.stroke = fibermesh.smooth_curve( self.stroke )
        #print 'self.stroke:', self.stroke
        
        fibermesh.save_stroke( self.stroke, 'last_stroke.py' )
        
        self.mesh = fibermesh.inflate( self.stroke )
        #self.mesh = fibermesh.test2()
        
        forward_events( self, self.rotator )
        self.Refresh()
    
    def doDraw( self ):
        glClearColor( .9,.9,.9,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        '''
        glColor3f( 1., 0., 0. )
        glBegin( GL_TRIANGLES )
        glVertex2f( -.5, -.5 )
        glVertex2f( .5, -.5 )
        glVertex2f( 0, .5 )
        glEnd()
        '''
        
        #print 'self.stroke:', self.stroke
        if len( self.stroke ) >= 2: self.drawStroke()
        if self.mesh is not None: self.drawMesh()
        
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
        
        self.camera.Apply()
        ## TODO: Something with the light.  Right now, we're getting OpenGL's default light.
        
        set_headlight( self.camera )
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE )
        glShadeModel( GL_SMOOTH )
    
    def drawStroke( self ):
        self.prepare2D()
        
        glLineWidth( 3 )
        glColor4ub( 174, 174, 0, 200 )
        glBegin( GL_LINE_STRIP )
        for pt in self.stroke:
            glVertex2f( *pt )
        glEnd()
    
    def drawMesh( self, draw_faces = True, draw_edges = False ):
        self.prepare3D()
        
        if draw_faces:
            glEnable( GL_LIGHTING )
            
            glEnable( GL_COLOR_MATERIAL )
            glColor3ub( 200, 200, 255 )
            draw_mesh_faces( self.mesh, draw_edges )

def main():
    import wx
    # app = wx.App() ## <-- on OS X tracebacks and stdout/stderr appear in a dialog
    app = wx.PySimpleApp()
    
    ## load the PyCrust shell
    ## NOTE: This generates the following errors and is very slow switching between the app window
    ##       and the PyCrust window.
    ## Sun Jun 15 20:19:36 yotronagon.local Python[32825] <Error>: CGContextRestoreGState: invalid context
    # crust = load_PyCrust( app )
    
    frame = wx.Frame( None, -1, 'simple view', size=(512,512) )
    canvas = FiberCanvas( frame )
    frame.Show(True)
    
    #frame = ExampleFrame( None, -1, 'example view', size = (512,512) )
    
    app.MainLoop()
    app.Destroy()

if __name__ == '__main__': main()
