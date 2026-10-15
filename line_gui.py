#!/usr/bin/env python

import sys
from GLUTWindow import GLUTWindow
from myarray import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from glhelpers import draw_linestrip, draw_lineloop

try:
    from OpenGL.GLUT import *
    from OpenGL.GL import *
    from OpenGL.GLU import *
except:
    print('''
ERROR: PyOpenGL not installed properly. 
        ''')
    sys.exit()

from mesh_gui import Camera

class LineWindow( GLUTWindow ):
    def __init__( self, **kwargs ):
        GLUTWindow.__init__( self )
        
        kwargs.setdefault( 'background_color', ( .3, .3, .3 ) )
        kwargs.setdefault( 'line', [] )
        default_camera = Camera()
        default_camera.proj = Camera.Ortho( -3, 3, -3, 3 )
        kwargs.setdefault( 'camera', default_camera )
        kwargs.setdefault( 'line_is_closed', False )
        
        self.background_color = kwargs[ 'background_color' ]
        self.line = kwargs[ 'line' ]
        self.camera = kwargs[ 'camera' ]
        self.line_is_closed = kwargs[ 'line_is_closed' ]
        
        ## I don't want dynamic binding here because
        ## subclasses' constructors haven't run yet
        LineWindow.reset( self )
    
    def reset( self ):
        print('Resetting')
        
        self.postRedisplay()
    
    def displayFunc( self ):
        col = self.background_color
        glClearColor( col[0], col[1], col[2], 1.0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        glEnable( GL_DEPTH_TEST )
        
        self.camera.Apply()
        glDisable( GL_LIGHTING )
        
        glLineWidth( 8 )
        glColor3ub( 200, 200, 255 )
        if self.line_is_closed:
            draw_lineloop( self.line )
        else:
            draw_linestrip( self.line )
    
    def motionFunc( self, x, y ):
        ## Very simply map the camera onto +z half of the sphere x^2 + y^2 + z^2 = 1
        
        view = glGetIntegerv( GL_VIEWPORT )
        assert view[0] == 0
        assert view[1] == 0
        
        wx = view[2]/2.
        wy = view[3]/2.
        
        ## Flip y so that 0,0 is in the bottom-left.
        y = view[3] - y - 1
        
        radius = min( view[2], view[3] )/2.
        sphere_x = ( x - wx ) / radius
        sphere_y = ( y - wy ) / radius
        
        ## Project onto the closest point on the sphere
        if sphere_x*sphere_x + sphere_y*sphere_y > 1:
            denom = sqrt( sphere_x*sphere_x + sphere_y*sphere_y )
            sphere_x /= denom
            sphere_y /= denom
        
        sphere_z = sqrt( max( 0, 1 - sphere_x*sphere_x - sphere_y*sphere_y ) )
        
        eye_mag = mag( self.camera.eye )
        self.camera.eye = eye_mag * dir( arrayf( ( sphere_x, sphere_y, sphere_z ) ) )
        # print 'New eye position: ', self.camera.eye
        
        glutPostRedisplay()
    
    def mouseFunc( self, button, state, x,y ):
        if state == GLUT_DOWN: pass
        pass
    
    def passiveMotionFunc( self, x,y ):
        pass
    
    def keyboardFunc( self, key, x, y ):
        if 'r' == key: self.reset()
        else: GLUTWindow.keyboardFunc( self, key, x, y )
        '''
        if key == 'q' or key == 'Q': sys.exit(0)
        if key == 'c': self.captureScreen()
        if key == '\\': self.toggleFullScreen()
        '''
        
        glutPostRedisplay()

def main():
    glutInit( sys.argv )
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH )
    
    w1 = LineWindow()
    w1.setWindowTitle( "Helix" )
    ## Make a Helix with two twists.  We'll make its axis parallel to the y axis,
    ## so it appears vertically oriented in the default view.
    w1.line = [
        ( cos( 2*pi*t ), 2*t, sin( 2*pi*t ) )
        for t in linspace( -1, 1, 100 )
        ]
    
    glutMainLoop()

# call main() if the script is being run
if __name__ == '__main__': main()
