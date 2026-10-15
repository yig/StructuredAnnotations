#!/usr/bin/env python

import sys

try:
    from OpenGL.GLUT import *
    from OpenGL.GL import *
    from OpenGL.GLU import *
except:
    print('''
ERROR: PyOpenGL not installed properly. 
        ''')
    sys.exit()

from GLUTWindow import GLUTWindow

class StrokeWindow( GLUTWindow ):
    def __init__( self, stroke ):
        GLUTWindow.__init__( self )
        
        self.stroke = stroke
    
    def displayFunc( self ):
        col = self.background_color
        
        glDisable( GL_DEPTH_TEST )
        
        glClearColor( col[0], col[1], col[2], 1.0 )
        glClear( GL_COLOR_BUFFER_BIT )
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluOrtho2D( -1, 1, -1, 1 )
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        
        glLineWidth( 3 )
        glColor3f( 0, 1, 1 )
        
        glBegin( GL_LINE_LOOP )
        for v in self.stroke: glVertex2f( *v )
        glEnd()
        
        glutSwapBuffers()

def test():
    glutInit( sys.argv )
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB )
    
    exec('from %s import stroke' % (sys.argv[1],))
    window = StrokeWindow( stroke )
    
    glutMainLoop()


# call main() if the script is being run
if __name__ == '__main__':
    test()
