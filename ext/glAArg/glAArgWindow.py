#!/usr/bin/env python3

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

from glAArg import *

class glAArgWindow( object ):
    def __init__( self, **kwargs ):
        
        kwargs.setdefault( 'size', ( 512, 512 ) )
        kwargs.setdefault( 'title', 'GLUT Window' )
        kwargs.setdefault( 'background_color', ( 0., .3, .3 ) )
        
        self.size = kwargs[ 'size' ]
        self.title = kwargs[ 'title' ]
        self.background_color = kwargs[ 'background_color' ]
        
        glutInitWindowSize( self.size[0], self.size[1] )
        self.id = glutCreateWindow( self.title )
        
        glutKeyboardFunc( self.keyboardFunc )
        glutDisplayFunc( self.displayFunc )
        
        glAAPrepareBlendingState()
        glAAPrepareTextureState()
        
        glAAInit()
        glAAGenerateAATex(0.0, 0, 0)
        glAAGenerateAATex(0.5, 1, 0)
        glAAEnable(GLAA_VERTEX_ARRAY)
        
        ## I don't want dynamic binding here because
        ## subclasses' constructors haven't run yet
        glAArgWindow.reset( self )
    
    def reset( self ):
        print('Resetting')
        
        glutPostRedisplay()
    
    def displayFunc( self ):
        col = self.background_color
        glClearColor( col[0], col[1], col[2], 1.0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        glDisable( GL_DEPTH_TEST )
        
        #glAAPrepareBlendingState()
        #glAAPrepareTextureState()
        
        ## Set up the viewing transform.
        glMatrixMode( GL_PROJECTION )
        glLoadIdentity()
        glOrtho( -128, 128, -128, 128, -128, 128 )
        
        glMatrixMode( GL_MODELVIEW )
        glLoadIdentity()
        
        ## Use the first anti-aliasing parameters.
        glBindTexture(GL_TEXTURE_2D, glAA_texture[0])
        
        glAAColor3f( 200, 200, 255 )
        #glAAColor4f( 200, 200, 255, 100 )
        
        ## Draw a line.
        glAALineWidth( 10 )
        glAABegin(GL_LINES)
        glAAVertex2f(123,     45)
        glAAVertex2f(123.4, -111)
        glAAEnd()
        glAAFlush()
        
        
        glAAColor3f( 200, 0, 255 )
        
        ## Translate.
        glPushMatrix()
        glTranslatef(-100, 0, 0)
        
        ## Draw a point.
        glAAPointSize( 20 )
        glAABegin(GL_POINTS)
        glAAVertex2f(50,      50)
        glAAEnd()
        glAAFlush()
        
        ## Un-translate.
        glPopMatrix()
        
        ## Swap buffers.
        glutSwapBuffers()
    
    def keyboardFunc( self, key, x, y ):
        if 'r' == key: self.reset()
        if key == 'q' or key == 'Q': sys.exit(0)
        
        glutPostRedisplay()

def main():
    glutInit( sys.argv )
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH )
    
    w1 = glAArgWindow( title = "glAArg test" )
    
    glutMainLoop()

# call main() if the script is being run
if __name__ == '__main__': main()
