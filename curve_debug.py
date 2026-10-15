## To use: output a 2d array to inspect using array.tofile( 'foo' ).
## then call this script with 'foo' as the argument.

from OpenGL.GLUT import *
from line_gui import LineWindow
from mesh_gui import Camera
from numpy import fromfile

glutInit( sys.argv )
glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH )

w1 = LineWindow( line_is_closed = True )
w1.setWindowTitle( "curve" )
w1.camera.proj = Camera.Ortho()

curve = fromfile( sys.argv[1] )
curve = curve.reshape( ( curve.shape[0]//2, 2 ) )
w1.line = curve

glutMainLoop()
