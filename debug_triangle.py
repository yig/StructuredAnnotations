## To use: output a 2d array to inspect using array.tofile( 'foo' ).
## then call this script with 'foo' as the argument.

from OpenGL.GLUT import *
from mesh_gui import TriMeshWindow
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.

glutInit( sys.argv )
glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH )

w1 = TriMeshWindow()
w1.setWindowTitle( "mesh" )

w1.mesh.vs.extend( [
    ( 0.502082109, -0.243724272, -0.0352211744 ),
    ( 0.495359987, -0.234970853, -0.0346866958 ),
    ( 0.498197854, -0.239017472, -0.040660087 )
    ] )
w1.mesh.vs.extend( [
    ( 0.471493304, -0.232644111, -0.0345747136 ),
    ( 0.486419261, -0.239733845, -0.035977155 ),
    ( 0.477539539, -0.230001435, -0.0398192294 )
    ] )
w1.mesh.faces.extend( [ (0,1,2), (3,4,5) ] )

w1.camera.center = average( asarray( w1.mesh.vs ), axis = 0 )
w1.camera.proj.fovy /= 100.

glutMainLoop()
