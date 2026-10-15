from OpenGL.GL import *
from OpenGL.GLU import *

from numpy import asarray, zeros, append

class gluTesselated( object ):
    def __init__( self, planar_line_loop ):
        planar_line_loop = asarray( planar_line_loop )
        if planar_line_loop.shape[1] == 2:
            planar_line_loop = append( planar_line_loop, zeros( ( planar_line_loop.shape[0], 1 ) ), axis = 1 )
        assert planar_line_loop.shape[1] == 3
        
        ## The list of simple contours; each contour is a tuple ( mode, sequence of 3d vertices ).
        self.contours = []
        ## Temporary storage for combined vertices.
        combined_vertices = []
        
        def myBegin( mode ):
            self.contours.append( [ mode, [] ] )
        def myEnd():
            #pass
            self.contours[-1][1] = asarray( self.contours[-1][1] )
            assert self.contours[-1][1].shape[1] == 3
        def myVertex3dv( vertex ):
            self.contours[-1][1].append( tuple( vertex ) )
        def myCombine( coords, d, w ):
            #print "myCombine:", coords
            combined_vertices.append( tuple( coords ) )
            return combined_vertices[-1]
        def myError( errno ):
            print("GLU Error:", gluErrorString( errno ))
        
        tess = gluNewTess()
        
        gluTessCallback( tess, GLU_TESS_BEGIN, myBegin )
        gluTessCallback( tess, GLU_TESS_VERTEX, myVertex3dv )
        gluTessCallback( tess, GLU_TESS_END, myEnd )
        gluTessCallback( tess, GLU_TESS_COMBINE, myCombine )
        gluTessCallback( tess, GLU_TESS_ERROR, myError )
        gluTessProperty( tess, GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_NONZERO )
        
        gluTessBeginPolygon( tess, None )
        gluTessBeginContour( tess )
        for pt in planar_line_loop:
            gluTessVertex( tess, pt, pt )
        gluTessEndContour( tess )
        gluTessEndPolygon( tess )
        
        gluDeleteTess( tess )
    
    def draw( self ):
        for mode, vertices in self.contours:
            if len( vertices ) == 0: continue
            
            glVertexPointer( 3, GL_DOUBLE, 0, asarray( vertices ) )
            glEnableClientState( GL_VERTEX_ARRAY )
            glDrawArrays( mode, 0, len( vertices ) )
            glDisableClientState( GL_VERTEX_ARRAY )

def test_loop( loop ):
    print("loop:")
    print(loop)
    tess = gluTesselated( loop )
    print("tesselated contours:")
    print(tess.contours)

def test():
    ## simple square
    simple = [
        ( 0,0 ),
        ( 1,0 ),
        ( 1,1 ),
        ( 0,1 )
        ]
    test_loop( simple )
    
    ## self-intersecting square
    self_intersecting = [
        ( 0,0 ),
        ( 1,0 ),
        ( 1,1 ),
        ( 0,1 ),
        ( 0,0 ),
        ( .25, 0 ),
        ( .25, .75 ),
        ( .75, .75 ),
        ( .75, .25 ),
        ( 0, .25 )
        ]
    test_loop( self_intersecting )

def main():
    import OpenGL
    print("PyOpenGL version:", OpenGL.__version__)
    
    test()

if __name__ == '__main__': main()
