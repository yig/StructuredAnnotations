from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.

try:
    import ctypes
except:
    print('''
ERROR: ctypes not installed properly. (ctypes needed for laplacian_editing_1d())
        ''')
    import sys
    sys.exit()

def platform_shared_library_suffix():
    import sys
    result = 'so'
    if 'win' in sys.platform.lower(): result = 'dll'
    ## No else if, because we want darwin to override win (which is a substring of darwin)
    if 'darwin' in sys.platform.lower(): result = 'dylib'
    return result

def to_ctypes_array( sequence, ctypes_type ):
    '''
    Returns a contiguous ctypes array suitable for extracting
    a pointer (with .ctypes.data_as() or ctypes_array_ref()) to pass to a C function.
    '''
    return ascontiguousarray( asarray( sequence, dtype = ctypes_type ) )
def ctypes_array_ref( ctypes_array, ctypes_type ):
    '''
    Returns a ctypes.POINTER to the internal data of 'ctypes_array' of type 'ctypes_type'.
    
    WARNING: Do not chain the calls to to_ctypes_array() and ctypes_array_ref() as in
               ctypes_array_ref( to_ctypes_array( sequence, ctypes_type ), ctypes_type ).
             ctypes_array_ref returns a pointer to the internal numpy storage of the array and
             cannot keep a reference to the array, so it will be garbage collected if these
             functions are chained together.
             Instead, keep a temporary reference to the result of to_ctypes_array():
               contig = to_ctypes_array( sequence, ctypes_type )
               ctypes_array_ref( contig, ctypes_type )
    '''
    return ctypes_array.ctypes.data_as( ctypes.POINTER( ctypes_type ) )

real_t = ctypes.c_float
lib = ctypes.cdll.LoadLibrary( 'pnpoly.' + platform_shared_library_suffix() )

# int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy);
__pnpoly = lib.pnpoly
__pnpoly.argtypes = \
    [
        ctypes.c_int,
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        real_t,
        real_t
    ]
__pnpoly.restype = ctypes.c_int
def pnpoly( line_loop_2D, point_2D ):
    '''
    Given a 2D line loop as a sequence of 2D values 'line_loop_2D' and
    a 2D point 'point_2D',
    returns 1 if the point is inside the polygon,
    and returns 0 otherwise.
    '''
    
    assert len( line_loop_2D ) > 2
    assert len( point_2D ) == 2
    
    line_loop_2D = asarray( line_loop_2D )
    line_loop_2D_xs = to_ctypes_array( line_loop_2D[:,0], real_t )
    line_loop_2D_ys = to_ctypes_array( line_loop_2D[:,1], real_t )
    
    result = __pnpoly(
        len( line_loop_2D ),
        ctypes_array_ref( line_loop_2D_xs, real_t ),
        ctypes_array_ref( line_loop_2D_ys, real_t ),
        point_2D[0],
        point_2D[1]
        )
    
    ## Convert to an 'int' so we don't leak a ctypes type.
    return int( result )

def line_loop_contains_point( line_loop_2D, point_2D ):
    '''
    Given a 2D line loop as a sequence of 2D values 'line_loop_2D' and
    a 2D point 'point_2D',
    returns True if 'point_2D' lies strictly inside 'line_loop_2D'
    (not on an edge) and False otherwise.
    '''
    return bool( pnpoly( line_loop_2D, point_2D ) )

def test():
    line_loop = [ ( -1,-1 ), ( 1, -1 ), ( 1, 1 ), ( -1, 1 ) ]
    test_pts = [ ( 0, 0 ), ( -1, -1 ), ( 2, 0 ) ]
    
    print('line loop:', line_loop)
    for pt in test_pts:
        print('test point:', pt)
        print('pnpoly():', pnpoly( line_loop, pt ))
        print('line_loop_contains_point():', line_loop_contains_point( line_loop, pt ))

def main():
    test()

if __name__ == "__main__": main()
