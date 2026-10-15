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

real_t = ctypes.c_double
lib = ctypes.cdll.LoadLibrary( 'edge_distances.' + platform_shared_library_suffix() )

__trace_point = lib.trace_point
__trace_point.argtypes = \
    [
        ctypes.c_int,
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.c_int,
        real_t,
        real_t,
        ctypes.POINTER( real_t )
    ]
__trace_point.restype = None
import edge_distances
def trace_point( point, linestrip, distance, dt = edge_distances.kDefaultTraceDt ):
    '''
    Given an N-d point 'point' and an N-d sequence of points forming a line strip 'linestrip',
    traces 'point' in the smoothed distance gradient field implied by 'linestrip' a
    distance 'distance'.
    Optional parameter 'dt' controls the step size of each step through the smoothed distance
    gradient field.
    
    Similar to edge_distances.trace_pts_in_smooth_distance_field(), except traces a single
    point at a time and there is no optional trajectory output variable.
    '''
    
    point = to_ctypes_array( point, real_t )
    linestrip = to_ctypes_array( linestrip, real_t )
    
    assert len( point ) == linestrip.shape[1]
    
    trace = to_ctypes_array( zeros( point.shape ), real_t )
    
    __trace_point(
        len( point ),
        ctypes_array_ref( point, real_t ),
        ctypes_array_ref( linestrip, real_t ),
        len( linestrip ),
        distance,
        dt,
        ctypes_array_ref( trace, real_t )
        )
    
    return trace


def test():
    cases = []
    
    ## a square
    edges = [(.25, .25), (.75, .25), (.75, .75), (.25, .75), (.25, .25)]
    pts = [(.3,.3)]
    cases.append( (pts, edges) )
    
    for pts, edges in cases:
        edges = array( edges, real_t )
        pts = array( pts, real_t )
        
        print('pts:', pts)
        print('edges:', edges)
        
        ## this shouldn't go anywhere (we're asking it to return immediately)
        traces = asarray( [ trace_point( point, edges, 0. ) for point in pts ] )
        assert ( (traces - pts)**2 ).sum() <= 1e-5
        
        print('traces:', traces)
        print('distance travelled:', sqrt( ( (traces - pts)**2 ).sum(-2) ))
        
        ## ideally this should travel distance 1
        traces = asarray( [ trace_point( point, edges, 1. ) for point in pts ] )
        #assert ( ( ( (traces - pts)**2 ).sum(-2) - ones( pts.shape ) )**2 ).sum() <= 1e-5
        
        print('traces:', traces)
        print('distance travelled:', sqrt( ( (traces - pts)**2 ).sum(-2) ))

def main():
    test()

if __name__ == "__main__": main()
