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
lib = ctypes.cdll.LoadLibrary( 'face_stitch_pattern.' + platform_shared_library_suffix() )

__face_stitch_pattern = lib.face_stitch_pattern
__face_stitch_pattern.argtypes = \
    [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER( ctypes.c_int )
    ]
__face_stitch_pattern.restype = None
def face_stitch_pattern( num_cross_sections, num_cross_section_verts ):
    '''
    Given a number of cross section 'num_cross_sections' and number of vertices in each cross
    section 'num_cross_section_verts', returns an array with shape
        ( ( num_cross_section_verts ) * ( num_cross_sections - 1 ) * 2, 3 )
    containing stitched triangle indices exactly the same way sweep_surface.stitch_cross_sections()
    would have done it.
    '''
    
    assert num_cross_sections >= 0
    assert num_cross_section_verts >= 0
    
    faces = to_ctypes_array( empty( ( ( num_cross_section_verts ) * ( num_cross_sections - 1 ) * 2, 3 ) ), ctypes.c_int )
    
    __face_stitch_pattern(
        num_cross_sections,
        num_cross_section_verts,
        ctypes_array_ref( faces, ctypes.c_int )
        )
    
    return faces

def test():
    pass

def main():
    test()

if __name__ == "__main__": main()
