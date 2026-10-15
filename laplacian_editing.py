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


'''
typedef double Real;
typedef unsigned int uint;
typedef unsigned char BOOL;

void cotangentWeights(
    const Real* center,
    const Real* one_ring, const uint valence,
    const BOOL has_boundary,
    Real* cot_alpha_out, Real* cot_beta_out, Real* area_out
    );

void Hn_coefficients(
    const int center,
    const int* one_ring, const uint valence,
    const BOOL has_boundary,
    const Real* cot_alpha, const Real* cot_beta, const Real area,
    int* indices_out, Real* values_out
    );
'''

Real = ctypes.c_double
uint = ctypes.c_uint
BOOL = ctypes.c_ubyte
liblaplacian_editing = ctypes.cdll.LoadLibrary( 'liblaplacian_editing.dylib' )

__cotangentWeights = liblaplacian_editing.cotangentWeights
__cotangentWeights.argtypes = \
    [
        ctypes.POINTER( Real ),
        ctypes.POINTER( Real ), uint,
        BOOL,
        ctypes.POINTER( Real ), ctypes.POINTER( Real ), ctypes.POINTER( Real )
    ]
__cotangentWeights.restype = None

__Hn_coefficients = liblaplacian_editing.Hn_coefficients
__Hn_coefficients.argtypes = \
    [
        ctypes.c_int,
        ctypes.POINTER( ctypes.c_int ), uint,
        BOOL,
        ctypes.POINTER( Real ), ctypes.POINTER( Real ), Real,
        ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( Real )
    ]
__Hn_coefficients.restype = None

def cotangentWeights( center, one_ring, has_boundary ):
    '''
    Same as "laplacian_editing_wrapper.h" cotangentWeights(), except valence is omitted.
    Returns 3 values: cot_alpha_out (array), cot_beta_out (array), area_out (float).
    '''
    
    center = ascontiguousarray( asarray( center, dtype = Real ) )
    valence = len( one_ring )
    one_ring = ascontiguousarray( asarray( one_ring, dtype = Real ).flatten() )
    
    cot_alpha_out = ascontiguousarray( zeros( valence, dtype = Real ) )
    cot_beta_out = ascontiguousarray( zeros( valence, dtype = Real ) )
    area_out = Real()
    
    __cotangentWeights(
        center.ctypes.data_as( ctypes.POINTER( Real ) ),
        one_ring.ctypes.data_as( ctypes.POINTER( Real ) ),
        valence,
        has_boundary,
        cot_alpha_out.ctypes.data_as( ctypes.POINTER( Real ) ),
        cot_beta_out.ctypes.data_as( ctypes.POINTER( Real ) ),
        ctypes.byref( area_out )
        )
    
    return cot_alpha_out, cot_beta_out, area_out.value

def Hn_coefficients( center, one_ring, has_boundary, cot_alpha, cot_beta, area ):
    '''
    Same as "laplacian_editing_wrapper.h" Hn_coefficients(), except valence is omitted.
    Returns 2 values: indices_out, values_out.
    '''
    
    valence = len( one_ring )
    assert valence == len( cot_alpha )
    assert valence == len( cot_beta )
    one_ring = ascontiguousarray( asarray( one_ring, dtype = ctypes.c_int ) )
    cot_alpha = ascontiguousarray( asarray( cot_alpha, dtype = Real ) )
    cot_beta = ascontiguousarray( asarray( cot_beta, dtype = Real ) )
    
    #from pydb import debugger
    #debugger()
    
    indices_out = ascontiguousarray( zeros( valence+1, dtype = ctypes.c_int ) )
    values_out = ascontiguousarray( zeros( valence+1, dtype = Real ) )
    
    __Hn_coefficients(
        center,
        one_ring.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) ),
        valence,
        has_boundary,
        cot_alpha.ctypes.data_as( ctypes.POINTER( Real ) ),
        cot_beta.ctypes.data_as( ctypes.POINTER( Real ) ),
        area,
        indices_out.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) ),
        values_out.ctypes.data_as( ctypes.POINTER( Real ) )
        )
    
    return indices_out, values_out


def main():
    center = (0,0,0)
    center_index = 0
    one_ring = [(1,0,0), (1,1,0), (0,1,0), (-1,1,0), (-1,0,0), (-1,-1,0), (0,-1,0), (1,-1,0)]
    one_ring_indices = arange(len(one_ring)) + 1
    has_boundary = False
    
    print('center index:', center_index)
    print('one ring indices:', one_ring_indices)
    print('center:', center)
    print('one_ring:', one_ring)
    cot_alpha, cot_beta, A = cotangentWeights( center, one_ring, has_boundary )
    print('cot_alpha:', cot_alpha)
    print('cot_beta:', cot_beta)
    print('A:', A)
    indices, values = Hn_coefficients( center_index, one_ring_indices, has_boundary, cot_alpha, cot_beta, A )
    print('indices:', indices)
    print('values:', values)

if __name__ == '__main__': main()
