from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
#from tictoc import tic, toc
import os

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

'''
typedef double real_t;

void E_ijk(
    real_t* p_0, real_t* p_1, real_t* p_2,
    real_t* E_out,
    real_t* dE_dp_out,
    real_t* d2E_dp2_out
    );
'''

real_t = ctypes.c_double
## This dirname(__file__) prefix doesn't work after py2app puts this Python file inside a zip.
## py2app puts the dylib in the current working directory, so just access it.
# Eline_ad = ctypes.cdll.LoadLibrary( os.path.join( os.path.dirname( __file__ ), 'Eline_ad.' + platform_shared_library_suffix() ) )
Eline_ad = ctypes.cdll.LoadLibrary( 'Eline_ad.' + platform_shared_library_suffix() )

__Efunc_R3x2_argtypes = \
    [
        ## Two R3 vectors:
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ## Output parameters: energy, gradient, hessian
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t )
    ]
__Efunc_R3x2_restype = None
def Efunc_R3x2_and_derivatives( Efunc, p_0, p_1, compute_hessian = False ):
    '''
    Same as Efunc, except the output parameters are returned as a tuple of values:
    If 'compute_hessian' is false, returns a tuple with 2 values:
        ( the energy scalar, the gradient vector (with respect to p_i[2]) )
    else, returns a tuple with 3 values:
        ( the energy scalar, the gradient vector (with respect to p_i[2]), the hessian matrix (with respect to p_i[2], p_j[2]) )
    '''
    
    Efunc.argtypes = __Efunc_R3x2_argtypes
    Efunc.restype = __Efunc_R3x2_restype
    
    p_0 = to_ctypes_array( p_0, real_t )
    p_1 = to_ctypes_array( p_1, real_t )
    assert len( p_0 ) == 3
    assert len( p_1 ) == 3
    
    E = real_t(0.)
    dE_dp = to_ctypes_array( zeros( 6 ), real_t )
    if compute_hessian:
        d2E_dp2 = to_ctypes_array( zeros( (6,6) ), real_t )
        d2E_dp2_array_ref = ctypes_array_ref( d2E_dp2, real_t )
    else:
        d2E_dp2 = None
        d2E_dp2_array_ref = None
    
    #tic()
    Efunc(
        ctypes_array_ref( p_0, real_t ),
        ctypes_array_ref( p_1, real_t ),
        ctypes.byref( E ),
        ctypes_array_ref( dE_dp, real_t ),
        d2E_dp2_array_ref
        )
    #toc()
    
    if compute_hessian:
        #print ( E.value, dE_dp.round(2), d2E_dp2.round(2) )
        return ( E.value, dE_dp, d2E_dp2 )
    else:
        #print ( E.value, dE_dp.round(2) )
        return ( E.value, dE_dp )

__Efunc_R3x3_argtypes = \
    [
        ## Three R3 vectors:
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ## Output parameters: energy, gradient, hessian
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t )
    ]
__Efunc_R3x3_restype = None
def Efunc_R3x3_and_derivatives( Efunc, p_0, p_1, p_2, compute_hessian = False ):
    '''
    Same as Efunc, except the output parameters are returned as a tuple of values:
    If 'compute_hessian' is false, returns a tuple with 2 values:
        ( the energy scalar, the gradient vector (with respect to p_i[2]) )
    else, returns a tuple with 3 values:
        ( the energy scalar, the gradient vector (with respect to p_i[2]), the hessian matrix (with respect to p_i[2], p_j[2]) )
    '''
    
    Efunc.argtypes = __Efunc_R3x3_argtypes
    Efunc.restype = __Efunc_R3x3_restype
    
    p_0 = to_ctypes_array( p_0, real_t )
    p_1 = to_ctypes_array( p_1, real_t )
    p_2 = to_ctypes_array( p_2, real_t )
    assert len( p_0 ) == 3
    assert len( p_1 ) == 3
    assert len( p_2 ) == 3
    
    E = real_t(0.)
    dE_dp = to_ctypes_array( zeros( 9 ), real_t )
    if compute_hessian:
        d2E_dp2 = to_ctypes_array( zeros( (9,9) ), real_t )
        d2E_dp2_array_ref = ctypes_array_ref( d2E_dp2, real_t )
    else:
        d2E_dp2 = None
        d2E_dp2_array_ref = None
    
    #tic()
    Efunc(
        ctypes_array_ref( p_0, real_t ),
        ctypes_array_ref( p_1, real_t ),
        ctypes_array_ref( p_2, real_t ),
        ctypes.byref( E ),
        ctypes_array_ref( dE_dp, real_t ),
        d2E_dp2_array_ref
        )
    #toc()
    
    if compute_hessian:
        #print ( E.value, dE_dp.round(2), d2E_dp2.round(2) )
        return ( E.value, dE_dp, d2E_dp2 )
    else:
        #print ( E.value, dE_dp.round(2) )
        return ( E.value, dE_dp )

__Efunc_R3x6_argtypes = \
    [
        ## Six R3 vectors:
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ## Output parameters: energy, gradient, hessian
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t )
    ]
__Efunc_R3x6_restype = None
def Efunc_R3x6_and_derivatives( Efunc, p_0, p_1, p_2, u_0, u_1, u_2, compute_hessian = False ):
    '''
    Same as Efunc, except the output parameters are returned as a tuple of values:
    If 'compute_hessian' is false, returns a tuple with 2 values:
        ( the energy scalar, the gradient vector (with respect to p_i[2]) )
    else, returns a tuple with 3 values:
        ( the energy scalar, the gradient vector (with respect to p_i[2]), the hessian matrix (with respect to p_i[2], p_j[2]) )
    '''
    
    Efunc.argtypes = __Efunc_R3x6_argtypes
    Efunc.restype = __Efunc_R3x6_restype
    
    p_0 = to_ctypes_array( p_0, real_t )
    p_1 = to_ctypes_array( p_1, real_t )
    p_2 = to_ctypes_array( p_2, real_t )
    assert len( p_0 ) == 3
    assert len( p_1 ) == 3
    assert len( p_2 ) == 3
    
    u_0 = to_ctypes_array( u_0, real_t )
    u_1 = to_ctypes_array( u_1, real_t )
    u_2 = to_ctypes_array( u_2, real_t )
    assert len( u_0 ) == 3
    assert len( u_1 ) == 3
    assert len( u_2 ) == 3
    
    E = real_t(0.)
    dE_dp = to_ctypes_array( zeros( 9 ), real_t )
    if compute_hessian:
        d2E_dp2 = to_ctypes_array( zeros( (9,9) ), real_t )
        d2E_dp2_array_ref = ctypes_array_ref( d2E_dp2, real_t )
    else:
        d2E_dp2 = None
        d2E_dp2_array_ref = None
    
    #tic()
    Efunc(
        ctypes_array_ref( p_0, real_t ),
        ctypes_array_ref( p_1, real_t ),
        ctypes_array_ref( p_2, real_t ),
        ctypes_array_ref( u_0, real_t ),
        ctypes_array_ref( u_1, real_t ),
        ctypes_array_ref( u_2, real_t ),
        ctypes.byref( E ),
        ctypes_array_ref( dE_dp, real_t ),
        d2E_dp2_array_ref
        )
    #toc()
    
    if compute_hessian:
        #print ( E.value, dE_dp.round(2), d2E_dp2.round(2) )
        return ( E.value, dE_dp, d2E_dp2 )
    else:
        #print ( E.value, dE_dp.round(2) )
        return ( E.value, dE_dp )

def Efunc_R3x6_Noneable_and_derivatives( Efunc, p_0, p_1, p_2, u_0, u_1, u_2, compute_hessian = False ):
    '''
    Same as Efunc_R3x6_and_derivatives, except the input parameters can be None and will be passed
    through as NULL.
    '''
    
    Efunc.argtypes = __Efunc_R3x6_argtypes
    Efunc.restype = __Efunc_R3x6_restype
    
    if p_0 is not None: p_0 = to_ctypes_array( p_0, real_t )
    if p_1 is not None: p_1 = to_ctypes_array( p_1, real_t )
    if p_2 is not None: p_2 = to_ctypes_array( p_2, real_t )
    if p_0 is not None: assert len( p_0 ) == 3
    if p_1 is not None: assert len( p_1 ) == 3
    if p_2 is not None: assert len( p_2 ) == 3
    
    if u_0 is not None: u_0 = to_ctypes_array( u_0, real_t )
    if u_1 is not None: u_1 = to_ctypes_array( u_1, real_t )
    if u_2 is not None: u_2 = to_ctypes_array( u_2, real_t )
    if u_0 is not None: assert len( u_0 ) == 3
    if u_1 is not None: assert len( u_1 ) == 3
    if u_2 is not None: assert len( u_2 ) == 3
    
    E = real_t(0.)
    dE_dp = to_ctypes_array( zeros( 9 ), real_t )
    if compute_hessian:
        d2E_dp2 = to_ctypes_array( zeros( (9,9) ), real_t )
        d2E_dp2_array_ref = ctypes_array_ref( d2E_dp2, real_t )
    else:
        d2E_dp2 = None
        d2E_dp2_array_ref = None
    
    #tic()
    Efunc(
        ctypes_array_ref( p_0, real_t ) if p_0 is not None else None,
        ctypes_array_ref( p_1, real_t ) if p_1 is not None else None,
        ctypes_array_ref( p_2, real_t ) if p_2 is not None else None,
        ctypes_array_ref( u_0, real_t ) if u_0 is not None else None,
        ctypes_array_ref( u_1, real_t ) if u_1 is not None else None,
        ctypes_array_ref( u_2, real_t ) if u_2 is not None else None,
        ctypes.byref( E ),
        ctypes_array_ref( dE_dp, real_t ),
        d2E_dp2_array_ref
        )
    #toc()
    
    if compute_hessian:
        #print ( E.value, dE_dp.round(2), d2E_dp2.round(2) )
        return ( E.value, dE_dp, d2E_dp2 )
    else:
        #print ( E.value, dE_dp.round(2) )
        return ( E.value, dE_dp )

__Emag = Eline_ad.Emag
def Emag_and_derivatives( *args, **kwargs ):
    '''
    Same as Efunc_R3x2_and_derivatives, but the energy function is "Eline_adt.cpp" Emag().
    '''
    return Efunc_R3x2_and_derivatives( __Emag, *args, **kwargs )

__Ecurvature = Eline_ad.Ecurvature
def Ecurvature_and_derivatives( *args, **kwargs ):
    '''
    Same as Efunc_R3x3_and_derivatives, but the energy function is "Eline_adt.cpp" Ecurvature().
    '''
    return Efunc_R3x3_and_derivatives( __Ecurvature, *args, **kwargs )

__Ecurvature_undeformed = Eline_ad.Ecurvature_undeformed
def Ecurvature_and_derivatives_undeformed( *args, **kwargs ):
    '''
    Same as Efunc_R3x6_and_derivatives, but the energy function is "Eline_adt.cpp" Ecurvature_undeformed().
    '''
    return Efunc_R3x6_and_derivatives( __Ecurvature_undeformed, *args, **kwargs )

__Elaplacian_editing = Eline_ad.Elaplacian_editing
def Elaplacian_editing_and_derivatives( *args, **kwargs ):
    '''
    Same as Efunc_R3x3_and_derivatives, but the energy function is "Eline_adt.cpp" Elaplacian_editing().
    '''
    return Efunc_R3x3_and_derivatives( __Elaplacian_editing, *args, **kwargs )

__Elaplacian_editing_undeformed = Eline_ad.Elaplacian_editing_undeformed
def Elaplacian_editing_and_derivatives_undeformed( *args, **kwargs ):
    '''
    Same as Efunc_R3x6_and_derivatives, but the energy function is "Eline_adt.cpp" Elaplacian_editing_undeformed().
    '''
    return Efunc_R3x6_Noneable_and_derivatives( __Elaplacian_editing_undeformed, *args, **kwargs )

__Etest = Eline_ad.Etest
def Etest_and_derivatives( *args, **kwargs ):
    '''
    Same as Efunc_R3x3_and_derivatives, but the energy function is "Eline_adt.cpp" Etest().
    '''
    return Efunc_R3x3_and_derivatives( __Etest, *args, **kwargs )

def test_assembly( N = 5 ):
    print("=== test_assembly() ===")
    
    ## A horizontal line along the x axis from 0 to 1.
    line = zeros( (N,3) )
    line[:,0] = linspace( 0, 1, N )
    ## Actually, a sine wave in the xz plane.
    #line[:,2] = sin( 2*pi* line[:,0] )
    
    print("line:", line)
    
    E_total = 0.
    grad_total = zeros( 3*N )
    hess_total = zeros( (3*N,3*N) )
    for i in range( 1, N-1 ):
        
        E, grad, hess = Etest_and_derivatives( line[ i-1 ], line[ i ], line[ i+1 ], compute_hessian = True )
        
        print("E(%s):" % (i,), E)
        print("grad(%s):" % (i,))
        print(grad.round(3))
        print("hess(%s):" % (i,))
        print(hess.round(3))
        
        E_total += E
        
        grad_total[ (i-1)*3 : (i+2)*3 ] += grad
        
        hess_total[ (i-1)*3:(i+2)*3, (i-1)*3:(i+2)*3 ] += hess
    
    print("E:", E_total)
    print("grad:")
    print(grad_total.round(3))
    print("hess:")
    print(hess_total.round(3))

def test_curvature_3D( N = 5 ):
    print("=== test_curvature_3D() ===")
    
    ## A horizontal line along the x axis from 0 to 1.
    line = zeros( (N,3) )
    line[:,0] = linspace( 0, 1, N )
    ## Actually, a sine wave in the xz plane.
    #line[:,2] = sin( 2*pi* line[:,0] )
    
    print("line:", line)
    
    E_total = 0.
    grad_total = zeros( 3*N )
    hess_total = zeros( (3*N,3*N) )
    for i in range( 1, N-1 ):
        
        E, grad, hess = Ecurvature_and_derivatives( line[ i-1 ], line[ i ], line[ i+1 ], compute_hessian = True )
        
        print("E(%s):" % (i,), E)
        print("grad(%s):" % (i,))
        print(grad.round(3))
        print("hess(%s):" % (i,))
        print(hess.round(3))
        
        E_total += E
        
        grad_total[ (i-1)*3 : (i+2)*3 ] += grad
        
        hess_total[ (i-1)*3:(i+2)*3, (i-1)*3:(i+2)*3 ] += hess
    
    print("E:", E_total)
    print("grad:")
    print(grad_total.round(3))
    print("hess:")
    print(hess_total.round(3))

def test_curvature_z( N = 5 ):
    print("=== test_curvature_z() ===")
    
    ## A horizontal line along the x axis from 0 to 1.
    line = zeros( (N,3) )
    line[:,0] = linspace( 0, 1, N )
    ## Actually, a sine wave in the xz plane.
    #line[:,2] = sin( 2*pi* line[:,0] )
    
    print("line:", line)
    
    E_total = 0.
    grad_total = zeros( N )
    hess_total = zeros( (N,N) )
    for i in range( 1, N-1 ):
        
        E, grad, hess = Ecurvature_and_derivatives( line[ i-1 ], line[ i ], line[ i+1 ], compute_hessian = True )
        
        print("E(%s):" % (i,), E)
        print("grad(%s):" % (i,))
        print(grad.round(3))
        print("hess(%s):" % (i,))
        print(hess.round(3))
        
        E_total += E
        
        grad_total[ i-1 : i+2 ] += grad[2::3]
        
        hess_total[ i-1:i+2, i-1:i+2 ] += hess[2::3,2::3]
    
    print("E:", E_total)
    print("grad:")
    print(grad_total.round(3))
    print("hess:")
    print(hess_total.round(3))

def test_mag( N = 5 ):
    print("=== test_mag() ===")
    
    ## A horizontal line along the x axis from 0 to 1.
    line = zeros( (N,3) )
    line[:,0] = linspace( 0, 1, N )
    ## Actually, a sine wave in the xz plane.
    #line[:,2] = sin( 2*pi* line[:,0] )
    
    print("line:", line)
    
    E_total = 0.
    grad_total = zeros( 3*N )
    hess_total = zeros( (3*N,3*N) )
    for i in range( 0, N-1 ):
        
        E, grad, hess = Emag_and_derivatives( line[ i ], line[ i+1 ], compute_hessian = True )
        
        print("E(%s):" % (i,), E)
        print("grad(%s):" % (i,))
        print(grad.round(3))
        print("hess(%s):" % (i,))
        print(hess.round(3))
        
        E_total += E
        
        grad_total[ i*3 : (i+2)*3 ] += grad
        
        hess_total[ i*3:(i+2)*3, i*3:(i+2)*3 ] += hess
    
    print("E:", E_total)
    print("grad:")
    print(grad_total.round(3))
    print("hess:")
    print(hess_total.round(3))

def main():
    #test_assembly(4)
    
    #test_curvature_3D(4)
    #test_curvature_z(3)
    test_mag()

if __name__ == "__main__": main()
