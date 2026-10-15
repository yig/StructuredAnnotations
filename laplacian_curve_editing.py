from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.

kDebugging = False

def L_for_curve( curve, curve_is_closed, extra = None ):
    '''
    Given a 2d curve (sequence of 2d points) 'curve', and
    a boolean 'curve_is_closed' specifying whether or not 'curve' is closed,
    returns L from Sorkine et al 2004, defined in section 3 just after equation 2.
    (In fact, it is the operator whose action is cursive-L in equation 2.)
    
    L is the uniform laplacian operator, such that dot( L, asarray( curve ) ) is the
    array of laplacian vectors.
    
    If optional parameter 'extra' is not None, it is expected to be a dictionary
    that will be filled with the following output parameters:
        'A': the adjacency matrix
        'd': the array of number of neighbors for each point in the curve
    
    tested
    '''
    
    N = len( curve )
    A = zeros( (N, N), dtype = int )
    d = zeros( N, dtype = int )
    
    ## endpoints
    if curve_is_closed:
        A[0,-1] = 1
        A[0,1] = 1
        A[-1,0] = 1
        A[-1,-2] = 1
        
        d[0] = 2
        d[-1] = 2
    else:
        A[0,1] = 1
        A[-1,-2] = 1
        
        d[0] = 1
        d[-1] = 1
    
    ## interior points
    for i in range( 1, N-1 ):
        A[ i, i-1 ] = 1
        A[ i, i+1 ] = 1
        
        d[i] = 2
    
    Dinv = diag( 1./d )
    
    L = identity( N ) - dot( Dinv, A )
    
    if extra is not None:
        extra['d'] = d
        extra['A'] = A
    return L

def new_curve_from_curve_and_constraints( curve, curve_is_closed, constraints ):
    '''
    Given a 2d curve (sequence of 2d points) 'curve',
    a boolean 'curve_is_closed' specifying whether or not 'curve' is closed,
    and a sequence of constraint tuples 'constraints' where each tuple is of the form
        ( curve point index, constrained location ),
    returns a new sequence of curve points satisfying the given constraints.
    
    NOTE: There must be at least two constraints.
    
    tested
        with constraints at curve's current location, returns an identical curve.
        with real curves and two and three constraints (both end-points always constrained).
    '''
    
    curve = asarray( curve )
    N = curve.shape[0]
    
    assert curve.shape[1] == 2
    assert len( constraints ) >= 2
    assert False not in [ index >= 0 and index < N for index, pos in constraints ]
    
    extra = {}
    L = L_for_curve( curve, curve_is_closed, extra = extra )
    d = extra['d']
    neighbors = [
        extra['A'][i,:].nonzero()[0] for i in range(N)
        ]
    del extra
    
    laplacian_coordinates = dot( L, curve )
    
    ## I will organize the vector of degrees of freedom as follows:
    ## x1 ... xN y1 ... yN a0 ... aN w0 ... wN tx0 ... txN ty0 ... tyN
    ## A dictionary of offsets for accessing DOF indices:
    ## (To go from a DOF and its index to an index in the flattened vector of DOFs
    ## or system matrix, do o['x'] + index.)
    o = {
        'x': 0*N,
        'y': 1*N,
        'a': 2*N,
        'w': 3*N,
        'tx': 4*N,
        'ty': 5*N
        }
    
    #Delta = laplacian_coordinates.T.ravel()
    
    ## The following energy terms are minimized all at once, where the degrees of freedom
    ## are new_curve_i and T_i.
    ## The weight for E2 should be large, since it's a position constraint.
    ## I'm not sure about the weight for E3.
    ## E1 = \sum ( T_i * Delta_i - L_i*new_curve_i )^2
    ## E2 = \sum_{ j in constrained vertices } ( new_curve_j - constraints[j] )^2
    ## E3 = \sum ( T_i * curve_i - new_curve_i )^2 + \sum_{ j in neighbors(i) } ( T_i * curve_j - new_curve_j )^2
    
    system = zeros( ( 6*N, 6*N ) )
    rhs = zeros( 6*N )
    
    w_E1 = 1e1
    w_E2 = 1e7
    ## This needs to be large or we get a sharp kink at constrained points.
    w_E3 = 1e7
    
    ## Following the explicit translations of E1, E2, E3 into matrices found in my
    ## black notebook under "curve dragging - laplacian curve editing"
    
    ## E1
    A = zeros( ( 2*N, 6*N ) )
    for i in range(N):
        A[ 2*i + 0, o['a'] + i ] = laplacian_coordinates[i,0]
        A[ 2*i + 0, o['w'] + i ] = laplacian_coordinates[i,1]
        
        A[ 2*i + 0, o['tx'] + i ] = 1
        
        A[ 2*i + 0, o['x'] : o['x'] + N ] = -L[i,:]
        
        
        A[ 2*i + 1, o['a'] + i ] = laplacian_coordinates[i,1]
        A[ 2*i + 1, o['w'] + i ] = -laplacian_coordinates[i,0]
        
        A[ 2*i + 1, o['ty'] + i ] = 1
        
        A[ 2*i + 1, o['y'] : o['y'] + N ] = -L[i,:]
    
    system += w_E1 * dot( A.T, A )
    
    ## E2
    C = zeros( ( 2*len( constraints ), 6*N ) )
    Crhs = zeros( C.shape[0] )
    for i, (constraint_index, constraint_pos) in enumerate( constraints ):
        C[ 2*i + 0, o['x'] + constraint_index ] = 1
        C[ 2*i + 1, o['y'] + constraint_index ] = 1
        Crhs[ 2*i : 2*(i+1) ] = constraint_pos
    
    system += w_E2 * dot( C.T, C )
    rhs += w_E2 * dot( C.T, Crhs )
    
    ## E3
    R = zeros( ( 2*( N + d.sum() ), 6*N ) )
    def fill_rows_of_R_for_vertex_index( row0, T_index, v_index ):
        R[ row0 + 0, o['a'] + T_index ] = curve[v_index,0]
        R[ row0 + 0, o['w'] + T_index ] = curve[v_index,1]
        
        R[ row0 + 0, o['tx'] + T_index ] = 1
        
        R[ row0 + 0, o['x'] + v_index ] = -1
        
        
        R[ row0 + 1, o['a'] + T_index ] = curve[v_index,1]
        R[ row0 + 1, o['w'] + T_index ] = -curve[v_index,0]
        
        R[ row0 + 1, o['ty'] + T_index ] = 1
        
        R[ row0 + 1, o['y'] + v_index ] = -1
    
    row = 0
    for i in range(N):
        fill_rows_of_R_for_vertex_index( row, i, i )
        row += 2
        
        for j in neighbors[ i ]:
            assert i != j
            fill_rows_of_R_for_vertex_index( row, i, j )
            row += 2
    
    assert row == R.shape[0]
    
    system += w_E3 * dot( R.T, R )
    
    
    if kDebugging:
        ## Test: constrain T to be identity.
        ## Result: curve doesn't change when two points are constrained to not move.
        for i in range(N):
            system[ o['a'] + i, : ] = 0.
            system[ o['a'] + i, o['a'] + i ] = 1.
            rhs[ o['a'] + i ] = 1.
            
            system[ o['w'] + i, : ] = 0.
            system[ o['w'] + i, o['w'] + i ] = 1.
            rhs[ o['w'] + i ] = 0.
            
            system[ o['tx'] + i, : ] = 0.
            system[ o['tx'] + i, o['tx'] + i ] = 1.
            rhs[ o['tx'] + i ] = 0.
            
            system[ o['ty'] + i, : ] = 0.
            system[ o['ty'] + i, o['ty'] + i ] = 1.
            rhs[ o['ty'] + i ] = 0.
    
    if kDebugging:
        print("det( system ):", linalg.det( system ))
    
    x_bar = linalg.solve( system, rhs )
    new_curve = zeros( curve.shape, curve.dtype )
    new_curve[:,0] = x_bar[ o['x'] : o['x'] + N ]
    new_curve[:,1] = x_bar[ o['y'] : o['y'] + N ]
    
    if kDebugging:
        ## print out T matrices
        T = zeros( ( 3,3 ) )
        T[2,2] = 1.
        for i in range(N):
            a = x_bar[ o['a'] + i ]
            w = x_bar[ o['w'] + i ]
            tx = x_bar[ o['tx'] + i ]
            ty = x_bar[ o['ty'] + i ]
            
            T[0,0] = a
            T[1,1] = a
            T[0,1] = w
            T[1,0] = -w
            T[0,2] = tx
            T[1,2] = ty
            
            print("T_%d:" % (i,))
            print(T.round(3))
    
    #assert max( [ ( ( new_curve[ci] - asarray( cpos ) )**2 ).sum() for ci, cpos in constraints ] ) < 1e-5
    #print max( [ sqrt( ( ( new_curve[ci] - asarray( cpos ) )**2 ) ).sum() for ci, cpos in constraints ] )
    ## If I adjust the constraint weight such that constraints are loosely obeyed, then this
    ## guaranteed snapping will look dicontinuous.
    #for ci, cpos in constraints: new_curve[ci] = cpos
    
    return new_curve


def test_curve( curve, curve_is_closed ):
    print("curve:", curve)
    L = L_for_curve( curve, curve_is_closed = curve_is_closed )
    print("L_for_curve:", L)
    new_curve = new_curve_from_curve_and_constraints(
        curve,
        curve_is_closed = curve_is_closed,
        ## constrain the endpoints where they already are; new_curve should be the same as curve
        constraints = [ ( 0, curve[0] ), ( len(curve)-1, curve[-1] ) ]
        )
    print("new_curve_from_curve_and_constraints:", new_curve)
    print("| curve - new_curve |:", sqrt( ( ( curve - new_curve )**2 ).sum() ))

def test():
    ts = linspace( 0, 1, 10 )
    
    ## A horizontal line from x = 0...1.
    curve = list(zip( ts, 0*ts ))
    closed = False
    test_curve( curve, closed )
    
    print('====================')
    
    ## A circle.
    curve = list(zip( cos( 2*pi*ts ), sin( 2*pi*ts ) ))
    closed = True
    test_curve( curve, closed )

def main():
    test()

if __name__ == "__main__": main()
