## Why can't I do this anymore?
from numpy import array, asarray, linalg, int32, float32, float64, matrix, asmatrix, zeros, prod, sum, arange
#import pydb

def svd( int_arrays, float_arrays, double_arrays ):
    print('eureka')
    
    mat = array( double_arrays[0] ).reshape( int_arrays[0] )
    print('mat:', mat)
    
    u,w,vt = linalg.svd( mat )
    print('u:', u)
    print('w:', w)
    print('vt:', vt)

def reduced_rows( int_arrays, float_arrays, double_arrays ):
    #pydb.debugger()
    
    mat = array( double_arrays[0] ).reshape( int_arrays[0] )
    
    result_lm = get_reduced_row_lhs_xform( mat )
    result_prod = result_lm @ mat
    
    int_arrays[1][:2] = result_lm.shape
    double_arrays[1][ : prod( result_lm.shape ) ] = result_lm.flat
    
    int_arrays[2][:2] = result_prod.shape
    double_arrays[2][ : prod( result_prod.shape ) ] = result_prod.flat
    
    print('result_left_multiply:')
    print(result_lm)
    print('result_left_multiply @ input_matrix:')
    print(result_prod)

def get_reduced_row_lhs_xform( A ):
    u,w,vt = linalg.svd( A )
    
    #print 'w:', w
    
    eps = 1e-7
    ranks = w > eps
    rank = sum( ranks )
    wp = array( w )
    wp[ ranks ] = 1 / wp[ ranks ]
    
    Sinv = zeros( (rank, rank) )
    Sinv[ ( arange( rank ), arange( rank ) ) ] = wp[ :rank ]
    result_lm = Sinv @ u.T[ :rank, : ]
    
    #Sinv = zeros( int_arrays[0] )
    #Sinv[ ( arange( Sinv.shape[0] ), arange( Sinv.shape[1] ) ) ] = wp
    #result_lm = matrix( Sinv ) * matrix( u ).T
    
    return result_lm

def hard( int_arrays, float_arrays, double_arrays ):
    ia = int_arrays[0]
    ja = int_arrays[1]
    Acon_rows = int_arrays[2][0]
    A_rows = int_arrays[2][1]
    
    values = double_arrays[0]
    bcon = double_arrays[1]
    x = double_arrays[2]
    
    assert ia[ -1 ] == len( values )
    assert ia[ -1 ] == len( ja )
    
    from scipy import sparse
    from scipy.sparse import linalg as linsolve
    Acon = sparse.csr_matrix( (values, ja, ia) )
    # pydb.debugger()
    x[:] = linsolve.spsolve( Acon, bcon )[ : A_rows ]


last_Acon_shape = (-1,-1)
last_system = None
last_C = None
last_reduce_lm = None
last_factorized = None
def lsq2( int_arrays, float_arrays, double_arrays ):
    global last_Acon_shape, last_system, last_C, last_reduce_lm, last_factorized
    
    ia = int_arrays[0]
    ja = int_arrays[1]
    Acon_rows = int_arrays[2][0]
    A_rows = int_arrays[2][1]
    
    values = double_arrays[0]
    bcon = double_arrays[1]
    x = double_arrays[2]
    
    assert ia[ -1 ] == len( values )
    assert ia[ -1 ] == len( ja )
    
    from scipy import sparse
    from scipy.sparse import linalg as linsolve
    Acon = sparse.csr_matrix( (values, ja, ia) )
    
    if False: # Acon.shape == last_Acon_shape:
        system = last_system
        C = last_C
        reduce_lm = last_reduce_lm
        factorized = last_factorized
    else:
        ## NO LONGER TRUE for lil_matrix: sparse matrices don't support writable slices
        Acon = Acon.tolil()
        ## a dense matrix for slicing
        # Acon = Acon.todense()
        A = Acon[ :A_rows, :A_rows ]
        C = Acon[ A_rows:, :A_rows ]
        
        reduce_lm = get_reduced_row_lhs_xform( C.T @ C )
        #pydb.debugger()
        Cp = (reduce_lm @ C.T) @ C
        
        ## NO LONGER TRUE for lil_matrix: sparse matrices don't support writable slices
        system = sparse.lil_matrix( (A.shape[0] + Cp.shape[0], A.shape[0] + Cp.shape[0]) )
        ## a dense matrix for slicing
        # system = matrix( zeros( (A.shape[0] + Cp.shape[0], A.shape[0] + Cp.shape[0]) ) )
        system[ :A_rows, :A_rows ] = A
        system[ A_rows:, :A_rows ] = Cp
        system[ :A_rows, A_rows: ] = Cp.T
        system = sparse.csr_matrix( system ) #, dtype = 'F' )
        factorized = linsolve.factorized( system )
        
        last_Acon_shape = Acon.shape
        last_system = system
        last_C = C
        last_reduce_lm = reduce_lm
        last_factorized = factorized
    
    C_rhs = bcon[ A_rows: ]
    Cp_rhs = reduce_lm @ (C.T @ matrix( C_rhs ).T)
    
    sys_rhs = zeros( A_rows + Cp_rhs.shape[0] )
    sys_rhs[ :A_rows ] = bcon[ :A_rows ]
    sys_rhs[ A_rows: ] = Cp_rhs.flat
    
    #pydb.debugger()
    # linsolve.use_solver( use = {'useUmfpack' : True} )
    # x[:] = linsolve.spsolve( system, sys_rhs )[ : A_rows ]
    x[:] = factorized( sys_rhs )[ : A_rows ]


def set_row( A, row, indices, values ):
    from scipy import sparse
    A[ row, : ] = sparse.csr_matrix( ( values, array( list(zip( [0]*len(indices), indices )) ).T ), shape = (1,A.shape[1]), dtype = A.dtype )
    ## Scipy 0.8 (or maybe 0.7+)'s lil matrices no longer allow setting a slice from a single scalar value.
    ## UPDATE: Or even sometimes from a sequence.
    #A[ row, : ] = 0.
    #A[ row, : ] *= 0.
    #A[ row, tuple( indices ) ] = values
def insert_col( A, col, indices, values ):
    
    ## row major matrix representations can't handle column slicing
    for index, val in zip( indices, values ): A[ index, col ] = val
    #print A
    #print col
    #print indices
    #print values
    ## This fails for some reason
    #A[ tuple( indices ), col ] = values
    ## This silently fails
    #A.T[ col, tuple( indices ) ] = values

## TODO: Make a factor_system_with_constraints() function which return a functor that operates on a passed RHS.
def solve_system_with_constraints( A, b, pos_constraints = None, linear_constraints = None, linear_constraint_type = 'lsq', weight = 1000., debugging = False ):
    '''
    == Input parameters ==
    A is the basic, unconstrained system matrix (NxN).
    b is the N-vector right-hand-side of the system.
    pos_constraints is a list of (index, value) tuples,
        where the degree-of-freedom at 'index' is constrained to have value 'value'.
    linear_constraints is a list of (indices, values, rhs) tuples,
        where 'indices' is a list of integer degree-of-freedom indices,
        'values' is a list of coefficients for the indices, and 'rhs' is a scalar value for the rhs
        of the equation.  That is, \sum values[i] * dofs[ indices[i] ] = rhs.
    
    linear_constraint_type (optional) takes one of 3 values:
        if 'lsq', then 'linear_constraints' are treated as least squares constraints.
        if 'hard', then 'linear_constraints' are treated as hard constraints (lagrange multipliers).
        if 'hard_lsq', then 'linear_constraints' are treated as hard constraints that must be
            satisfied in a least squares sense (if a solution exists, then they are solved exactly).
    weight (optional) is the weight on the least squares constraints if 'linear_constraint_type' is 'lsq'.
    
    NOTE: If the linear_constraint_type is 'lsq', then positional constraints trump linear constraints.
          If the linear_constraint_type is 'hard' or 'hard_lsq', then linear constraints trump positional constraints.
    
    == Output parameters ==
    Returns an N-vector of values for the dofs.
    '''
    
    from scipy import sparse
    from scipy.sparse import linalg as linsolve
    from numpy import where
    
    
    Real = b.dtype
    assert A.dtype == b.dtype
    
    assert A.shape[0] == A.shape[1]
    assert b.shape[0] == A.shape[1]
    
    if pos_constraints is None: pos_constraints = []
    if linear_constraints is None: linear_constraints = []
    
    if linear_constraint_type == 'lsq':
        
        #import pydb
        #pydb.debugger()
        
        if linear_constraints:
            ## Set-up linear constraints matrix
            C = sparse.lil_matrix( ( len( linear_constraints ), A.shape[0] ), dtype = Real )
            C_rhs = zeros( len( linear_constraints ), dtype = Real )
            for i, (indices, values, rhs) in enumerate( linear_constraints ):
                set_row( C, i, indices, values )
                C_rhs[ i ] = rhs
            
            if debugging:
                print("C:")
                print(C)
                print("C.T * C:")
                print(C.T * C)
            
            system = A + weight*( C.T * C )
            rhs = b + weight*( C.T * C_rhs )
            
        else:
            
            system = sparse.lil_matrix( A, copy = True )
            rhs = array( b )
        
        if pos_constraints:
            ## Set hard position constraints after we setup the least squares system,
            ## so that they're not affected by the least squares constraints.
            system = sparse.lil_matrix( system )
            for index, value in pos_constraints:
                ## Set the corresponding row to the identity matrix, and set the right-hand-side to the desired value.
                set_row( system, index, [index], [1.] )
                rhs[ index ] = value
        
        ## For solver efficiency, convert to csr()
        system = sparse.csr_matrix( system )
        
        if debugging:
            print('b', b.round(5))
            print('b constrained', rhs.round(5))
            print('A:')
            print(A)
            print('A constrained:')
            print(system)
        
        # factorized_system =linsolve.factorized( system )
        # x = factorized_system( rhs )
        x = linsolve.spsolve( system, rhs )
        if debugging: print('solution:', x.round(5))
    
    elif linear_constraint_type == 'hard':
        
        ## Create the big, constrained matrix
        Acon = sparse.lil_matrix( ( A.shape[0] + len(linear_constraints), A.shape[0] + len(linear_constraints) ), dtype = Real )
        ## copy A to the upper-left of Acon
        Acon[ :A.shape[0], :A.shape[1] ] = A
        ## copy b to bcon
        bcon = zeros( len(b) + len(linear_constraints) )
        bcon[ : len(b) ] = b
        
        ## Set hard position constraints
        for index, value in pos_constraints:
            ## Set the corresponding row to the identity matrix, and set the right-hand-side to the desired value.
            set_row( Acon, index, [index], [1.] )
            bcon[ index ] = value
        
        ## Set-up linear constraint
        for i, (indices, values, rhs) in enumerate( linear_constraints ):
            ## Set the row and column
            set_row( Acon, A.shape[0] + i, indices, values )
            insert_col( Acon, A.shape[0] + i, indices, values )
            bcon[ A.shape[0] + i ] = rhs
        
        ## For solver efficiency, convert to csr()
        Acon = Acon.tocsr()
        
        if debugging:
            print('b', b.round(5))
            print('b constrained', bcon.round(5))
            print('A:')
            print(A)
            print('A constrained:')
            print(Acon)
        
        x = linsolve.spsolve( Acon, bcon )[ : len(b) ]
        if debugging: print('solution:', x.round(5))
    
    elif linear_constraint_type == 'hard_lsq':
        
        lcp = linear_constraints
        
        if linear_constraints:
            
            ## Set-up linear constraints matrix
            C = sparse.lil_matrix( ( len( linear_constraints ), A.shape[0] ), dtype = Real )
            C_rhs = zeros( len( linear_constraints ), dtype = Real )
            for i, (indices, values, rhs) in enumerate( linear_constraints ):
                set_row( C, i, indices, values )
                C_rhs[ i ] = rhs
            
            #import pydb
            #pydb.debugger()
            ## TODO: a sparse version of this
            reduce_lm = sparse.csr_matrix( get_reduced_row_lhs_xform( ( C.T @ C ).todense() ) )
            Cp = (reduce_lm @ C.T) @ C
            Cp_rhs = (reduce_lm @ C.T) @ C_rhs
            
            ## eqs is a dictionary mapping rows to (list of indices, list of coefficients for indices, rhs)
            eqs = {}
            Cp = sparse.coo_matrix( Cp )
            for row, col, val in zip( Cp.row, Cp.col, Cp.data ):
                
                indices, coefs, rhs = eqs.setdefault( row, ([], [], Cp_rhs[ row ]) )
                indices.append( col )
                coefs.append( val )
            
            #from pprint import pprint
            #pprint( eqs )
            
            lcp = list(eqs.values())
        
        return solve_system_with_constraints( A, b, pos_constraints, lcp, 'hard' )
        
    else:
        raise RuntimeError('Unknown linear constraint type: ' + linear_constraint_type)
    
    return x


def main():
    #test1()
    test2()

def test1():
    svd( [ array( [3, 3], int32 ) ], [], [ array( [1, 0, 0, 0, 1, 0, 0, 0, 1], float64 ) ] )

def test2():
    ## Test solve_system_with_constraints()
    from scipy import sparse
    from numpy import zeros, ones
    import mymath
    # to easily reload
    import importlib
    importlib.reload( mymath )
    A = sparse.eye( 5,5 ).tocsr()
    b = ones( A.shape[0] )
    contypes = ['lsq', 'hard', 'hard_lsq']
    for contype in contypes: print(mymath.solve_system_with_constraints( A, b, [], [], contype ).tolist(), contype)
    print('===')
    for contype in contypes: print(mymath.solve_system_with_constraints( A, b, [(3,31337.)], [], contype ).tolist(), contype)
    print('===')
    for contype in contypes: print(mymath.solve_system_with_constraints( A, b, [], [ ([3], [1.], 31337.) ], contype ).tolist(), contype)
    print('===')
    for contype in contypes: print(mymath.solve_system_with_constraints( A, b, [], [ ([3], [2.], 31337.) ], contype ).tolist(), contype)
    print('===')
    for contype in contypes: print(mymath.solve_system_with_constraints( A, b, [(3,8008.)], [ ([3], [2.], 31337.) ], contype ).tolist(), contype)
    print('===')
    for contype in contypes: print(mymath.solve_system_with_constraints( A, b, [(4,8008.)], [ ([3], [2.], 31337.) ], contype ).tolist(), contype)

if __name__ == '__main__': main()
