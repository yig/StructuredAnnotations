from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *
from vecutil import *
#from helpers import *

from pprint import pprint

## If kEndpoints is True, then energy will be calculated at the end-points of the lines.
## This is usually a good thing to do, except that something is wrong with my end-point
## energy, so I'm leaving it off.
## The (only) symptom I notice is that with only one tilted cross section constraint, the tilt
## doesn't propogate along the spine (it is damped).
## I think then that energy is leaking off of the ends of the spine.
kEndpoints = False

def optimize_sweep_spine2d( spine2d, constraints ):
    '''
    Given a 2d chain of vertices (sequence of 2d points) and a list of constraints (pairs of
    the vertex index to constrain and the 2d point to constrain it to),
    returns a beautified chain of vertices.
    '''
    
    from springs_backwards import static_forwards_solution, Environment
    import Eline_ad
    #Eline_ad.Ecurvature_and_derivatives = Eline_ad.Emag_and_derivatives
    #Eline_ad.Ecurvature_and_derivatives = Eline_ad.Elaplacian_editing_and_derivatives
    
    env = Environment()
    env.DisableAllForces()
    #env.kBending = True
    #env.b = 1000.
    #env.kCurvature = True
    #env.kCurvatureUndeformed = True
    #env.c = 1.
    #env.kLaplacianEditing = True
    env.kLaplacianEditingUndeformed = True
    env.l = 1.
    #env.kGravity = True
    #env.g = 10.
    
    N = len(spine2d)
    edges = list(zip( list(range(N-1)), list(range(1,N)) ))
    
    #pprint( "spine2d:" )
    #pprint( spine2d )
    newpos, converged = static_forwards_solution(
        spine2d, edges, spine2d, [], constraints, env,
        verbose = False
        )
    if converged:
        #pprint( "newpos:" )
        #pprint( spine2d )
        return newpos
    else:
        return spine2d

def optimize_sweep_spine3d( spine3d, constraints ):
    '''
    Given a 3d chain of vertices (sequence of 3d points) and a list of constraints (pairs of
    the vertex index to constrain and the 3d point to constrain it to),
    returns a beautified chain of vertices.
    
    used
    '''
    
    ## Spine must be non-empty.
    assert len( spine3d ) > 0
    ## Constraints must be non-empty.
    assert len( constraints ) > 0
    ## Constraint indices must be valid (in [0,len(spine3d)).
    assert False not in [ index >= 0 and index < len( spine3d ) for index, val in constraints ]
    
    ## Bake in the constraints as hard constraints, so we can set the corresponding gradients
    ## to zero.
    spine3d = array( spine3d )
    for constraint_vertex_index, constraint_point in constraints:
        spine3d[constraint_vertex_index] = constraint_point
    
    def constrain_grad_hess( grad, hess ):
        ## Apply constraints, by setting constrained degrees-of-freedom to have rows in the hessian
        ## from the identity matrix and gradients of 0.
        ## NOTE: Settings gradients to zero assumes that the current values already have the
        ##       constraints baked in, as hard constraints.
        for constraint_vertex_index, constraint_point in constraints:
            dim = 3
            for c in range(dim):
                if hess is not None:
                    hess[ constraint_vertex_index*dim + c, : ] = 0.
                    hess[ constraint_vertex_index*dim + c, constraint_vertex_index*dim + c ] = 1.
            grad[ constraint_vertex_index*dim : ( constraint_vertex_index + 1 )*dim ] = 0.
    
    return optimize_sweep_spine3d_constraint_function( spine3d, constrain_grad_hess )


def optimize_sweep_spine3d_constraint_function( spine3d, constrain_grad_hess ):
    '''
    Given a 3d chain of vertices (sequence of 3d points) 'spine3d' and
    a function 'constrain_grad_hess' to apply constraints to the gradient and hessian,
    returns a beautified chain of vertices.
    Raises a RuntimeError if the optimization procedure diverges.
    
    The function 'constrain_grad_hess' takes two parameters which it is expected to modify:
        a 3N array corresponding to the gradient, and
        a 3N-by-3N matrix corresponding to the hessian.
    NOTE: The hessian parameter may be None, in which case only the gradient should be modified.
    
    used
    '''
    
    return optimize_shell_with_energy_and_constraint_functions(
        spine3d,
        [ (1., Elaplacian_editing_undeformed( spine3d )) ],
        [ constrain_grad_hess ]
        )

def Elaplacian_editing_undeformed( spine3d ):
    '''
    Given an undeformed 3d chain of vertices (sequence of 3d points) 'spine3d',
    returns a function suitable for passing to optimize_shell_with_energy_and_constraint_functions()
    that minimizes the laplacian along the chain.
    
    used
    '''
    
    N = len( spine3d )
    
    ## For the line energy
    import Eline_ad
    
    def e_grad_hess( xvec, compute_hessian ):
        ## Reshape the values into a 3d array.
        xvec = xvec.reshape( (N,3) )
        
        ## Initialize the energy, gradient, and hessian.
        ## We need to return a scalar, a 3N array, and a 3N x 3N array.
        E = 0.
        grad = zerosf( 3*N )
        if compute_hessian:
            hess = zerosf( (3*N,3*N) )
        
        for i in range( 1, N-1 ):
            egh_i = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                xvec[ i-1 ], xvec[ i ], xvec[ i+1 ],
                spine3d[ i-1 ], spine3d[ i ], spine3d[ i+1 ],
                compute_hessian = compute_hessian
                )
            
            E += egh_i[0]
            grad[ (i-1)*3 : (i+2)*3 ] += egh_i[1]
            if compute_hessian:
                hess[ (i-1)*3:(i+2)*3, (i-1)*3:(i+2)*3 ] += egh_i[2]
        
        
        if kEndpoints:
            ## i == 0
            i = 0
            egh_i = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                None, xvec[ i ], xvec[ i+1 ],
                None, spine3d[ i ], spine3d[ i+1 ],
                compute_hessian = compute_hessian
                )
            E += egh_i[0]
            grad[ (i)*3 : (i+2)*3 ] += egh_i[1][3:]
            if compute_hessian:
                hess[ (i)*3:(i+2)*3, (i)*3:(i+2)*3 ] += egh_i[2][3:,3:]
            
            ## i == N-1
            i = N-1
            egh_i = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                xvec[ i-1 ], xvec[ i ], None,
                spine3d[ i-1 ], spine3d[ i ], None,
                compute_hessian = compute_hessian
                )
            E += egh_i[0]
            grad[ (i-1)*3 : (i+1)*3 ] += egh_i[1][:-3]
            if compute_hessian:
                hess[ (i-1)*3:(i+1)*3, (i-1)*3:(i+1)*3 ] += egh_i[2][:-3,:-3]
        
        if compute_hessian:
            return E, grad, hess
        else:
            return E, grad
    
    return e_grad_hess

def Elinear_lsq_constraints( eqs ):
    '''
    Given a sequence 'eqs' of equations, where each equation contains two items,
    a sequence of ( index, coefficient ) tuples and a right-hand side scalar, where
        \sum coefficient * x[index] = right-hand side scalar
    returns a function suitable for passing to optimize_shell_with_energy_and_constraint_functions()
    that minimizes the squared difference of each equation's left-hand side from its right-hand side.
    '''
    
    M = len( eqs )
    
    def e_grad_hess( xvec, compute_hessian ):
        N = len( xvec )
        
        ## Let's build an M-by-N matrix of all M equations' left-hand side coefficients,
        ## and an M-vector of all M equations' right-hand side values.
        C = zerosf( (M,N) )
        rhss = zerosf( M )
        
        for i, (lhs, rhs) in enumerate( eqs ):
            rhss[i] = rhs
            for index, coeff in lhs:
                C[ i, index ] = coeff
        
        ## This expression is used twice, so store it.
        diff = dot( C, xvec ) - rhss
        
        E = dot( diff, diff )
        grad = dot( diff, C )
        if compute_hessian:
            hess = dot( C.T, C )
            return E, grad, hess
        else:
            return E, grad
    
    return e_grad_hess

def optimize_shell_with_energy_and_constraint_functions( x0, energy_and_derivatives_funcs, constraint_funcs, linear = False ):
    '''
    Given an initial guess configuration 'x0' of N values in any shape,
    a sequence of one or more ( scalar weight, function ) pairs 'energy_and_derivatives_funcs' to compute the energy and its derivatives,
    and a sequence of one or more functions 'constrain_grad_hess' to apply constraints to the gradient and hessian,
    returns an optimized configuration in the original shape of 'x0'.
    Raises a RuntimeError if the optimization procedure diverges.
    
    The functions in 'energy_and_derivatives_funcs' take two parameters,
        the configuration 'x' to evaluate ('x' is a passed as a flat N array)
        and
        a boolean 'compute_hessian' determining whether or not to compute the hessian
    and return
        two items if compute_hessian is False: the energy as a scalar and the gradient as an N array,
        or
        three items if compute_hessian is True: the energy as a scalar, the gradient as an N array, and the hessian as an N-by-N array.
    
    The functions in 'constraint_funcs' takes two parameters which they are expected to modify:
        an N array corresponding to the gradient, and
        an N-by-N matrix corresponding to the hessian.
    NOTE: The hessian parameter may be None, in which case only the gradient should be modified.
    
    If optional parameter 'linear' is True, then the system is assumed to be linear and solved directly.
    
    used
    '''
    
    assert len( energy_and_derivatives_funcs ) > 0
    assert len( constraint_funcs ) > 0
    
    ## 1 Hessian, gradient callback shell function.
    ## 2 Call one of the scipy.optimize functions:
    ##      http://www.scipy.org/doc/api_docs/SciPy.optimize.html
    ##      http://www.scipy.org/doc/api_docs/SciPy.optimize.optimize.html#fmin_ncg
    
    x0 = asarray( x0 )
    N = prod( x0.shape )
    
    ## 1
    def e_grad_hess( xvec, compute_hessian ):
        ## Initialize the energy, gradient, and hessian.
        ## We need to return a scalar, a 3N array, and a 3N x 3N array.
        E = 0.
        grad = zerosf( N )
        if compute_hessian:
            hess = zerosf( (N,N) )
        else:
            hess = None
        
        ## Sum the energy, gradient, and derivatives.
        for weight, energy_and_derivative_func in energy_and_derivatives_funcs:
            energy_e_grad_hess = energy_and_derivative_func( xvec, compute_hessian )
            E += weight * energy_e_grad_hess[0]
            grad += weight * energy_e_grad_hess[1]
            if compute_hessian:
                hess += weight * energy_e_grad_hess[2]
        
        ## Apply the constraints.
        for constraint_func in constraint_funcs: constraint_func( grad, hess )
        if compute_hessian:
            return E, grad, hess
        else:
            return E, grad
    
    
    ## 2
    if linear:
        ## Actually, no need for optimize, this is a linear system.
        b, a = e_grad_hess( x0.ravel(), compute_hessian = True )[1:]
        return x0 - linalg.solve( a, b ).reshape( x0.shape )
    
    
    ## For optimize
    import scipy.optimize
    
    def callback( xvec ): print("iteration complete")
    #opt_result = scipy.optimize.fmin_cg( 
    #opt_result = scipy.optimize.fmin(
    opt_result = scipy.optimize.fmin_ncg(
        ## f
        lambda xvec: e_grad_hess( xvec, compute_hessian = False )[0],
        ## x0
        x0.ravel(),
        fprime = lambda xvec: e_grad_hess( xvec, compute_hessian = False )[1],
        fhess = lambda xvec: e_grad_hess( xvec, compute_hessian = True )[2],
        ## Don't output 'Optimization terminated successfully' messages and statistics.
        #disp = 0,
        ## The default relative tolerance is 1e-5, but it's too large and
        ## optimization often terminates early.
        avextol = 1e-10,
        #callback = callback,
        ## We need full output to determine if it terminated successfully.
        full_output = True
        )
    #pprint( "opt_result" )
    #pprint( opt_result )
    ## Full output item 5 is the termination error code; if it's zero, optimization converged.
    converged = not opt_result[5]
    
    if converged:
        newpos = opt_result[0].reshape( x0.shape )
        #pprint( "new positions:" )
        #pprint( newpos )
        return newpos
    else:
        pprint( "diverged!" )
        raise RuntimeError("optimize_shell diverged!")
        return x0

def diffuse_values_along_sweep_spine( spine, value_constraints ):
    '''
    Given a 1 or 2 or 3-dimensional chain of vertices (sequence of 1 or 2 or 3-dimensional points)
    and a list of value constraints (pairs of the vertex index to constrain and the
    1 or 2 or 3-dimensional value to constrain it to),
    returns an array of values, one for each point, where the values at unconstrained vertices
    smoothly interpolates the constrained values.
    
    NOTE: 1 or 2 or 3d values means each value is a sequence of that many scalars.
          In particular, this means that using 1d values means that 1-element sequences are used,
          and the result will be an N-by-1 array, not an N-vector.
          Call reshape() or ravel() on the resulting N-by-1 array to turn it into an N-vector.
    '''
    return diffuse_values_along_sweep_spine_factored( spine )( value_constraints )
    
    ## Spine must be non-empty.
    assert len( spine ) > 0
    ## Constraints must be non-empty.
    assert len( value_constraints ) > 0
    ## Constraint indices must be valid (in [0,len(spine)).
    assert False not in [ index >= 0 and index < len( spine ) for index, val in value_constraints ]
    ## Constraint values must all have the same dimension.
    assert 1 == len( set( [ len( val ) for index, val in value_constraints ] ) )
    
    ## If there is only one value in value_constraints, we don't have enough constraints for
    ## our laplacian-based energy.
    ## In this case, let all values be identical.
    if len( value_constraints ) == 1:
        N = len( spine )
        K = len( value_constraints[0][1] )
        values = emptyf( (N,K) )
        values[:,:] = value_constraints[0][1]
        return values
    
    ## We will optimize for an array of 1, 2, or 3-dimensional values using laplacian editing energy.
    ## We will always pass the spine as the undeformed positions --- laplacian editing energy is
    ## only going to look at the lengths of undeformed edges.
    ## Our laplacian editing energy functions always take 3d points, which is where our limitation
    ## of dimensions <= 3 comes from.
    
    ## 1 Make a 3d array of undeformed positions from the spine.
    ## 2 Make an array of values to solve for, initializing them to the constraints and zero elsewhere.
    ## 3 Write a callback function similar to the one in optimize_sweep_spine3d(), but only returning
    ##   gradient and hessian for the k dimensions of the values we care about.
    
    ## For optimize
    import scipy.optimize
    ## For the line energy
    import Eline_ad
    
    ## 1 Make a 3d array of undeformed positions from the spine.
    N = len(spine)
    spine3d = zerosf( (N,3) )
    spine3d[ :, :len( spine[0] ) ] = spine
    
    ## 2
    K = len( value_constraints[0][1] )
    values = zerosf( (N,K) )
    ## Bake in the constraints as hard constraints, so we can set the corresponding gradients
    ## to zero.
    for constraint_vertex_index, constraint_point in value_constraints:
        values[constraint_vertex_index][:K] = constraint_point
    
    
    ## 1
    def e_grad_hess( xvec_in, compute_hessian ):
        ## Reshape the vector into an array of 3d points.
        xvec = zerosf( (N,3) )
        xvec[:,:K] = xvec_in.reshape( (N,K) )
        
        ## Initialize the energy, gradient, and hessian.
        ## We need to return a scalar, a 3N array, and a 3N x 3N array.
        E = 0.
        grad = zerosf( K*N )
        if compute_hessian:
            hess = zerosf( (K*N,K*N) )
        
        for i in range( 1, N-1 ):
            egh_i = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                xvec[ i-1 ], xvec[ i ], xvec[ i+1 ],
                spine3d[ i-1 ], spine3d[ i ], spine3d[ i+1 ],
                compute_hessian = compute_hessian
                )
            
            E += egh_i[0]
            for c in range(K):
                grad[ c + (i-1)*K : (i+2)*K : K ] += egh_i[1][c::3]
            if compute_hessian:
                for ci in range(K):
                    for cj in range(K):
                        hess[ ci + (i-1)*K : (i+2)*K : K, cj + (i-1)*K : (i+2)*K : K ] += egh_i[2][ci::3,cj::3]
        
        
        if kEndpoints:
            ## i == 0
            i = 0
            egh_i = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                None, xvec[ i ], xvec[ i+1 ],
                None, spine3d[ i ], spine3d[ i+1 ],
                compute_hessian = compute_hessian
                )
            E += egh_i[0]
            for c in range(K):
                grad[ c + (i)*K : (i+2)*K : K ] += egh_i[1][c+3::3]
            if compute_hessian:
                for ci in range(K):
                    for cj in range(K):
                        hess[ ci + (i)*K : (i+2)*K : K, cj + (i)*K : (i+2)*K : K ] += egh_i[2][ci+3::3,cj+3::3]
            
            ## i == N-1
            i = N-1
            egh_i = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                xvec[ i-1 ], xvec[ i ], None,
                spine3d[ i-1 ], spine3d[ i ], None,
                compute_hessian = compute_hessian
                )
            E += egh_i[0]
            for c in range(K):
                grad[ c + (i-1)*K : (i+1)*K : K ] += egh_i[1][c:-3:3]
            if compute_hessian:
                for ci in range(K):
                    for cj in range(K):
                        hess[ ci + (i-1)*K : (i+1)*K : K, cj + (i-1)*K : (i+1)*K : K ] += egh_i[2][ci:-3:3,cj:-3:3]
        
        ## Apply constraints, by setting constrained degrees-of-freedom to have rows in the hessian
        ## from the identity matrix and gradients of 0.
        ## NOTE: Settings gradients to zero assumes that the current values (xvec) already have the
        ##       constraints baked in, as hard constraints.
        for constraint_vertex_index, constraint_point in value_constraints:
            dim = K
            for c in range(dim):
                if compute_hessian:
                    hess[ constraint_vertex_index*dim + c, : ] = 0.
                    hess[ constraint_vertex_index*dim + c, constraint_vertex_index*dim + c ] = 1.
            grad[ constraint_vertex_index*dim : ( constraint_vertex_index + 1 )*dim ] = 0.
        
        if compute_hessian:
            return E, grad, hess
        else:
            return E, grad
    
    
    ## 2
    def callback( xvec ): print("iteration complete")
    #opt_result = scipy.optimize.fmin_cg( 
    #opt_result = scipy.optimize.fmin(
    opt_result = scipy.optimize.fmin_ncg(
        ## f
        lambda xvec: e_grad_hess( xvec, compute_hessian = False )[0],
        ## x0
        values.ravel(),
        fprime = lambda xvec: e_grad_hess( xvec, compute_hessian = False )[1],
        fhess = lambda xvec: e_grad_hess( xvec, compute_hessian = True )[2],
        ## Don't output 'Optimization terminated successfully' messages and statistics.
        disp = 0,
        ## The default relative tolerance is 1e-5, but I changed it in
        ## optimize_shell_with_energy_and_constraint_functions() so I'm also changing it here.
        avextol = 1e-10,
        #callback = callback,
        ## We need full output to determine if it terminated successfully.
        full_output = True
        )
    #pprint( "opt_result" )
    #pprint( opt_result )
    ## Full output item 5 is the termination error code; if it's zero, optimization converged.
    converged = not opt_result[5]
    
    if converged:
        new_values = opt_result[0].reshape( (N,K) )
        #pprint( "values:" )
        #pprint( new_values )
        return new_values
    else:
        pprint( "diverged!" )
        return values

def diffuse_values_along_sweep_spine_factored( spine ):
    '''
    Given a 1 or 2 or 3-dimensional chain of vertices (sequence of 1 or 2 or 3-dimensional points),
    returns an object which can be called with 'value_constraints' as in
    diffuse_values_along_sweep_spine().
    This is useful for repeatedly diffusing values along the same spine.
    Also, 'value_constraints' may have any dimension (unlike diffuse_values_along_sweep_spine()).
    '''
    
    ## Spine must be non-empty.
    assert len( spine ) > 0
    
    ## For the line energy
    import Eline_ad
    
    ## 1 Make a 3d array of undeformed positions from the spine.
    N = len(spine)
    spine3d = zerosf( (N,3) )
    spine3d[ :, :len( spine[0] ) ] = spine
    
    K = 1
    def e_grad_hess( xvec, compute_hessian ):
        ## Initialize the energy, gradient, and hessian.
        ## We need to return a scalar, a 3N array, and a 3N x 3N array.
        E = 0.
        grad = zerosf( K*N )
        if compute_hessian:
            hess = zerosf( (K*N,K*N) )
        
        for i in range( 1, N-1 ):
            egh_i = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                xvec[ i-1 ], xvec[ i ], xvec[ i+1 ],
                spine3d[ i-1 ], spine3d[ i ], spine3d[ i+1 ],
                compute_hessian = compute_hessian
                )
            
            E += egh_i[0]
            for c in range(K):
                grad[ c + (i-1)*K : (i+2)*K : K ] += egh_i[1][c::3]
            if compute_hessian:
                for ci in range(K):
                    for cj in range(K):
                        hess[ ci + (i-1)*K : (i+2)*K : K, cj + (i-1)*K : (i+2)*K : K ] += egh_i[2][ci::3,cj::3]
        
        
        if kEndpoints:
            ## i == 0
            i = 0
            egh_i = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                None, xvec[ i ], xvec[ i+1 ],
                None, spine3d[ i ], spine3d[ i+1 ],
                compute_hessian = compute_hessian
                )
            E += egh_i[0]
            for c in range(K):
                grad[ c + (i)*K : (i+2)*K : K ] += egh_i[1][c+3::3]
            if compute_hessian:
                for ci in range(K):
                    for cj in range(K):
                        hess[ ci + (i)*K : (i+2)*K : K, cj + (i)*K : (i+2)*K : K ] += egh_i[2][ci+3::3,cj+3::3]
            
            ## i == N-1
            i = N-1
            egh_i = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                xvec[ i-1 ], xvec[ i ], None,
                spine3d[ i-1 ], spine3d[ i ], None,
                compute_hessian = compute_hessian
                )
            E += egh_i[0]
            for c in range(K):
                grad[ c + (i-1)*K : (i+1)*K : K ] += egh_i[1][c:-3:3]
            if compute_hessian:
                for ci in range(K):
                    for cj in range(K):
                        hess[ ci + (i-1)*K : (i+1)*K : K, cj + (i-1)*K : (i+1)*K : K ] += egh_i[2][ci:-3:3,cj:-3:3]
        
        if compute_hessian:
            return E, grad, hess
        else:
            return E, grad
    
    hess = e_grad_hess( spine3d, True )[2]
    
    class DiffuseValuesAlongSpineFactored( object ):
        def __init__( self, hessian ):
            ## hessian should be a square NxN matrix.
            #hessian = asarray( hessian )
            assert len( hessian.shape ) == 2
            assert hessian.shape[0] == hessian.shape[1]
            
            self.hessian = hessian
        
        def __call__( self, value_constraints ):
            N = self.hessian.shape[0]
            
            ## Constraints must be non-empty.
            assert len( value_constraints ) > 0
            ## Constraint indices must be valid (in [0,N)).
            assert False not in [ index >= 0 and index < N for index, val in value_constraints ]
            ## Constraint values must all have the same dimension.
            assert 1 == len( set( [ len( val ) for index, val in value_constraints ] ) )
            
            K = len( value_constraints[0][1] )
            
            A = array( self.hessian )
            b = zerosf( (N,K) )
            
            if len( value_constraints ) == 1:
                b[:,:] = value_constraints[0][1]
                return b
            
            for index, val in value_constraints:
                A[ index, : ] = 0.
                A[ index, index ] = 1.
                
                b[ index ] = val
            
            return linalg.solve( A, b )
    
    return DiffuseValuesAlongSpineFactored( hess )

def test2d():
    N = 21
    sweep_spine2d = array( [ (x, 0.) for x in linspace( 0, 1, N ) ] )
    dx = abs( sweep_spine2d[1][0] - sweep_spine2d[0][0] )
    
    test_type = 2
    if test_type == 1:
        ## Recovery test; perturb a value in the middle.
        sweep_spine2d[N//2][1] = 0.01
        #constraints = [ (0, sweep_spine2d[0]), (N-1, sweep_spine2d[N-1]) ]
        constraints = [
            (0, sweep_spine2d[0]),
            (1, sweep_spine2d[1]),
            (N-2, sweep_spine2d[N-2]),
            (N-1, sweep_spine2d[N-1])
            ]
    
    elif test_type == 2:
        constraints = {}
        
        assert N//2 > 0 and N//2 + 1 < N-1 and N//2 - 1 > 0
        
        ## Clamp one vertex at each end.
        constraints[0] = sweep_spine2d[0].copy()
        constraints[N-1] = sweep_spine2d[N-1].copy()
        ## Clamp the middle 3 vertices in a vertical line.
        constraints[N//2] = sweep_spine2d[N//2].copy()
        constraints[ N//2 - 1 ] = sweep_spine2d[ N//2 ].copy()
        constraints[ N//2 - 1 ][1] = -dx
        constraints[ N//2 + 1 ] = sweep_spine2d[ N//2 ].copy()
        constraints[ N//2 + 1 ][1] = dx
        
        constraints = list(constraints.items())
    
    elif test_type == 3:
        constraints = {}
        
        assert N//2 > 1 and N//2 + 1 < N-2
        
        ## Clamp two vertices at each end, vertically.
        constraints[0] = sweep_spine2d[0].copy()
        constraints[1] = sweep_spine2d[0].copy()
        constraints[1][1] = dx
        constraints[N-1] = sweep_spine2d[N-1].copy()
        constraints[N-2] = sweep_spine2d[N-1].copy()
        constraints[N-2][1] = dx
        
        ## NOTE: This isn't symmetric; I'd have to set N//2's y to -dx (and have an even # of points)
        constraints[N//2] = sweep_spine2d[ N//2 ].copy()
        constraints[ N//2 + 1 ] = sweep_spine2d[ N//2 ].copy()
        constraints[ N//2 + 1 ][1] = dx
        
        constraints = list(constraints.items())
    
    else:
        raise RuntimeError("Unknown test type")
    
    
    print("sweep_spine2d:")
    pprint( sweep_spine2d )
    
    print("constraints:")
    pprint( constraints )
    
    opt_spine = optimize_sweep_spine2d( sweep_spine2d, constraints )
    
    print("opt_spine:")
    pprint( opt_spine )
    
    test_opt_vals( opt_spine )
    
    display_line( opt_spine, update_line_func = lambda line: optimize_sweep_spine2d( line, constraints ) )

def test3d():
    N = 21
    sweep_spine3d = array( [ (x, 0., 0.) for x in linspace( 0, 1, N ) ] )
    dx = abs( sweep_spine3d[1][0] - sweep_spine3d[0][0] )
    
    test_type = 4
    if test_type == 1:
        ## Recovery test; perturb a value in the middle.
        sweep_spine3d[N//2][1] = 0.01
        #constraints = [ (0, sweep_spine2d[0]), (N-1, sweep_spine2d[N-1]) ]
        constraints = [
            (0, sweep_spine3d[0]),
            (1, sweep_spine3d[1]),
            (N-2, sweep_spine3d[N-2]),
            (N-1, sweep_spine3d[N-1])
            ]
    
    elif test_type == 2:
        constraints = {}
        
        assert N//2 > 0 and N//2 + 1 < N-1 and N//2 - 1 > 0
        
        ## Clamp one vertex at each end.
        constraints[0] = sweep_spine3d[0].copy()
        constraints[N-1] = sweep_spine3d[N-1].copy()
        ## Clamp the middle 3 vertices in a vertical line.
        constraints[N//2] = sweep_spine3d[N//2].copy()
        constraints[ N//2 - 1 ] = sweep_spine3d[ N//2 ].copy()
        constraints[ N//2 - 1 ][1] = -dx
        constraints[ N//2 + 1 ] = sweep_spine3d[ N//2 ].copy()
        constraints[ N//2 + 1 ][1] = dx
        
        constraints = list(constraints.items())
    
    elif test_type == 3:
        constraints = {}
        
        assert N//2 > 1 and N//2 + 1 < N-2
        
        ## Clamp two vertices at each end, vertically.
        constraints[0] = sweep_spine3d[0].copy()
        constraints[1] = sweep_spine3d[0].copy()
        constraints[1][1] = dx
        constraints[N-1] = sweep_spine3d[N-1].copy()
        constraints[N-2] = sweep_spine3d[N-1].copy()
        constraints[N-2][1] = dx
        
        ## NOTE: This isn't symmetric; I'd have to set N//2's y to -dx (and have an even # of points)
        constraints[N//2] = sweep_spine3d[ N//2 ].copy()
        constraints[ N//2 + 1 ] = sweep_spine3d[ N//2 ].copy()
        constraints[ N//2 + 1 ][1] = dx
        
        constraints = list(constraints.items())
    
    elif test_type == 4:
        constraints = {}
        
        assert N//2 > 0 and N//2 + 1 < N-1 and N//2 - 1 > 0
        
        ## Clamp one vertex at each end.
        constraints[0] = sweep_spine3d[0].copy()
        constraints[N-1] = sweep_spine3d[N-1].copy()
        ## Clamp the middle 3 vertices in a diagonal line.
        constraints[N//2] = sweep_spine3d[N//2].copy()
        constraints[ N//2 - 1 ] = sweep_spine3d[ N//2 ].copy()
        constraints[ N//2 - 1 ][1] = -dx
        constraints[ N//2 - 1 ][2] = -dx
        constraints[ N//2 + 1 ] = sweep_spine3d[ N//2 ].copy()
        constraints[ N//2 + 1 ][1] = dx
        constraints[ N//2 + 1 ][2] = dx
        
        constraints = list(constraints.items())
    
    else:
        raise RuntimeError("Unknown test type")
    
    
    print("sweep_spine3d:")
    pprint( sweep_spine3d )
    
    print("constraints:")
    pprint( constraints )
    
    opt_spine = optimize_sweep_spine3d( sweep_spine3d, constraints )
    
    print("opt_spine:")
    pprint( opt_spine )
    
    test_opt_vals( opt_spine )
    
    display_line( opt_spine )

def test_opt_vals( spine ):
    N = len(spine)
    
    test_type = 2
    if test_type == 1:
        constraints = [ (0,(0.,)), (N-1,(1.,)) ]
    elif test_type == 2:
        constraints = [ (0,(0.,)), (N-1,(0.,)), (N//2,(1.,)) ]
    else:
        raise RuntimeError("Unknown test type")
    
    print("test_opt_vals constraints:")
    print(constraints)
    results = diffuse_values_along_sweep_spine( spine, constraints )
    print("test_opt_vals results:")
    print(results)

def display_line( line, update_line_func = None ):
    from OpenGL.GLUT import glutInit, glutInitDisplayMode, glutMainLoop, glutPostRedisplay
    from OpenGL.GLUT import GLUT_RGBA, GLUT_DOUBLE, GLUT_DEPTH, GLUT_UP
    glutInit()
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH )
    from line_gui import LineWindow
    
    class LineWindowUpdatable( LineWindow ):
        def keyboardFunc( self, key, x, y ):
            if key in ('\n','\r',' ') and update_line_func is not None:
                self.line = list( update_line_func( self.line ) )
                glutPostRedisplay()
            else: LineWindow.keyboardFunc( self, key, x, y )
    
    w = LineWindowUpdatable()
    from mesh_gui import Camera
    w.camera.proj = Camera.Ortho( -.5, 1.5, -.5, 1.5 )
    w.line = list( line )
    
    glutMainLoop()

def main():
    #test2d()
    test3d()

if __name__ == "__main__": main()
