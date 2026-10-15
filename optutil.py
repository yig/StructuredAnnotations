from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.

from helpers import x_intercept

def find_zero( func, initial_value, initial_step, **kwargs ):
    '''
    Given an objective function 'func' and an initial scalar value 'initial_value',
    solves for 'func' = 0 by searching for a better scalar value.
    The first step taken is 'initial_step' (func( initial_value + initial_step ) will be evaluated).
    Returns a tuple containing ( the optimized scalar value, and whether or not the search
    converged ).
    Note that if the search diverges, the last step taken will be returned.
    
    Optional parameters:
        
    
    tested
    '''
    
    kwargs.setdefault( 'verbose', False )
    kwargs.setdefault( 'xthresh', 1e-5 )
    kwargs.setdefault( 'fthresh', 1e-8 )
    kwargs.setdefault( 'max_iter', 20 )
    
    verbose = kwargs['verbose']
    
    XSTEP_THRESHOLD = kwargs['xthresh']
    F_THRESHOLD = kwargs['fthresh']
    MAX_ITERATIONS = kwargs['max_iter']
    
    x_step = initial_step
    x_n = initial_value
    f_n = func( x_n )
    if verbose: print('func( %s ): %s' % (x_n, f_n))
    x_n_p_1 = x_n + x_step
    
    iteration = 0
    diverged = False
    while True:
        if verbose: print('-- iteration %s, x: %s --' % (iteration, x_n_p_1))
        iteration += 1
        
        xdiff = abs( x_n_p_1 - x_n )
        if verbose: print('| x_n+1 - x_n |:', xdiff)
        
        f_n_p_1 = func( x_n_p_1 )
        if verbose: print('func( %s ): %s' % (x_n_p_1, f_n_p_1))
        fdiff = abs( f_n_p_1 - f_n )
        if verbose: print('| f_n+1 - f_n |:', fdiff)
        
        if xdiff < XSTEP_THRESHOLD: break
        if fdiff < F_THRESHOLD: break
        
        if iteration >= MAX_ITERATIONS:
            print('Diverged.')
            diverged = True
            break
        
        ## Linear estimate for correct 'x' based on previous iteration
        ## The two points of the line are
        ## (x_n, f_n) and (x_n_p_1, f_n_p_1).
        x_n_p_2 = x_intercept( (x_n, f_n), (x_n_p_1, f_n_p_1) )
        x_n, f_n = x_n_p_1, f_n_p_1
        x_n_p_1 = x_n_p_2
        del x_n_p_2
    
    return x_n_p_1, not diverged

def binary_search_zero( func, lower_value, upper_value, **kwargs ):
    '''
    Given an objective scalar function 'func' and lower, upper values [lower_value, upper_value],
    finds 'func' = 0 via binary search between lower/upper_value.
    Returns a tuple containing ( the optimized scalar value, and whether or not the search
    converged ).
    Note that if the search diverges, the last step taken will be returned.
    
    Optional parameters:
        
    
    untested
    '''
    
    assert lower_value < upper_value
    
    kwargs.setdefault( 'verbose', False )
    kwargs.setdefault( 'xthresh', 1e-5 )
    kwargs.setdefault( 'fthresh', 1e-8 )
    kwargs.setdefault( 'max_iter', 20 )
    
    verbose = kwargs['verbose']
    
    XSTEP_THRESHOLD = kwargs['xthresh']
    F_THRESHOLD = kwargs['fthresh']
    MAX_ITERATIONS = kwargs['max_iter']
    
    lower_func = func( lower_value )
    upper_func = func( upper_value )
    if verbose: print('func( %s ): %s' % (lower_value, lower_func))
    if verbose: print('func( %s ): %s' % (upper_value, upper_func))
    
    middle_value = lower_value
    middle_func = lower_func
    
    iteration = 0
    diverged = False
    while True:
        if verbose: print('-- iteration %s, [lower,upper] values: [%s,%s] --' % (iteration, lower_value, upper_value))
        iteration += 1
        
        xdiff = upper_value - lower_value
        if verbose: print('| upper value - lower value |:', xdiff)
        
        fdiff = abs( upper_func - lower_func )
        if verbose: print('| f( upper value ) - f( lower value ) |:', fdiff)
        
        if xdiff < XSTEP_THRESHOLD: break
        if fdiff < F_THRESHOLD: break
        
        if iteration >= MAX_ITERATIONS:
            print('Diverged.')
            diverged = True
            break
        
        middle_value = .5 * ( lower_value + upper_value )
        middle_func = func( middle_value )
        
        if middle_func < 0.:
            #old_value, old_func = lower_value, lower_func
            lower_value, lower_func = middle_value, middle_func
        else:
            #old_value, old_func = upper_value, upper_func
            upper_value, upper_func = middle_value, middle_func
    
    return middle_value, not diverged

def find_minimum( func, initial_value, initial_step, **kwargs ):
    '''
    Given an objective function 'func' and an initial scalar value 'initial_value',
    finds the minimum value of 'func' by searching for a better scalar value.
    The first step taken is 'initial_step' (func( initial_value + initial_step ) will be evaluated).
    Returns a tuple containing ( the optimized scalar value, and whether or not the search
    converged ).
    Note that if the search diverges, the last step taken will be returned.
    
    Optional parameters:
        
    
    tested
    '''
    
    kwargs.setdefault( 'verbose', False )
    kwargs.setdefault( 'xthresh', 1e-5 )
    kwargs.setdefault( 'fthresh', 1e-8 )
    kwargs.setdefault( 'max_iter', 20 )
    
    verbose = kwargs['verbose']
    
    XSTEP_THRESHOLD = kwargs['xthresh']
    F_THRESHOLD = kwargs['fthresh']
    MAX_ITERATIONS = kwargs['max_iter']
    
    x_ns = []
    f_ns = []
    
    x_ns.append( initial_value )
    f_ns.append( func( x_ns[-1] ) )
    if verbose: print('func( %s ): %s' % (x_ns[-1], f_ns[-1]))
    
    x_ns.append( x_ns[-1] + initial_step )
    f_ns.append( func( x_ns[-1] ) )
    if verbose: print('func( %s ): %s' % (x_ns[-1], f_ns[-1]))
    
    iteration = 0
    diverged = False
    while True:
        assert len( x_ns ) == len( f_ns )
        
        if verbose: print('-- iteration %s --' % (iteration,))
        iteration += 1
        
        xdiff = abs( x_ns[-1] - x_ns[-2] )
        if verbose: print('| x_n+1 - x_n |:', xdiff)
        fdiff = abs( f_ns[-1] - f_ns[-2] )
        if verbose: print('| f_n+1 - f_n |:', fdiff)
        
        if xdiff < XSTEP_THRESHOLD: break
        if fdiff < F_THRESHOLD: break
        
        if iteration >= MAX_ITERATIONS:
            print('Diverged.')
            diverged = True
            break
        
        ## Otherwise, the gradient is not zero.
        while len( f_ns ) < 3:
            assert iteration == 1
            ## Evaluate again, because we need three values
            x_ns.append( x_ns[-1] + initial_step )
            f_ns.append( func( x_ns[-1] ) )
            if verbose: print('func( %s ): %s' % (x_ns[-1], f_ns[-1]))
        
        assert len( x_ns ) == len( f_ns )
        ## Linear estimate for 'x' based on backwards differencing of previous iterations.
        ## The two points of the line are
        ## (x_ns[-2], f_ns[-2] - f_ns[-3]) and (x_ns[-1], f_ns[-1] - f_ns[-2]).
        x_ns.append( x_intercept( (x_ns[-2], f_ns[-2] - f_ns[-3]), (x_ns[-1], f_ns[-1] - f_ns[-2]) ) )
        f_ns.append( func( x_ns[-1] ) )
        if verbose: print('func( %s ): %s' % (x_ns[-1], f_ns[-1]))
        del x_ns[:-3]
        del f_ns[:-3]
    
    return x_ns[-1], not diverged

def minimize( dofs, forces, constraints ):
    '''
    Mini meshopt framework.
    '''
    
    x0 = dofs.ravel.copy()
    N = x0.shape[0]
    
    ## TODO: test for convergence
    while True:
        J = spzero( (N,N) )
        f = zerosf( N )
        
        for force in forces:
            for element in forces.elements( dofs ):
                ## TODO: gather dofs, compute force, hessian, scatter into the sparse matrices
                forces.compute_gradient_and_hessian()
                pass
        
        for dof_index in constraints:
            J[ dof_index, : ] = 0.
            ## TODO: Can I zero the column like this without shifting
            ##       values in the column to the right-hand side?
            J[ :, dof_index ] = 0.
            J[ dof_index, dof_index ] = 1.
            f[ dof_index ] = 0.
        
        ## TODO: take a Newton step
        ## TODO: should I use tnc or TAO?
