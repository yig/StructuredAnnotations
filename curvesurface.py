from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *

from pprint import pprint

try:
    from pydb import debugger
    
    ## Also add an exception hook.
    import pydb, sys
    sys.excepthook = pydb.exception_hook
    
except ImportError:
    import pdb
    def debugger():
        pdb.set_trace()


import scipy.sparse
assembly_matrix_t = scipy.sparse.lil_matrix
## The dok_matrix hung.  Maybe it was working, but at glacially speed.
#assembly_matrix_t = scipy.sparse.dok_matrix
solve_matrix_t = scipy.sparse.csc_matrix

kDefaultIterations = 4

def optimize_surface_linear_together( mesh, constraints, **kwargs ):
    '''
    Given a TriMesh 'mesh', Constraints 'constraints',
    minimizes linearized stretching and bending energy to find new mesh vertex positions.
    If 'graph_not_cotangent_weights' is False (the default is false), cotangent weights are used
    instead of graph laplacian weights.
    If 'vertex_normal_motion_only' is True (the default is false), vertices are only allowed to move
    in their normal directions.
    If 'stretching_weight' is None, there will be no stretching term.
    If 'bending_weight' is None, there will be no bending term.
    NOTE: At least one of 'stretching_weight' or 'bending_weight' must be non-None.
    
    Returns 0 (to be compatible with fibermesh.optimize_surface2()).
    
    NOTE: Modifies 'mesh' in place.
    '''
    
    kwargs.setdefault( "graph_not_cotangent_weights", True )
    kwargs.setdefault( "vertex_normal_motion_only", False )
    kwargs.setdefault( "stretching_weight", 1.0 )
    kwargs.setdefault( "bending_weight", None )
    kwargs.setdefault( "name", "" )
    
    graph_not_cotangent_weights = kwargs["graph_not_cotangent_weights"]
    vertex_normal_motion_only = kwargs["vertex_normal_motion_only"]
    name = kwargs["name"]
    
    from scipy import sparse
    from tictoc import tic, toc
    
    ## Mesh vertices should be unique (must they be, or only positions of constraints?)
    assert len(set([ tuple(v) for v in mesh.vs ])) == len( mesh.vs )
    ## We want our mesh vertices to be a 2d array.
    mesh.vs = asarrayf( mesh.vs )
    
    constraints.assert_self( mesh )
    
    #debugger()
    
    from fibermesh import gen_weighted_graph_laplacian_matrix, gen_cotangent_laplacian_matrix, Constraints
    
    if graph_not_cotangent_weights:
        gen_laplacian = gen_weighted_graph_laplacian_matrix
    else:
        ## Cotangent weights makes some faces turn extremely long and narrow and eventually collapse,
        ## and the system diverges.
        ## And that problem doesn't go away if we constrain some edges (average_edge_lengths).
        ## It does go away if we constrain vertices to only move in the normal direction, though
        ## (vertex_normal_motion_only).
        gen_laplacian = gen_cotangent_laplacian_matrix
    
    w_position = 1e7
    
    ## We can re-use the sparsity pattern of L
    N = len( mesh.vs )
    ## TODO Q: What am I doing wrong with the areas?
    ##         Either I should multiply the diagonal matrix areas^{-1} on the left of L2,
    ##         or else I should multiply it as a diagonal matrix between L*L.
    ##         And then what about L?
    tic( 'generate laplacian matrix' )
    L, vert_areas = gen_laplacian( mesh, None )
    toc()
    
    w_s = 0.
    if kwargs['stretching_weight'] is not None: w_s = kwargs['stretching_weight']
    assert w_s >= 0.
    
    w_b = 0.
    if kwargs['bending_weight'] is not None: w_b = kwargs['bending_weight']
    assert w_b >= 0.
    
    print('stretching weight:', w_s)
    print('bending weight:', w_b)
    assert w_s + w_b > 0.
    
    E = assembly_matrix_t( L.shape )
    
    if w_s > 0.:
        E = E + w_s*L
    
    if w_b > 0.:
        L2 = L*L
        E = E + w_b*L2
    
    import mymath
    xyz_hard_linear_constraints = list( constraints.mymath_xy_hard_and_linear_constraints() )
    xyz_hard_linear_constraints.append( constraints.mymath_z_hard_and_linear_constraints() )
    
    #print 'xyz_hard_linear_constraints:'
    #pprint( xyz_hard_linear_constraints )
    #debugger()
    
    #print 'general_linear_constraints:'
    #pprint( constraints.general_linear_constraints )
    
    ## We're going to assume that the system is coupled.
    ## So, put all degrees of freedom into a 3N system.
    
    tic( 'building A3' )
    A3 = assembly_matrix_t( ( 3*N, 3*N ), dtype = E.dtype )
    A3[:N,:N] = E
    A3[N:2*N,N:2*N] = E
    A3[2*N:3*N,2*N:3*N] = E
    toc()
    
    b3 = zeros( 3*N, A3.dtype )
    
    tic( 'accumulate constraints' )
    ## Accumulate all hard and linear constraints, adjusting indices.
    hard_constraints = []
    linear_constraints = []
    for coord in (0,1,2):
        for vert, rhs in xyz_hard_linear_constraints[ coord ][0]:
            hard_constraints.append( ( vert + coord*N, rhs ) )
        for verts, weights, rhs in xyz_hard_linear_constraints[ coord ][1]:
                linear_constraints.append( (
                    [ vi + coord*N for vi in verts ],
                    weights,
                    rhs
                    ) )
    
    ## Add the general linear constraints.
    #print 'before:', constraints.general_linear_constraints
    #print 'after:', map( Constraints.general_linear_constraint_collect_terms, constraints.general_linear_constraints )
    for verts_coords, weights, rhs in map( Constraints.general_linear_constraint_collect_terms, constraints.general_linear_constraints ):
        linear_constraints.append( (
            [ vi + coord*N for vi, coord in verts_coords ],
            weights,
            rhs
            ) )
    toc()
    
    tic( 'mymath.solve_system_with_constraints()' )
    solution3 = mymath.solve_system_with_constraints(
        A3, b3,
        hard_constraints,
        linear_constraints,
        'lsq',
        w_position
        )
    toc()
    
    vs = arrayf( mesh.vs )
    assert vs.shape == ( len( mesh.vs ), 3 )
    for i in range(N): vs[i,:] = solution3[ i ], solution3[ i + N ], solution3[ i + 2*N ]
    
    ## Verify that our hard constraints were met.
    for coord in (0,1,2):
        assert False not in [ ( vs[ vi ][ coord ] - pos )**2 < 1e-5 for vi, pos in xyz_hard_linear_constraints[coord][0] ]
        #for vi, pos in xyz_hard_linear_constraints[i][0]: vs[ vi, coord ] = pos
    
    if vertex_normal_motion_only:
        mesh.vs += ( ( vs - mesh.vs ) * mesh.vertex_normals ).sum(1).reshape((N,1)) * mesh.vertex_normals
    else:
        mesh.vs = vs
    
    mesh.positions_changed()
    
    mesh.write_OBJ( name + "optimized-linear.obj" )

def optimize_surface_linear_separate( mesh, constraints, **kwargs ):
    '''
    Given a TriMesh 'mesh', Constraints 'constraints',
    minimizes linearized stretching and bending energy to find new mesh vertex positions.
    If 'graph_not_cotangent_weights' is False (the default is false), cotangent weights are used
    instead of graph laplacian weights.
    If 'vertex_normal_motion_only' is True (the default is false), vertices are only allowed to move
    in their normal directions.
    If 'stretching_weight' is None, there will be no stretching term.
    If 'bending_weight' is None, there will be no bending term.
    NOTE: At least one of 'stretching_weight' or 'bending_weight' must be non-None.
    
    Returns 0 (to be compatible with fibermesh.optimize_surface2()).
    
    NOTE: Modifies 'mesh' in place.
    '''
    
    kwargs.setdefault( "graph_not_cotangent_weights", True )
    kwargs.setdefault( "vertex_normal_motion_only", False )
    kwargs.setdefault( "stretching_weight", 1.0 )
    kwargs.setdefault( "bending_weight", None )
    kwargs.setdefault( "name", "" )
    
    graph_not_cotangent_weights = kwargs["graph_not_cotangent_weights"]
    vertex_normal_motion_only = kwargs["vertex_normal_motion_only"]
    name = kwargs["name"]
    
    from scipy import sparse
    from tictoc import tic, toc
    
    ## Mesh vertices should be unique (must they be, or only positions of constraints?)
    assert len(set([ tuple(v) for v in mesh.vs ])) == len( mesh.vs )
    ## We want our mesh vertices to be a 2d array.
    mesh.vs = asarrayf( mesh.vs )
    
    constraints.assert_self( mesh )
    
    #debugger()
    
    from fibermesh import gen_weighted_graph_laplacian_matrix, gen_cotangent_laplacian_matrix, Constraints
    
    if graph_not_cotangent_weights:
        gen_laplacian = gen_weighted_graph_laplacian_matrix
    else:
        ## Cotangent weights makes some faces turn extremely long and narrow and eventually collapse,
        ## and the system diverges.
        ## And that problem doesn't go away if we constrain some edges (average_edge_lengths).
        ## It does go away if we constrain vertices to only move in the normal direction, though
        ## (vertex_normal_motion_only).
        gen_laplacian = gen_cotangent_laplacian_matrix
    
    w_position = 1e7
    
    ## We can re-use the sparsity pattern of L
    N = len( mesh.vs )
    ## TODO Q: What am I doing wrong with the areas?
    ##         Either I should multiply the diagonal matrix areas^{-1} on the left of L2,
    ##         or else I should multiply it as a diagonal matrix between L*L.
    ##         And then what about L?
    tic( 'generate laplacian matrix' )
    L, vert_areas = gen_laplacian( mesh, None )
    toc()
    
    w_s = 0.
    if kwargs['stretching_weight'] is not None: w_s = kwargs['stretching_weight']
    assert w_s >= 0.
    
    w_b = 0.
    if kwargs['bending_weight'] is not None: w_b = kwargs['bending_weight']
    assert w_b >= 0.
    
    print('stretching weight:', w_s)
    print('bending weight:', w_b)
    assert w_s + w_b > 0.
    
    E = assembly_matrix_t( L.shape )
    
    if w_s > 0.:
        E = E + w_s*L
    
    if w_b > 0.:
        L2 = L*L
        E = E + w_b*L2
    
    import mymath
    xyz_hard_linear_constraints = list( constraints.mymath_xy_hard_and_linear_constraints() )
    xyz_hard_linear_constraints.append( constraints.mymath_z_hard_and_linear_constraints() )
    
    #print 'xyz_hard_linear_constraints:'
    #pprint( xyz_hard_linear_constraints )
    #debugger()
    
    #print 'general_linear_constraints:'
    #pprint( constraints.general_linear_constraints )
    
    ## We're going to assume that the system is NOT coupled.
    assert len( constraints.general_linear_constraints ) == 0
    ## So, put all degrees of freedom into an N-by-N system.
    
    b = zeros( N, E.dtype )
    solution3 = zeros( (N,3), E.dtype )
    
    for coord in (0,1,2):
        tic( 'accumulate constraints' )
        ## Accumulate all hard and linear constraints, adjusting indices.
        hard_constraints = []
        linear_constraints = []
        for vert, rhs in xyz_hard_linear_constraints[ coord ][0]:
            hard_constraints.append( ( vert, rhs ) )
        for verts, weights, rhs in xyz_hard_linear_constraints[ coord ][1]:
                linear_constraints.append( (
                    [ vi for vi in verts ],
                    weights,
                    rhs
                    ) )
        toc()
        
        tic( 'mymath.solve_system_with_constraints()' )
        solution3[:,coord] = mymath.solve_system_with_constraints(
            E, b,
            hard_constraints,
            linear_constraints,
            'lsq',
            w_position
            )
        toc()
    
    vs = arrayf( mesh.vs )
    assert vs.shape == ( len( mesh.vs ), 3 )
    vs[:] = solution3[:]
    
    ## Verify that our hard constraints were met.
    for coord in (0,1,2):
        assert False not in [ ( vs[ vi ][ coord ] - pos )**2 < 1e-5 for vi, pos in xyz_hard_linear_constraints[coord][0] ]
        #for vi, pos in xyz_hard_linear_constraints[i][0]: vs[ vi, coord ] = pos
    
    if vertex_normal_motion_only:
        mesh.vs += ( ( vs - mesh.vs ) * mesh.vertex_normals ).sum(1).reshape((N,1)) * mesh.vertex_normals
    else:
        mesh.vs = vs
    
    mesh.positions_changed()
    
    mesh.write_OBJ( name + "optimized-linear.obj" )

def optimize_surface_simple( mesh, xyz_constraints, **kwargs ):
    '''
    Given a TriMesh 'mesh', a list of ( index, xyz ) tuples 'xyz_constraints',
    minimizes linearized stretching and bending energy to find new mesh vertex positions.
    If 'graph_not_cotangent_weights' is False (the default is false), cotangent weights are used
    instead of graph laplacian weights.
    If 'vertex_normal_motion_only' is True (the default is false), vertices are only allowed to move
    in their normal directions.
    If 'stretching_weight' is None, there will be no stretching term.
    If 'bending_weight' is None, there will be no bending term.
    NOTE: At least one of 'stretching_weight' or 'bending_weight' must be non-None.
    
    The process will run in a loop, iterating 'max_iterations' times (default 1).
    
    Returns the maximum difference in vertex positions
    due to the last iteration performed.
    
    NOTE: Modifies 'mesh' in place.
    '''
    
    from tictoc import tic, toc
    tic( 'optimize_surface_simple()' )
    
    kwargs.setdefault( "graph_not_cotangent_weights", True )
    kwargs.setdefault( "vertex_normal_motion_only", False )
    kwargs.setdefault( "stretching_weight", 1.0 )
    kwargs.setdefault( "bending_weight", None )
    kwargs.setdefault( "max_iterations", 1 )
    kwargs.setdefault( "absolute_convergence_threshold", 0. )
    kwargs.setdefault( "name", "" )
    
    print('kwargs:')
    pprint( kwargs )
    
    max_iterations = kwargs["max_iterations"]
    absolute_convergence_threshold = kwargs["absolute_convergence_threshold"]
    graph_not_cotangent_weights = kwargs["graph_not_cotangent_weights"]
    vertex_normal_motion_only = kwargs["vertex_normal_motion_only"]
    name = kwargs["name"]
    
    from scipy import sparse
    
    ## Mesh vertices should be unique (must they be, or only positions of constraints?)
    assert len(set([ tuple(v) for v in mesh.vs ])) == len( mesh.vs )
    ## We want our mesh vertices to be a 2d array.
    mesh.vs = asarrayf( mesh.vs )
    
    #debugger()
    
    from fibermesh import gen_weighted_graph_laplacian_matrix, gen_cotangent_laplacian_matrix
    
    if graph_not_cotangent_weights:
        gen_laplacian = gen_weighted_graph_laplacian_matrix
    else:
        ## Cotangent weights makes some faces turn extremely long and narrow and eventually collapse,
        ## and the system diverges.
        ## And that problem doesn't go away if we constrain some edges (average_edge_lengths).
        ## It does go away if we constrain vertices to only move in the normal direction, though
        ## (vertex_normal_motion_only).
        gen_laplacian = gen_cotangent_laplacian_matrix
    
    positions4laplacian = optimize_positions_for_mass_laplacian_simple_constraints
    
    w_position = 1e7
    
    w_s = 0.
    if kwargs['stretching_weight'] is not None: w_s = kwargs['stretching_weight']
    assert w_s >= 0.
    
    w_b = 0.
    if kwargs['bending_weight'] is not None: w_b = kwargs['bending_weight']
    assert w_b >= 0.
    
    print('stretching weight:', w_s)
    print('bending weight:', w_b)
    assert w_s + w_b > 0.
    
    
    N = len( mesh.vs )
    #L = assembly_matrix_t( (N,N) )
    L = None
    
    last_vs = arrayf( mesh.vs )
    max_diff = absolute_convergence_threshold + 1.
    for iter in range(max_iterations):
        if max_diff < absolute_convergence_threshold: break
        tic( 'optimize_surface_simple() iteration' )
        
        ## TODO Q: What am I doing wrong with the areas?
        ##         Either I should multiply the diagonal matrix areas^{-1} on the left of L2,
        ##         or else I should multiply it as a diagonal matrix between L*L.
        ##         And then what about L?
        tic( 'generate laplacian matrix' )
        L, vert_areas = gen_laplacian( mesh, L )
        toc()
        
        vs = positions4laplacian( mesh, L, vert_areas, xyz_constraints, w_position, w_s, w_b )
        if vertex_normal_motion_only:
            mesh.vs += ( ( vs - mesh.vs ) * mesh.vertex_normals ).sum(1).reshape((N,1)) * mesh.vertex_normals
        else:
            mesh.vs = vs
        
        mesh.positions_changed()
        
        mesh.write_OBJ( name + 'optimized%s.obj' % (iter,) )
        max_diff = sqrt( max( ( ( last_vs - mesh.vs )**2 ).sum(1) ) )
        print('max( | mesh.vs - last_vs | ):', max_diff)
        last_vs = arrayf( mesh.vs )
        
        toc()
    
    toc()
    return max_diff
    
def optimize_positions_for_mass_laplacian_simple_constraints( mesh, L, vert_areas, xyz_constraints, w_position, w_s, w_b ):
    #debugger()
    
    from tictoc import tic, toc
    tic( 'optimize_positions_for_mass_laplacian_simple_constraints' )
    
    from scipy import sparse
    import scipy.sparse.linalg.dsolve as linsolve
    
    ## Make sure our dimensions are correct.
    N = len( mesh.vs )
    assert L.shape[0] == N
    assert L.shape[1] == N
    assert len( vert_areas ) == N
    
    tic( 'build E from L' )
    E = assembly_matrix_t( L.shape )
    
    if w_s > 0.:
        E = E + w_s*L
    
    if w_b > 0.:
        L2 = L*L
        E = E + w_b*L2
    
    toc()
    
    
    vs = arrayf( mesh.vs ) #zerosf( rhs.shape )
    assert vs.shape == ( N, 3 )
    
    ## Fast constraint setting code comes from "recovery.py".
    tic( 'Add constraints' )
    A = E
    ## Add hard xyz constraints.
    
    ## Zero the constrained rows and set the diagonal to 1.
    xyz_constraint_indices = asarray([ index for index, xyz in xyz_constraints ])
    xyz_constraint_values = asarray([ xyz for index, xyz in xyz_constraints ])
    
    ## Zero the rows.
    diag = ones( N )
    diag[ xyz_constraint_indices ] = 0.
    #D = sparse.identity( system.shape[0] ).tocsr()
    #D.setdiag( diag )
    D = sparse.coo_matrix( ( diag, ( arange( N ), arange( N ) ) ) ).tocsr()
    A = D * A
    
    ## Set the constraints' diagonals to 1.
    #diag = zeros( rhs.shape[0] )
    #diag[ value_constraint_indices ] = 1.
    #D.setdiag( diag )
    #D.setdiag( 1. - diag )
    D = sparse.coo_matrix( ( ones( xyz_constraint_indices.shape[0] ), ( xyz_constraint_indices, xyz_constraint_indices ) ), shape = A.shape )
    A = A + D
    
    
    b = zeros( (N,3), A.dtype )
    
    ## Set corresponding entries of the right-hand-side vector to the constraint value.
    #for index, value in value_constraints_dict.iteritems(): rhs[ index ] = value
    b[ xyz_constraint_indices ] = xyz_constraint_values
    
    toc()
    
    
    tic( 'Factor and solve' )
    assert A.shape == (N,N)
    A = solve_matrix_t( A ) #, dtype = 'F' )
    # linsolve.use_solver( use = {'useUmfpack' : True} )
    factorized = linsolve.factorized( A )
    
    ## Solve for new xy positions.
    for i in range(3): vs[:,i] = factorized( b[:,i].flatten() )
    
    toc()
    
    
    ## Verify that our hard constraints were met.
    for coord in (0,1,2):
        assert False not in [ ( vs[ vi ][ coord ] - pos[ coord ] )**2 < 1e-5 for vi, pos in xyz_constraints ]
        #for vi, pos in xyz_constraints: vs[ vi, coord ] = pos
    
    toc()
    return vs

def optimize_surface_linear( mesh, constraints, **kwargs ):
    if( len( constraints.general_linear_constraints ) > 0 ):
        optimize_surface_linear_together( mesh, constraints, **kwargs )
    else:
        if constraints.has_non_xyz_position_constraints():
            optimize_surface_linear_separate( mesh, constraints, **kwargs )
        else:
            optimize_surface_simple( mesh, constraints.xyz_position_constraints(), **kwargs )
