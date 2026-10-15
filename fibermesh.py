from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *
from helpers import *
from geometry import line_segments_cross

from pprint import pprint

#from pydb import debugger

import sys, pdb
def do_debugger(type, value, tb):
    pdb.pm()
#sys.excepthook = do_debugger

import scipy.sparse
assembly_matrix_t = scipy.sparse.lil_matrix
## The dok_matrix hung.  Maybe it was working, but at glacially speed.
#assembly_matrix_t = scipy.sparse.dok_matrix
solve_matrix_t = scipy.sparse.csc_matrix

kDefaultIterations = 4
## This applies to gen_cotangent_laplacian_matrix() and gen_weighted_graph_laplacian_matrix().
## Setting it to True seems to generate nicer results (and still converges...).
## (This shouldn't make a difference beyond numerical stability, although there is some interaction
## with constraint weights.)
kScaleLaplacianWeights = False
## This doesn't seem to have much of an effect.
kClipNegativeCurvatureConstraints = True

def inflate( curve, iterations = kDefaultIterations, name = '' ):
    '''
    Given a curve (simple, closed list of (x,y) points), inflates the curve according to the
    FiberMesh algorithm.
    '''
    
    if len( name ) > 0: name = name + '_'
    
    ## Steps to FiberMesh initial mesh creation.
    ## 1 Generate a regular mesh.
    ## 2 For each point on the curve, find the closest point on the regular mesh.
    ## 3 Remove all points from the regular mesh outside of the curve (and not found in step 2).
    ## 4 Glue two copies of the regular grid to each other along the boundary (possibly offset the back copy away from the front in the look direction).
    ## 5 Generate the uniform Laplacian matrix for the mesh.
    ## 6 Find the current Laplacian magnitudes (LMs) and edge lengths along the boundary curve.
    ## 7 Using the matrix from Step 5, solve for LMs and edge lengths everywhere.  (In the first iteration, only set target curvature and edge lengths for vertices along the curve.)
    ## 8 Using the target edge lengths (from step 7) everywhere, find target edge vectors.
    ## 9 Using the target edge vectors (from step 8) near the boundary curve, and the target laplacians computed from laplacian magnitudes (from step 7) and normals, find new vertex positions.
    ## 10 Go to step 6 (until convergence).
    
    ## Alternate recipe, where the Laplacian matrix for the mesh is geometry dependent.
    ## 1
    ## 2
    ## 3
    ## 4
    ## 5b Generate the geometry-dependent Laplacian matrix for the mesh.
    ## 6b Find the current Laplacian magnitudes (LMs) along the boundary curve.
    ## 7b Using the matrix from Step 5, solve for LMs everywhere.  (In the first iteration, only set target curvature for vertices along the curve.)
    ## 8b no-op
    ## 9b Using the target laplacians computed from laplacian magnitudes (from step 7) and normals, find new vertex positions.
    ## 10
    
    #debugger()
    
    ## 1-4
    mesh, boundary_vertices, front_vertices, back_vertices = gen_glued_mesh( curve )
    mesh.write_OBJ( name + 'glued.obj' )
    #return mesh
    
    boundary_vertices, boundary_vi_to_curve_vi = attach_curve2d_to_mesh_vertices( curve, mesh, boundary_vertices )
    
    boundary_vi_to_pos = [ append( curve[ cvi ], [0.] ) for cvi in boundary_vi_to_curve_vi ]
    boundary_vi_to_curvature = [
        curvature_1d(
            curve[ (cvi + len(curve) - 1) % len(curve) ],
            curve[ cvi ],
            curve[ (cvi + 1) % len(curve) ]
            )
        for cvi in boundary_vi_to_curve_vi
        ]
    
    print('curve curvatures:', boundary_vi_to_curvature)
    
    #debugger()
    #mesh = optimize_surface( mesh, boundary_vertices, boundary_vi_to_pos, boundary_vi_to_curvature, iterations, name )
    constraints = Constraints()
    constraints.xy_constraint_vertices = boundary_vertices
    constraints.xy_constraint_positions = array( boundary_vi_to_pos )[:,:2]
    constraints.z_constraint_vertices = boundary_vertices
    constraints.z_constraint_values = array( boundary_vi_to_pos )[:,2]
    mesh = optimize_surface2( mesh, constraints, max_iterations = iterations, name = name )
    return mesh

    
def attach_curve2d_to_mesh_vertices( curve, mesh, candidate_vertices ):
    '''
    For each mesh vertex index in the list 'candidate_vertices'
    finds the index of the closest 2d point in curve.
    Returns two lists: the first is a subset of 'candidate_vertices' excluding vertex indices which
    map to the same curve point.  The second is the list of indices into 'curve' which are closest
    to the corresponding vertex indices in the modified 'candidate_vertices'.
    '''
    
    curve = asarray( curve )
    
    assert len(set( candidate_vertices )) == len( candidate_vertices )
    
    ## For each mesh vertex index in 'boundary_vertices', the corresponding vertex index of the
    ## closest point on 'curve'.
    candidate_vi_to_curve_vi = [ closest_point_index( mesh.vs[ mvi ][:2], curve ) for mvi in candidate_vertices ]
    ## We don't want two points on the mesh boundary to map to the same curve point.
    ## TODO: Something proper, like parameterizing the curve and boundary vertices
    ##       or deforming the curve into the boundary vertices.
    # while len(set( candidate_vi_to_curve_vi )) != len( candidate_vi_to_curve_vi ):
    from helpers import histogram
    hist = histogram( candidate_vi_to_curve_vi )
    
    ## One approach: delete curve points until mesh boundary vertices attach to
    ##               unique curve points.  This may remove all points of the curve!
    '''
    ## Iterate over indices to delete in reverse order, because deleting from back to front
    ## is stable.
    to_delete = [ index for (index,count) in hist.iteritems() if count > 1 ]
    to_delete.sort()
    print 'to_delete:', to_delete
    to_delete.reverse()
    curve = list( curve )
    for index in to_delete: del curve[ index ]
    assert len( curve ) > 0
    print 'len( curve ):', len( curve )
    '''
    
    ## Another approach: remove candidate vertices until mesh candidate vertices attach to
    ##                   unique curve points.
    ## Remove the M-1 farther candidate vertices for each M-covered curve point.
    deleted = 0
    original_count = len( candidate_vi_to_curve_vi )
    for curvi,count in hist.items():
        assert count >= 1
        if count == 1: continue
        
        deleted += count-1
        ## Collect the list of mesh vertex indices in candidate_vertices that map to the same
        ## curve point.
        ## Pair the mesh vertex indices with the distances to the point.
        dist_candis = [
            ( mag2( mesh.vs[mvi][:2] - curve[ curvi ] ), candi )
            for candi, mvi in enumerate( candidate_vertices )
            if candidate_vi_to_curve_vi[candi] == curvi
            ]
        dist_candis.sort()
        ## Keep the closest index
        del dist_candis[0]
        ## Iterate over indices to delete in reverse order, because deleting from back to front
        ## is stable.
        to_delete = [ candi for dist2, candi in dist_candis ]
        to_delete.sort()
        to_delete.reverse()
        for candi in to_delete:
            del candidate_vertices[ candi ]
            del candidate_vi_to_curve_vi[ candi ]
    print('Deleted %s/%s vertices from the set of attach-curve-candidate vertices.' % (deleted, original_count))
    
    ## Repeat.
    ## UPDATE: actually, I shouldn't have to.
    # candidate_vi_to_curve_vi = [ closest_point_index( mesh.vs[ mvi ][:2], curve ) for mvi in candidate_vertices ]
    assert candidate_vi_to_curve_vi == [ closest_point_index( mesh.vs[ mvi ][:2], curve ) for mvi in candidate_vertices ]
    # end while
    
    assert len(set( candidate_vi_to_curve_vi )) == len( candidate_vi_to_curve_vi )
    return candidate_vertices, candidate_vi_to_curve_vi

def gen_glued_mesh( curve, N = None ):
    '''
    Given a closed loop as a list of xy points 'curve',
    creates a mesh from two discs glued to each other along a boundary roughly in the shape of
    'curve'.
    Returns ( mesh, boundary_vertices, front_vertices, back_vertices ),
    where 'mesh' is the created TriMesh, 'boundary_vertics' is the list of vertex indices (in mesh.vs)
    which followed the glued boundary, 'front_vertices' is the list of vertex indices (in mesh.vs)
    on the disc offset in the positive z direction, and 'back_vertices' are the vertex indices
    of the disc offset in the negative z direction.  Note that the each vertex index in
    'front_vertics' is the duplicate vertex of the vertex index at the same offset in
    'back_vertices'.
    
    Optional parameter N specifies the number of vertices in the mesh.
    '''
    
    ## 1
    ## I increased 'eps' to 1 (from .01) because the gen_equilateral_grid() function creates
    ## half-triangles instead of a staggered grid, and I only want real equilateral triangles.
    eps = .01
    xmin = min( [ x for x,y in curve ] ) - eps
    xmax = max( [ x for x,y in curve ] ) + eps
    ymin = min( [ y for x,y in curve ] ) - eps
    ymax = max( [ y for x,y in curve ] ) + eps
    
    '''
    ## Keep things square.
    if xmax-xmin > ymax-ymin:
        ymax = ymin + xmax-xmin
    else:
        xmax = xmin + ymax-ymin
    '''
    
    ## We're scaling the grid to fit around the curve; this scaling is nonuniform, so the
    ## resulting triangles will not have 2 equals sides (the gen_regular_grid() function generates
    ## right angle, not equilateral, triangles).
    ## TODO: try meshopt/bin/generate_equilateral_triangular.py
    ## A 5x5 grid does work, and the system is barely small enough to keep a human from going
    ## insane when comparing all the values between, say, the C++ and python version of inflation.
    #mesh = gen_regular_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 5, yverts = 5 )
    #mesh = gen_regular_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 10, yverts = 10 )
    #mesh = gen_regular_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 15, yverts = 15 )
    #mesh = gen_regular_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 20, yverts = 20 )
    #mesh = gen_regular_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 3, yverts = 3 )
    #mesh = gen_regular_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 3, yverts = 5 )
    #mesh = gen_regular_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 5, yverts = 5 )
    #mesh = gen_equilateral_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 5, yverts = 5 )
    #mesh = gen_equilateral_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 10, yverts = 10 )
    #mesh = gen_equilateral_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 20, yverts = 20 )
    #mesh = gen_equilateral_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = 3, yverts = 3 )
    
    if N is None:
        ## N = 100 is barely small enough to make the python optimize_surface2() tolerable.
        #N = 100
        N = 200
        #N = 400
        #N = 2000
    NX = max( 3, sqrt( N / (ymax-ymin) * (xmax-xmin) ) )
    NY = max( 3, NX * (ymax-ymin)/(xmax-xmin) )
    ## Round to the nearest integer.
    NX = int( round( NX ) )
    NY = int( round( NY ) )
    mesh = gen_regular_grid( xmin = xmin, xmax = xmax, ymin =  ymin, ymax = ymax, xverts = NX, yverts = NY )
    
    #debugger()
    
    ## 2
    mesh_verts_2d = asarray( mesh.vs )[:,:2]
    ## For each point on the curve, the index of the closest vertex in the mesh.
    closest_pts = [ closest_point_index( pt, mesh_verts_2d ) for pt in curve ]
    
    ## 3
    ## For each edge of the mesh intersected by the curve, include the edge's vertices.
    intersected_vertices = set()
    for edge in mesh.edges:
        for c in range(len(curve)):
            if line_segments_cross(
                ( mesh.vs[ edge[0] ], mesh.vs[ edge[1] ] ),
                ( curve[c], curve[ (c+1) % len( curve ) ] )
                ):
                intersected_vertices.add( edge[0] )
                intersected_vertices.add( edge[1] )
                break
    ## For each vertex in the mesh, a boolean stating whether each vertex is inside the curve.
    is_inside = [ contour_contains( curve, vert ) for vert in mesh_verts_2d ]
    ## For each vertex in the mesh, a boolean stating whether to keep the vertex.
    ## (We keep vertices if they are inside the curve or were chosen as the closest point for a
    ## point on the curve.)
    keep_map = [ is_inside[i] or i in closest_pts or i in intersected_vertices for i in range( len( mesh.vs ) ) ]
    ## Remove points from the mesh.
    mesh.remove_vertex_indices( [ i for i in range( len( mesh.vs ) ) if not keep_map[i] ] )
    ## Also remove dangling ears, vertices with only two neighbors (which must be on the boundary).
    assert False not in [ mesh.vertex_valence( i ) >= 2 for i in range(len( mesh.vs ) ) ]
    mesh.remove_vertex_indices( [ i for i in mesh.boundary_vertices() if mesh.vertex_valence( i ) == 2 ] )
    assert False not in [ mesh.vertex_valence( i ) >= 3 for i in range(len( mesh.vs ) ) ]
    ## These variables are no longer valid.
    del closest_pts
    del mesh_verts_2d
    
    ## Return a 2D mesh (before duplicating the non-boundary vertices).
    #return mesh, None, None, None
    #debugger()
    
    ## 4
    ## Duplicate non-boundary vertices.
    boundary_vertices = mesh.boundary_vertices()
    duplicate_map = [True] * len(mesh.vs)
    for vi in boundary_vertices: duplicate_map[vi] = False
    ## For each vertex in the mesh, store the new, duplicate vertex index (and -1 if not duplicated).
    old2dup = [-1] * len( mesh.vs )
    mesh.vs = list( mesh.vs )
    for i in range(len(mesh.vs)):
        if duplicate_map[i]:
            old2dup[i] = len( mesh.vs )
            mesh.vs.append( array( mesh.vs[ i ] ) )
    assert False not in [ -1 != old2dup[i] or not duplicate_map[i] for i in range(len(old2dup)) ]
    ## Every remaining triangle must include a vertex that was duplicated.
    assert False not in [ True in [ duplicate_map[i] for i in tri ] for tri in mesh.faces ]
    ## Duplicate every triangle and use as many duplicate indices as possible (boundary vertices
    ## don't get duplicated).
    old_nfaces = len(mesh.faces)
    for tri_index in range(len(mesh.faces)):
        tri = mesh.faces[ tri_index ]
        assert True in [ duplicate_map[i] for i in tri ]
        newtri = list(tri)
        for i in range(3):
            if duplicate_map[ tri[i] ]: newtri[i] = old2dup[ tri[i] ]
        assert -1 not in newtri
        ## Also, flip the orientation of the glued faces.
        mesh.faces.append( ( newtri[0], newtri[2], newtri[1] ) )
        #mesh.faces.append( ( newtri[0], newtri[1], newtri[2] ) )
    
    #mesh.faces = mesh.faces[:old_nfaces]
    
    #debugger()
    
    ## Offset old faces in the positive z direction, and new faces in the negative z direction.
    front_vertices = []
    back_vertices = []
    ## Our offset will be 10% of the average side length (if it were a rectangle).
    z_off = .1 * ( fabs( xmax - xmin ) + fabs( ymax - ymin ) )/2
    for old, dup in enumerate( old2dup ):
        ## Skip boundary vertices
        if -1 == dup: continue
        mesh.vs[old][2] = z_off
        mesh.vs[dup][2] = -z_off
        
        front_vertices.append( old )
        back_vertices.append( dup )
    
    mesh.topology_changed()
    return mesh, boundary_vertices, front_vertices, back_vertices

def optimize_surface( mesh, boundary_vertices, boundary_vertices_positions, boundary_vertices_curvatures, iterations = kDefaultIterations, name = '' ):
    '''
    Given a TriMesh 'mesh', boundary silhouette vertex indices 'boundary_vertices',
    and corresponding positions 'boundary_vertices_positions' and
    scalar curvature values 'boundary_vertices_curvatures',
    perform steps 5-10 of inflate().
    
    NOTE: Modifies 'mesh' in place and also returns it.
    '''
    
    from scipy import sparse
    
    ## Boundary vertices must be unique
    assert len(set(boundary_vertices)) == len(boundary_vertices)
    ## Boundary vertex positions must be unique
    assert len(set([ tuple(v) for v in boundary_vertices_positions ])) == len( boundary_vertices_positions )
    ## Mesh vertices should be unique (must they be, or only positions of boundary vertices?)
    assert len(set([ tuple(v) for v in mesh.vs ])) == len( mesh.vs )
    
    #debugger()
    
    #gen_laplacian, compute_target_curvatures_func, positions4laplacians = gen_cotangent_laplacian_matrix, compute_target_curvatures_mass_laplacian, optimize_positions_for_mass_laplacians
    #gen_laplacian, compute_target_curvatures_func, positions4laplacians = gen_cotangent_laplacian_matrix, compute_target_curvatures_mass_laplacian, optimize_positions_for_mass_laplacians_silhouette
    #gen_laplacian, compute_target_curvatures_func, positions4laplacians = gen_weighted_graph_laplacian_matrix, compute_target_curvatures_mass_laplacian, optimize_positions_for_mass_laplacians
    #gen_laplacian, compute_target_curvatures_func, positions4laplacians = gen_weighted_graph_laplacian_matrix, compute_target_curvatures_mass_laplacian, optimize_positions_for_mass_laplacians_silhouette
    #gen_laplacian, compute_target_curvatures_func, positions4laplacians = gen_graph_laplacian_matrix, compute_target_curvatures_paper, optimize_positions_for_unweighted_laplacians
    
    compute_target_curvatures_func = compute_target_curvatures_mass_laplacian
    ## Cotangent weights makes some extremely long, narrow faces that eventually lead to divergence.
    #gen_laplacian = gen_cotangent_laplacian_matrix
    gen_laplacian = gen_weighted_graph_laplacian_matrix
    positions4laplacians = optimize_positions_for_mass_laplacians
    #positions4laplacians = optimize_positions_for_mass_laplacians_silhouette
    
    
    ## Takeo's FiberMesh seems to use curvature constraint weight .1 and position constraint weight 10000.
    w_curvature = .1
    w_position = 10000.
    w_position = 1e7
    
    ## We can re-use the sparsity pattern of L
    N = len( mesh.vs )
    #L = assembly_matrix_t( (N,N) )
    L = None
    
    last_vs = arrayf( mesh.vs )
    for iter in range(iterations):
        
        ## 5b
        #L = assembly_matrix_t( (N,N) )
        L, vert_areas = gen_laplacian( mesh, L )
        ## The sparsity pattern won't change, so convert L into a csc matrix.
        ## UPDATE: Except that csc matrices don't support writable slices, which are used in gen_cotangent_laplacian_matrix()
        #if iter == 0: L = solve_matrix_t( L )
        
        #print 'L symmetry:', ( (L.T - L)**2 ).sum().sum()
        #print 'vert_areas:'
        #print vert_areas
        
        ## 6b
        ## This is done before calling this function.
        ## I'm not sure my curvature values are correct though.
        
        ## 7b
        ## Extend L with len(boundary_vertices) rows which keep values close to the values
        ## at boundary_vertices.
        #LMs = compute_target_curvatures_func( mesh, L, range(len(mesh.vs)), LMs )
        if iter == -1:
            LMs = compute_target_curvatures_func( mesh, L, boundary_vertices, asarray( boundary_vertices_curvatures ), w_curvature )
        else:
            ## This matches FiberMesh LaplacianSmoothing4:
            LMs = compute_target_curvatures_func( mesh, L, list(range(len(mesh.vs))), get_current_curvatures( mesh, L ), w_curvature )
        #print 'LMs:'
        #print LMs
        
        ## 8b no-op
        
        ## 9b
        #mesh.vs = optimize_positions_for_laplacians( mesh, L, vert_areas, LMs, boundary_vertices, boundary_vertices_positions, w_position )
        #mesh.vs = optimize_positions_for_laplacians2( mesh, L, vert_areas, LMs, boundary_vertices, boundary_vertices_positions, w_position )
        #mesh.vs = optimize_positions_for_unweighted_laplacians( mesh, L, vert_areas, LMs, boundary_vertices, boundary_vertices_positions, w_position )
        #mesh.vs = optimize_positions_for_mass_laplacians( mesh, L, vert_areas, LMs, boundary_vertices, boundary_vertices_positions, w_position )
        mesh.vs = positions4laplacians( mesh, L, vert_areas, LMs, boundary_vertices, boundary_vertices_positions, w_position )
        
        mesh.positions_changed()
        
        mesh.write_OBJ( name + 'optimized%s.obj' % (iter,) )
        print('max( | mesh.vs - last_vs | ):', sqrt( max( ( ( last_vs - asarrayf( mesh.vs ) )**2 ).sum(1) ) ))
        last_vs = arrayf( mesh.vs )
    
    return mesh

class Constraints( object ):
    def __init__( self ):
        ## Constraints on the xy positions of vertices.
        ## xy_constraint_vertices is the list of vertex indices, and
        ## xy_constraint_positions is the corresponding list of xy positions.
        self.xy_constraint_vertices = []
        self.xy_constraint_positions = []
        
        ## Constraints on the z values of vertices.
        ## z_constraint_vertices is the list of vertex indices, and
        ## z_constraint_values is the corresponding list of z values.
        self.z_constraint_vertices = []
        self.z_constraint_values = []
        
        ## Relative constraints on the z values of vertices.
        ## z_relative_constrained_groups_vertices is a list of groups of vertices constrained
        ## relative to one another.  Each group is a list of vertex indices.
        ## z_relative_constrained_groups_values is the corresponding list of values
        ## (a list for each group).
        self.z_relative_constrained_groups_vertices = []
        self.z_relative_constrained_groups_values = []
        
        ## Constraints on the weighted sum of vertex positions.
        ## xy_vertex_weighted_sum_indices is a list of groups of vertices.
        ## xy_vertex_weighted_sum_weights is a list of groups of values, corresponding to the weight
        ## of each vertex in the corresponding group in vertex_weighted_sum_indices.
        ## xy_vertex_weighted_sum_rhss is a list of right-hand-side values (2-vectors),
        ## one for each group.
        ## For each group, the constraint is \sum mesh.vs[ xy_indices[i] ].xy * xy_weights[i] = xy_rhs.
        self.xy_vertex_weighted_sum_indices = []
        self.xy_vertex_weighted_sum_weights = []
        self.xy_vertex_weighted_sum_rhss = []
        ## Similar to the xy_vertex_weighted_sum fields, except limited to z:
        ## For each group, the constraint is \sum mesh.vs[ z_indices[i] ].z * z_weights[i] = z_rhs.
        self.z_vertex_weighted_sum_indices = []
        self.z_vertex_weighted_sum_weights = []
        self.z_vertex_weighted_sum_rhss = []
        
        ## General linear constraints.
        ## general_linear_constraints is a list of linear constraints.
        ## Each linear constraint is a tuple of three items.
        ## The first item is a list of pairs specifying degrees of freedom involved in the
        ## equation: ( vertex index, coordinate offset ), where coordinate offset
        ## is one of 0, 1, or 2 specifying x, y, or z.
        ## The second item is a list of coefficients for each degree of freedom.
        ## The third item is the right-hand-side value.
        ## Precisely, each linear constraint is a constraint of the form:
        ##   \sum mesh.vs[ vertex_index[i] ][ coordinate_offset[i] ] * coefficients[i] = rhs
        self.general_linear_constraints = []
    
    
    def assert_self( self, mesh ):
        #PAINT()
        
        ## Make sure we have a value for every constrained vertex.
        assert len( self.xy_constraint_vertices ) == len( self.xy_constraint_positions )
        assert len( self.z_constraint_vertices ) == len( self.z_constraint_values )
        assert len( self.z_relative_constrained_groups_vertices ) == len( self.z_relative_constrained_groups_values )
        assert len( self.z_vertex_weighted_sum_indices ) == len( self.z_vertex_weighted_sum_weights )
        assert len( self.z_vertex_weighted_sum_indices ) == len( self.z_vertex_weighted_sum_rhss )
        assert False not in [ len( verts ) == len( vals ) for verts, vals in zip( self.z_relative_constrained_groups_vertices, self.z_relative_constrained_groups_values ) ]
        assert False not in [ len( verts ) == len( vals ) for verts, vals in zip( self.z_vertex_weighted_sum_indices, self.z_vertex_weighted_sum_weights ) ]
        ## Make sure we don't have doubly constrained vertices.
        assert len( set( self.xy_constraint_vertices ) ) == len( self.xy_constraint_vertices )
        assert len( set( self.z_constraint_vertices ) ) == len( self.z_constraint_vertices )
        assert False not in [ len( set( verts ) ) == len( verts ) for verts in self.z_relative_constrained_groups_vertices ]
        ## Users of this class should call general_linear_constraint_collect_terms() on each constraint.
        #assert False not in [ len( set( dofs ) ) == len( dofs ) for dofs, coeffs, rhs in self.general_linear_constraints ]
        
        ## Make sure our vertex indices are in-bounds.
        assert False not in [ v >= 0 and v < len( mesh.vs ) for v in self.xy_constraint_vertices + self.z_constraint_vertices ]
        assert False not in [ False not in [ v >= 0 and v < len( mesh.vs ) for v in verts ] for verts in self.z_relative_constrained_groups_vertices ]
        assert False not in [ False not in [ v >= 0 and v < len( mesh.vs ) for v in verts ] for verts in self.xy_vertex_weighted_sum_indices ]
        assert False not in [ False not in [ v >= 0 and v < len( mesh.vs ) for v in verts ] for verts in self.z_vertex_weighted_sum_indices ]
        
        ## We need at least one xy and one z constraint (or a general linear constraint)
        assert len( self.xy_constraint_vertices ) > 0 or len( self.xy_vertex_weighted_sum_indices ) > 0 or len( self.general_linear_constraints ) > 0
        assert len( self.z_constraint_vertices ) > 0 or len( self.z_vertex_weighted_sum_indices ) > 0 or len( self.general_linear_constraints ) > 0
    
    
    def general_linear_constraint_collect_terms( general_linear_constraint ):
        '''
        Given a general linear constraint as found in Constraints.general_linear_constraints,
        returns the same constraint with terms collapsed.
        That is, if a dof appears multiple times, in the parameter 'general_linear_constraint',
        it will appear once in the return value (with coefficient equal to the sum of
        the dof's coefficients in 'general_linear_constraint').
        '''
        
        assert len( general_linear_constraint[0] ) == len( general_linear_constraint[1] )
        
        dof2coeff = {}
        for dof, coeff in zip( general_linear_constraint[0], general_linear_constraint[1] ):
            dof2coeff.setdefault( dof, 0. )
            dof2coeff[ dof ] += coeff
        
        result = ( [], [], general_linear_constraint[2] )
        for dof, coeff in dof2coeff.items():
            result[0].append( dof )
            result[1].append( coeff )
        return result
    general_linear_constraint_collect_terms = staticmethod( general_linear_constraint_collect_terms )
    
    
    def mymath_z_hard_and_linear_constraints( self ):
        '''
        Returns two lists, hard constraints and linear constraints suitable for passing to
        mymath.solve_system_with_constraints(), corresponding to the z constraints specified
        by this object (self).  This does not include any coupled constraints such as in
        general_linear_constraints.
        '''
        
        #PAINT()
        
        linear_constraints = []
        for verts, relvals in zip( self.z_relative_constrained_groups_vertices, self.z_relative_constrained_groups_values ):
            assert len( verts ) == len( relvals )
            if len( verts ) == 0: continue
            
            ## Choose the first vertex in the group as "relative 0"
            base_vi = verts[0]
            base_z = relvals[0]
            ## For each vertex in the group, its index - base_vi must equal its zval - base_z
            for vert, zoff in zip( verts[1:], relvals[1:] ):
                linear_constraints.append( ( [vert, base_vi], [1., -1.], zoff - base_z ) )
        
        ## Handle z_vertex_weighted_sum constraints.
        for (indices, weights), zrhs in zip( zip( self.z_vertex_weighted_sum_indices, self.z_vertex_weighted_sum_weights ), self.z_vertex_weighted_sum_rhss ):
            #print 'z_vertex_weighted_sum:', ( indices, weights, zrhs )
            linear_constraints.append( ( indices, weights, zrhs ) )
        
        #print 'z constraints:'
        #pprint( zip( z_constraint_vertices, z_constraint_values ) )
        
        ## Move z constraints into linear_constraints
        #linear_constraints = [ ([index], [1.], val) for index, val in izip( z_constraint_vertices, z_constraint_values ) ]
        #z_constraint_vertices, z_constraint_values = [], []
        
        #print 'linear_constraints:'
        #pprint( linear_constraints )
        return list(zip( self.z_constraint_vertices, self.z_constraint_values )), linear_constraints
    
    def mymath_xy_hard_and_linear_constraints( self ):
        '''
        Returns two pairs of lists,
         ( x hard constraints, x linear constraints ), ( y hard constraints, y linear constraints ),
        suitable for passing to mymath.solve_system_with_constraints(),
        corresponding to the xy constraints specified by this object (self).
        This does not include any coupled constraints such as in general_linear_constraints.
        '''
        
        #PAINT()
        
        result = ( [], [] ), ( [], [] )
        
        ## Handle xy_vertex_weighted_sum constraints.
        for xycoord in (0,1):
            hard_constraints = result[ xycoord ][0]
            linear_constraints = result[ xycoord ][1]
            
            for (indices, weights), xyrhs in zip( zip( self.xy_vertex_weighted_sum_indices, self.xy_vertex_weighted_sum_weights ), self.xy_vertex_weighted_sum_rhss ):
                #print 'xy_vertex_weighted_sum:', ( indices, weights, xyrhs )
                linear_constraints.append( ( indices, weights, xyrhs[xycoord] ) )
            
            hard_constraints[:] = list(zip( self.xy_constraint_vertices, [ xy[xycoord] for xy in self.xy_constraint_positions ] ))
        
        return result
    
    def xyz_position_constraints( self ):
        '''
        Returns a list of elements ( index, xyz ) where the index
        is the vertex index and xyz is the location for the vertex.
        Call has_non_xyz_position_constraints() to check if there
        are any other kinds of constraints.
        '''
        
        indices = []
        xyzs = []
        
        index2xy = dict([
            ( index, xy )
            for index, xy
            in zip( self.xy_constraint_vertices, self.xy_constraint_positions )
            ])
        
        index2z = dict([
            ( index, z )
            for index, z
            in zip( self.z_constraint_vertices, self.z_constraint_values )
            ])
        
        common_indices = set( index2xy.keys() ).intersection( set( index2z.keys() ) )
        
        result = [
            ( index, tuple( index2xy[ index ] ) + ( index2z[ index ], ) )
            for index in common_indices
            ]
        return result
    
    def has_non_xyz_position_constraints( self ):
        '''
        Returns True if this object contains constraints that cannot be represented
        by the result of 'xyz_position_constraints()'.
        Returns False otherwise.
        '''
        
        if len( self.z_relative_constrained_groups_vertices ) > 0: return True
        if len( self.xy_vertex_weighted_sum_indices ) > 0: return True
        if len( self.z_vertex_weighted_sum_indices ) > 0: return True
        if len( self.general_linear_constraints ) > 0: return True
        
        if set( self.xy_constraint_vertices ) != set( self.z_constraint_vertices ): return True
        
        return False

def optimize_surface2( mesh, constraints, max_iterations = kDefaultIterations, absolute_convergence_threshold = 0., graph_not_cotangent_weights = True, average_edge_lengths = False, vertex_normal_motion_only = False, name = '' ):
    '''
    Given a TriMesh 'mesh', Constraints 'constraints',
    performs 'max_iterations' iterations of the non-linear inflation process,
    unless the maximum displacement of any vertex falls below 'absolute_convergence_threshold'.
    By default, 4 iterations are performed and 'absolute_convergence_threshold' is set to the
    impossible to reach value of 0.
    If 'graph_not_cotangent_weights' is False (the default is false), cotangent weights are used
    instead of graph laplacian weights.
    If 'average_edge_lengths' is True (the default is False), some edges are constrained to keep
    a smoothly varying length.
    If 'vertex_normal_motion_only' is True (the default is false), vertices are only allowed to move
    in their normal directions.
    
    Returns the maximum displacement of any vertex in the last iteration.
    
    NOTE: Modifies 'mesh' in place.
    '''
    
    from scipy import sparse
    
    ## Mesh vertices should be unique (must they be, or only positions of constraints?)
    assert len(set([ tuple(v) for v in mesh.vs ])) == len( mesh.vs )
    ## We want our mesh vertices to be a 2d array.
    mesh.vs = asarrayf( mesh.vs )
    
    #debugger()
    
    if graph_not_cotangent_weights:
        gen_laplacian = gen_weighted_graph_laplacian_matrix
    else:
        ## Cotangent weights makes some faces turn extremely long and narrow and eventually collapse,
        ## and the system diverges.
        ## And that problem doesn't go away if we constrain some edges (average_edge_lengths).
        ## It does go away if we constrain vertices to only move in the normal direction, though
        ## (vertex_normal_motion_only).
        gen_laplacian = gen_cotangent_laplacian_matrix
    
    compute_target_curvatures_func = compute_target_curvatures_mass_laplacian
    positions4laplacians = optimize_positions_for_mass_laplacians_complex_constraints
    
    ## Takeo's FiberMesh seems to use curvature constraint weight .1 and position constraint weight 10000.
    w_curvature = .1
    #w_position = 10000.
    ## I need this larger number for my linear constraints.
    w_position = 1e7
    
    ## We can re-use the sparsity pattern of L
    N = len( mesh.vs )
    #L = assembly_matrix_t( (N,N) )
    L = None
    
    last_vs = arrayf( mesh.vs )
    max_diff = absolute_convergence_threshold + 1.
    for iter in range(max_iterations):
        if max_diff < absolute_convergence_threshold: break
        
        ## 5b
        #L = assembly_matrix_t( (N,N) )
        L, vert_areas = gen_laplacian( mesh, L )
        L2 = L*L
        #print 'L nnz:', L.getnnz()
        #print 'L2 nnz:', L2.getnnz()
        ## The sparsity pattern won't change, so convert L into a csc matrix.
        ## UPDATE: Except that csc matrices don't support writable slices, which are used in gen_cotangent_laplacian_matrix()
        #if iter == 0: L = solve_matrix_t( L )
        
        #print 'L symmetry:', ( (L.T - L)**2 ).sum().sum()
        #print 'vert_areas:'
        #print vert_areas
        
        ## 6b
        ## This is done before calling this function.
        ## I'm not sure my curvature values are correct though.
        
        ## 7b
        ## Extend L with len(boundary_vertices) rows which keep values close to the values
        ## at boundary_vertices.
        #LMs = compute_target_curvatures_func( mesh, L, range(len(mesh.vs)), LMs )
        ## This matches FiberMesh LaplacianSmoothing4:
        LMs = compute_target_curvatures_func( mesh, L, L2, vert_areas, list(range(len(mesh.vs))), get_current_curvatures( mesh, L ), w_curvature )
        #print 'LMs:'
        #print LMs
        
        ## 8b no-op
        if average_edge_lengths:
            es = compute_target_curvatures_func( mesh, L, L2, vert_areas, list(range(len(mesh.vs))), get_current_average_edge_lengths( mesh ), w_curvature )
            old_xy_vertex_weighted_sum_indices = list( constraints.xy_vertex_weighted_sum_indices )
            old_xy_vertex_weighted_sum_weights = list( constraints.xy_vertex_weighted_sum_weights )
            old_xy_vertex_weighted_sum_rhss = list( constraints.xy_vertex_weighted_sum_rhss )
            old_z_vertex_weighted_sum_indices = list( constraints.z_vertex_weighted_sum_indices )
            old_z_vertex_weighted_sum_weights = list( constraints.z_vertex_weighted_sum_weights )
            old_z_vertex_weighted_sum_rhss = list( constraints.z_vertex_weighted_sum_rhss )
            for vi in constraints.xy_constraint_vertices:
                neighs = mesh.vertex_vertex_neighbors( vi )
                for vj in neighs:
                    ## Skip edges where both vertices are constrained.
                    if vj in constraints.xy_constraint_vertices: continue
                    
                    print('constraining edge ( %s, %s )' % ( vi, vj ))
                    
                    nij = ( es[ vi ] + es[ vj ] )/2. * dir( mesh.vs[ vi ] - mesh.vs[ vj ] )
                    #nij = 100 * dir( nij )
                    
                    constraints.xy_vertex_weighted_sum_indices.append( [ vi, vj ] )
                    constraints.xy_vertex_weighted_sum_weights.append( [ 1., -1 ] )
                    constraints.xy_vertex_weighted_sum_rhss.append( nij[:2] )
                    constraints.z_vertex_weighted_sum_indices.append( [ vi, vj ] )
                    constraints.z_vertex_weighted_sum_weights.append( [ 1., -1 ] )
                    constraints.z_vertex_weighted_sum_rhss.append( nij[2] )
        
        ## 9b
        vs = positions4laplacians( mesh, L, L2, LMs, constraints, w_position )
        if vertex_normal_motion_only:
            mesh.vs += ( ( vs - mesh.vs ) * mesh.vertex_normals ).sum(1).reshape((N,1)) * mesh.vertex_normals
        else:
            mesh.vs = vs
        
        if average_edge_lengths:
            constraints.xy_vertex_weighted_sum_indices = old_xy_vertex_weighted_sum_indices
            constraints.xy_vertex_weighted_sum_weights = old_xy_vertex_weighted_sum_weights
            constraints.xy_vertex_weighted_sum_rhss = old_xy_vertex_weighted_sum_rhss
            constraints.z_vertex_weighted_sum_indices = old_z_vertex_weighted_sum_indices
            constraints.z_vertex_weighted_sum_weights = old_z_vertex_weighted_sum_weights
            constraints.z_vertex_weighted_sum_rhss = old_z_vertex_weighted_sum_rhss
        
        mesh.positions_changed()
        
        mesh.write_OBJ( name + 'optimized%s.obj' % (iter,) )
        max_diff = sqrt( max( ( ( last_vs - mesh.vs )**2 ).sum(1) ) )
        print('max( | mesh.vs - last_vs | ):', max_diff)
        last_vs = arrayf( mesh.vs )
    
    return max_diff

def optimize_positions_for_laplacians_paper( mesh, L, areas, LMs, constraint_vertices, constraint_positions, constraint_weight = 1. ):
    '''
    Returns an len( mesh.vs ) x 3 array of new vertex positions.
    Parameters are
        'mesh', the input TriMesh,
        'L', the laplacian matrix corresponding to 'mesh',
        'areas', the area of each vertex,
        'LMs', the target curvature value for each mesh vertex,
        'constraint_vertices', the list of indices of constrained vertices,
        'constraint_positions', positions corresponding to the constrained vertices,
        (optional) 'constraint_weight', the scalar weight given to constraints (defaults to 1, the same as 'L')
    '''
    
    assert len( mesh.vs ) == L.shape[0]
    assert len( mesh.vs ) == L.shape[1]
    assert len( constraint_vertices ) == len( constraint_positions )
    
    assert False not in [ v >= 0 and v < len( mesh.vs ) for v in constraint_vertices ]
    
    assert constraint_weight >= 0.
    
    from scipy import sparse
    import scipy.sparse.linalg.dsolve as linsolve
    
    
    N = len( mesh.vs )
    BV = len( constraint_vertices )
    M = N + BV
    A = assembly_matrix_t( (M,N) )
    b = zerosf( (M,3) )
    for i in range( len( constraint_vertices ) ):
        A[ N + i, constraint_vertices[i] ] = constraint_weight
        b[ N + i, : ] = constraint_positions[i] * constraint_weight
    
    ## Put target laplacians at the top of b
    b[ :N, : ] = asarray( mesh.vertex_normals ) * ( areas * LMs ).reshape( ( LMs.shape[0], 1 ) )
    #b[ :N, : ] = asarray( mesh.vertex_normals ) * LMs.reshape( ( LMs.shape[0], 1 ) )
    
    ## Put the Laplacian matrix at the top of A
    A[:N,:N] = L
    
    rhs = A.T * b
    system = A.T * A
    assert rhs.shape == (N,3)
    assert system.shape == (N,N)
    
    system = solve_matrix_t( system ) #, dtype = 'F' )
    # linsolve.use_solver( use = {'useUmfpack' : True} )
    factorized = linsolve.factorized( system )
    
    vs = arrayf( mesh.vs ) #zerosf( rhs.shape )
    assert vs.shape == ( len( mesh.vs ), 3 )
    for i in range(3): vs[:,i] = factorized( rhs[:,i].flatten() )
    for vi, pos in zip( constraint_vertices, constraint_positions ): vs[ vi ] = pos
    
    #vs = [ xyz for xyz in mesh.vs ]
    return vs


def optimize_positions_for_unweighted_laplacians( mesh, L, areas, LMs, constraint_vertices, constraint_positions, constraint_weight = 1. ):
    '''
    Returns an len( mesh.vs ) x 3 array of new vertex positions.
    Parameters are
        'mesh', the input TriMesh,
        'L', the laplacian matrix corresponding to 'mesh',
        'areas', the area of each vertex,
        'LMs', the target curvature value for each mesh vertex,
        'constraint_vertices', the list of indices of constrained vertices,
        'constraint_positions', positions corresponding to the constrained vertices,
        (optional) 'constraint_weight', the scalar weight given to constraints (defaults to 1, the same as 'L')
    '''
    
    assert len( mesh.vs ) == L.shape[0]
    assert len( mesh.vs ) == L.shape[1]
    assert len( constraint_vertices ) == len( constraint_positions )
    
    assert False not in [ v >= 0 and v < len( mesh.vs ) for v in constraint_vertices ]
    
    assert constraint_weight >= 0.
    
    from scipy import sparse
    import scipy.sparse.linalg.dsolve as linsolve
    
    
    N = len( mesh.vs )
    ## Put the Laplacian matrix at the top of A
    A = L.T * L
    #A = L * L
    ## Put target laplacians at the top of b
    b = asarray( mesh.vertex_normals ) * ( areas * LMs ).reshape( ( LMs.shape[0], 1 ) )
    #b = asarray( mesh.vertex_normals ) * LMs.reshape( ( LMs.shape[0], 1 ) )
    #b = L * asmatrix( asarray( mesh.vertex_normals ) * LMs.reshape( ( LMs.shape[0], 1 ) ) )
    
    ## Add hard constraints
    A = assembly_matrix_t( A )
    for i in range( len( constraint_vertices ) ):
        A[ constraint_vertices[i], : ] = 0.
        A[ constraint_vertices[i], constraint_vertices[i] ] = 1.
        b[ constraint_vertices[i], : ] = constraint_positions[i]
    
    rhs = asarray( b )
    system = A
    assert rhs.shape == (N,3)
    assert system.shape == (N,N)
    
    system = solve_matrix_t( system ) #, dtype = 'F' )
    # linsolve.use_solver( use = {'useUmfpack' : True} )
    factorized = linsolve.factorized( system )
    
    vs = arrayf( mesh.vs ) #zerosf( rhs.shape )
    assert vs.shape == ( len( mesh.vs ), 3 )
    for i in range(3): vs[:,i] = factorized( rhs[:,i].flatten() )
    assert False not in [ mag2( vs[ vi ] - pos ) < 1e-5 for vi, pos in zip( constraint_vertices, constraint_positions ) ]
    #for vi, pos in izip( constraint_vertices, constraint_positions ): vs[ vi ] = pos
    
    #vs = [ xyz for xyz in mesh.vs ]
    return vs


def optimize_positions_for_mass_laplacians( mesh, L, areas, LMs, constraint_vertices, constraint_positions, constraint_weight = 1. ):
    '''
    Returns an len( mesh.vs ) x 3 array of new vertex positions.
    Parameters are
        'mesh', the input TriMesh,
        'L', the laplacian matrix corresponding to 'mesh',
        'areas', the area of each vertex (unused),
        'LMs', the target curvature value for each mesh vertex,
        'constraint_vertices', the list of indices of constrained vertices,
        'constraint_positions', positions corresponding to the constrained vertices,
        (optional) 'constraint_weight', the scalar weight given to constraints (defaults to 1, the same as 'L')
    '''
    
    assert len( mesh.vs ) == L.shape[0]
    assert len( mesh.vs ) == L.shape[1]
    assert len( constraint_vertices ) == len( constraint_positions )
    
    assert False not in [ v >= 0 and v < len( mesh.vs ) for v in constraint_vertices ]
    
    assert constraint_weight >= 0.
    
    from scipy import sparse
    import scipy.sparse.linalg.dsolve as linsolve
    
    
    N = len( mesh.vs )
    ## Put the Laplacian matrix at the top of A
    #A = L.T * L
    A = L * L
    ## Put target laplacians at the top of b
    #b = asarray( mesh.vertex_normals ) * ( areas * LMs ).reshape( ( LMs.shape[0], 1 ) )
    #b = asarray( mesh.vertex_normals ) * LMs.reshape( ( LMs.shape[0], 1 ) )
    b = L * asmatrix( asarray( mesh.vertex_normals ) * LMs.reshape( ( LMs.shape[0], 1 ) ) )
    
    ## Add hard constraints
    A = assembly_matrix_t( A )
    for i in range( len( constraint_vertices ) ):
        A[ constraint_vertices[i], : ] = 0.
        A[ constraint_vertices[i], constraint_vertices[i] ] = 1.
        b[ constraint_vertices[i], : ] = constraint_positions[i]
    
    rhs = asarray( b )
    system = A
    assert rhs.shape == (N,3)
    assert system.shape == (N,N)
    
    system = solve_matrix_t( system ) #, dtype = 'F' )
    
    system, rhs = improve_condition( system, rhs )
    
    # linsolve.use_solver( use = {'useUmfpack' : True} )
    factorized = linsolve.factorized( system )
    
    vs = arrayf( mesh.vs ) #zerosf( rhs.shape )
    assert vs.shape == ( len( mesh.vs ), 3 )
    for i in range(3): vs[:,i] = factorized( rhs[:,i].flatten() )
    assert False not in [ mag2( vs[ vi ] - pos ) < 1e-5 for vi, pos in zip( constraint_vertices, constraint_positions ) ]
    #expect( False not in [ mag2( vs[ vi ] - pos ) < 1e-5 for vi, pos in izip( constraint_vertices, constraint_positions ) ] )
    for vi, pos in zip( constraint_vertices, constraint_positions ): vs[ vi ] = pos
    
    #vs = [ xyz for xyz in mesh.vs ]
    return vs

def optimize_positions_for_mass_laplacians_silhouette( mesh, L, areas, LMs, constraint_vertices, constraint_positions, constraint_weight = 1. ):
    ## Turn positional constraints into silhouette constraints (xy only, not z)
    xy_constraint_vertices = list( constraint_vertices[1:] )
    xy_constraint_positions = list( constraint_positions[1:] )
    
    ## Turn the first constraint into a full xyz constraint.
    ## This works fine, and produces the same result.
    xyz_constraint_vertices = [ constraint_vertices[0] ]
    xyz_constraint_positions = [ constraint_positions[0] ]
    
    ## As another test, add a second xyz constraint but fix its z value to be different.
    ## For the second constraints, choose the point with the farthest (in xy) distance from
    ## out first xyz constraint.
    from helpers import farthest_point_index
    two = farthest_point_index( xyz_constraint_positions[0][:2], asarray( xy_constraint_positions )[:,:2] )
    xyz_constraint_vertices.append( xy_constraint_vertices[ two ] )
    xyz_constraint_positions.append( xy_constraint_positions[ two ] )
    del xy_constraint_vertices[ two ]
    del xy_constraint_positions[ two ]
    ## Fix z's of the two constraints to be different.
    xyz_constraint_positions[0][2] = 0.
    xyz_constraint_positions[1][2] = -1.
    
    return optimize_positions_for_mass_laplacians_silhouette_real( mesh, L, areas, LMs, xy_constraint_vertices, xy_constraint_positions, xyz_constraint_vertices, xyz_constraint_positions, constraint_weight )

def optimize_positions_for_mass_laplacians_silhouette_real( mesh, L, areas, LMs, xy_constraint_vertices, xy_constraint_positions, xyz_constraint_vertices, xyz_constraint_positions, constraint_weight = 1. ):
    '''
    Returns an len( mesh.vs ) x 3 array of new vertex positions.
    Parameters are
        'mesh', the input TriMesh,
        'L', the laplacian matrix corresponding to 'mesh',
        'areas', the area of each vertex (unused),
        'LMs', the target curvature value for each mesh vertex,
        'xy_constraint_vertices', the list of indices of xy constrained vertices,
        'xy_constraint_positions', xy positions corresponding to the xy constrained vertices,
        'xyz_constraint_vertices', the list of indices of xyz constrained vertices,
        'xyz_constraint_positions', xyz positions corresponding to the xyz constrained vertices,
        (optional) 'constraint_weight', the scalar weight given to constraints (defaults to 1, the same as 'L')
    '''
    
    assert len( mesh.vs ) == L.shape[0]
    assert len( mesh.vs ) == L.shape[1]
    assert len( xy_constraint_vertices ) == len( xy_constraint_positions )
    assert len( xyz_constraint_vertices ) == len( xyz_constraint_positions )
    
    assert False not in [ v >= 0 and v < len( mesh.vs ) for v in xy_constraint_vertices + xyz_constraint_vertices ]
    assert len( set( xy_constraint_vertices + xyz_constraint_vertices ) ) == len( xy_constraint_vertices ) + len( xyz_constraint_vertices )
    
    assert constraint_weight >= 0.
    
    from scipy import sparse
    import scipy.sparse.linalg.dsolve as linsolve
    
    
    N = len( mesh.vs )
    ## Put the Laplacian matrix at the top of A
    #A = L.T * L
    A = L * L
    ## Put target laplacians at the top of b
    #b = asarray( mesh.vertex_normals ) * ( areas * LMs ).reshape( ( LMs.shape[0], 1 ) )
    #b = asarray( mesh.vertex_normals ) * LMs.reshape( ( LMs.shape[0], 1 ) )
    b = L * asmatrix( asarray( mesh.vertex_normals ) * LMs.reshape( ( LMs.shape[0], 1 ) ) )
    
    ## Add hard constraints
    Axy = assembly_matrix_t( A, copy = True )
    for i in range( len( xy_constraint_vertices ) ):
        Axy[ xy_constraint_vertices[i], : ] = 0.
        Axy[ xy_constraint_vertices[i], xy_constraint_vertices[i] ] = 1.
        b[ xy_constraint_vertices[i], :2 ] = xy_constraint_positions[i][:2]
    
    Az = assembly_matrix_t( A, copy = True )
    for i in range( len( xyz_constraint_vertices ) ):
        Axy[ xyz_constraint_vertices[i], : ] = 0.
        Axy[ xyz_constraint_vertices[i], xyz_constraint_vertices[i] ] = 1.
        
        Az[ xyz_constraint_vertices[i], : ] = 0.
        Az[ xyz_constraint_vertices[i], xyz_constraint_vertices[i] ] = 1.
        
        b[ xyz_constraint_vertices[i], : ] = xyz_constraint_positions[i]
    
    rhs = asarray( b )
    systemxy = Axy
    systemz = Az
    assert rhs.shape == (N,3)
    assert systemxy.shape == (N,N)
    assert systemz.shape == (N,N)
    
    systemxy = solve_matrix_t( systemxy ) #, dtype = 'F' )
    systemz = solve_matrix_t( systemz ) #, dtype = 'F' )
    
    # linsolve.use_solver( use = {'useUmfpack' : True} )
    factorizedxy = linsolve.factorized( systemxy )
    factorizedz = linsolve.factorized( systemz )
    
    vs = arrayf( mesh.vs ) #zerosf( rhs.shape )
    assert vs.shape == ( len( mesh.vs ), 3 )
    for i in range(2): vs[:,i] = factorizedxy( rhs[:,i].flatten() )
    for i in range(2,3): vs[:,i] = factorizedz( rhs[:,i].flatten() )
    assert False not in [ mag2( vs[ vi ][:2] - pos[:2] ) < 1e-5 for vi, pos in zip( xy_constraint_vertices, xy_constraint_positions ) ]
    assert False not in [ mag2( vs[ vi ] - pos ) < 1e-5 for vi, pos in zip( xyz_constraint_vertices, xyz_constraint_positions ) ]
    for vi, pos in zip( xy_constraint_vertices, xy_constraint_positions ): vs[ vi ][:2] = pos[:2]
    for vi, pos in zip( xyz_constraint_vertices, xyz_constraint_positions ): vs[ vi ] = pos
    
    #vs = [ xyz for xyz in mesh.vs ]
    return vs

def optimize_positions_for_mass_laplacians_complex_constraints( mesh, L, L2, LMs, constraints, constraint_weight = 1. ):
    '''
    Returns an len( mesh.vs ) x 3 array of new vertex positions.
    Parameters are
        'mesh', the input TriMesh,
        'L', the laplacian matrix corresponding to 'mesh',
        'LMs', the target curvature value for each mesh vertex,
        'constraints' containing the following fields:
            'xy_constraint_vertices', the list of indices of xy constrained vertices,
            'xy_constraint_positions', xy positions corresponding to the xy constrained vertices,
            'z_constraint_vertices', the list of indices of z constrained vertices,
            'z_constraint_values', z values corresponding to the z constrained vertices,
            'z_relative_constrained_groups_vertices', the list of lists of indices of vertices, each sublist's vertices are relatively constraint in z,
            'z_relative_constrained_groups_values', the list of lists of z values, each sublist's values correspond to the relative z values of constrained vertices in the group,
            and more!  See comments in 'class Constraints'
        (optional) 'constraint_weight', the scalar weight given to constraints (defaults to 1, the same as 'L')
    '''
    
    #debugger()
    
    # return optimize_positions_for_mass_laplacians( mesh, L, ones(len(mesh.vs)), LMs, constraints.xy_constraint_vertices, array( [ (x,y,0.) for x,y in constraints.xy_constraint_positions ] ), constraint_weight )
    
    xy_constraint_vertices = constraints.xy_constraint_vertices
    xy_constraint_positions = constraints.xy_constraint_positions
    z_constraint_vertices = constraints.z_constraint_vertices
    z_constraint_values = constraints.z_constraint_values
    z_relative_constrained_groups_vertices = constraints.z_relative_constrained_groups_vertices
    z_relative_constrained_groups_values = constraints.z_relative_constrained_groups_values
    xy_vertex_weighted_sum_indices = constraints.xy_vertex_weighted_sum_indices
    xy_vertex_weighted_sum_weights = constraints.xy_vertex_weighted_sum_weights
    xy_vertex_weighted_sum_rhss = constraints.xy_vertex_weighted_sum_rhss
    z_vertex_weighted_sum_indices = constraints.z_vertex_weighted_sum_indices
    z_vertex_weighted_sum_weights = constraints.z_vertex_weighted_sum_weights
    z_vertex_weighted_sum_rhss = constraints.z_vertex_weighted_sum_rhss
    general_linear_constraints = constraints.general_linear_constraints
    
    #print 'constraints:'
    #pprint( constraints.__dict__ )
    #print 'zrel constraints:'
    #for verts, vals in izip( z_relative_constrained_groups_vertices, z_relative_constrained_groups_values ): pprint( zip( verts, vals ) )
    
    ## Make sure our dimensions are correct.
    assert len( mesh.vs ) == L.shape[0]
    assert len( mesh.vs ) == L.shape[1]
    assert L.shape == L2.shape
    constraints.assert_self( mesh )
    
    assert constraint_weight >= 0.
    
    from scipy import sparse
    import scipy.sparse.linalg.dsolve as linsolve
    
    vs = arrayf( mesh.vs ) #zerosf( rhs.shape )
    assert vs.shape == ( len( mesh.vs ), 3 )
    
    N = len( mesh.vs )
    ## Put the Laplacian matrix at the top of A
    #A = L.T * L
    #A = L * L
    A = L2.copy()
    ## Put target laplacians at the top of b
    #b = asarray( mesh.vertex_normals ) * ( areas * LMs ).reshape( ( LMs.shape[0], 1 ) )
    #b = asarray( mesh.vertex_normals ) * LMs.reshape( ( LMs.shape[0], 1 ) )
    b = L * ( asmatrix( asarray( mesh.vertex_normals ) * LMs.reshape( ( LMs.shape[0], 1 ) ) ) )
    
    import mymath
    xyz_hard_linear_constraints = list( constraints.mymath_xy_hard_and_linear_constraints() )
    xyz_hard_linear_constraints.append( constraints.mymath_z_hard_and_linear_constraints() )
    
    print('xyz_hard_linear_constraints:')
    pprint( xyz_hard_linear_constraints )
    #debugger()
    
    print('general_linear_constraints:')
    pprint( general_linear_constraints )
    
    if len( constraints.general_linear_constraints ) > 0:
        #PAINT()
        
        ## If we have general_linear_constraints, we're going to assume that the system is coupled.
        ## So, put all degrees of freedom into a 3N system.
        
        A3 = assembly_matrix_t( ( 3*N, 3*N ), dtype = A.dtype )
        A3[:N,:N] = A
        A3[N:2*N,N:2*N] = A
        A3[2*N:3*N,2*N:3*N] = A
        
        b3 = empty( 3*N, b.dtype )
        b3[:N] = b[:,0].ravel()
        b3[N:2*N] = b[:,1].ravel()
        b3[2*N:3*N] = b[:,2].ravel()
        
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
        
        solution3 = mymath.solve_system_with_constraints(
            A3, b3,
            hard_constraints,
            linear_constraints,
            'lsq',
            constraint_weight
            )
        
        for i in range(N): vs[i,:] = solution3[ i ], solution3[ i + N ], solution3[ i + 2*N ]
        
    else:
        #PAINT()
        
        rhs3 = asarray( b )
        assert rhs3.shape == (N,3)
        
        coords = [0,1,2]
        ## Can we handle x and y together by directly constraining the matrix?
        ## This is a common case (for example, when we have only region_boundary constraints).
        if( len( xyz_hard_linear_constraints[0][1] ) == 0 and
            len( xyz_hard_linear_constraints[1][1] ) == 0 and
            ## x and y hard constraints operate on the same vertices
            [ vi for vi, rhs in xyz_hard_linear_constraints[0][0] ] == [ vi for vi, rhs in xyz_hard_linear_constraints[1][0] ]
            ):
            #PAINT()
            
            ## Add hard xy constraints.
            assert len( constraints.xy_vertex_weighted_sum_indices ) == 0
            Axy = assembly_matrix_t( A, copy = True )
            for i in range(len( xyz_hard_linear_constraints[0][0] )):
                vi = xyz_hard_linear_constraints[0][0][i][0]
                ## Scipy 0.8 (or maybe 0.7+)'s lil matrices no longer allow setting a slice from a single scalar value.
                #Axy[ vi, : ] = 0.
                Axy[ vi, : ] *= 0.
                Axy[ vi, vi ] = 1.
                b[ vi, :2 ] = xyz_hard_linear_constraints[0][0][i][1], xyz_hard_linear_constraints[1][0][i][1]
            
            systemxy = Axy
            assert systemxy.shape == (N,N)
            systemxy = solve_matrix_t( systemxy ) #, dtype = 'F' )
            # linsolve.use_solver( use = {'useUmfpack' : True} )
            factorizedxy = linsolve.factorized( systemxy )
            
            ## Solve for new xy positions.
            for i in range(2): vs[:,i] = factorizedxy( rhs3[:,i].flatten() )
            #assert False not in [ mag2( vs[ vi ][:2] - asarray((posx,posy)) ) < 1e-5 for (vi, posx), (viy, posy) in izip( xyz_hard_linear_constraints[0][0], xyz_hard_linear_constraints[1][0] ) ]
            #for vi, pos in izip( xy_constraint_vertices, xy_constraint_positions ): vs[ vi ][:2] = pos[:2]
            
            ## We've taken care of x and y, all that's left is z.
            coords = [2]
        
        ## The uncoupled general case.
        for coord in coords:
            #PAINT()
            vs[:,coord] = mymath.solve_system_with_constraints(
                A, rhs3[:,coord].flatten(),
                xyz_hard_linear_constraints[coord][0],
                xyz_hard_linear_constraints[coord][1],
                'lsq',
                constraint_weight
                )
    
    ## Verify that our hard constraints were met.
    for coord in (0,1,2):
        assert False not in [ ( vs[ vi ][ coord ] - pos )**2 < 1e-5 for vi, pos in xyz_hard_linear_constraints[coord][0] ]
        #for vi, pos in xyz_hard_linear_constraints[i][0]: vs[ vi, coord ] = pos
    
    ## Pass through vertices unchanged, to debug something else.
    #vs = [ xyz for xyz in mesh.vs ]
    return vs

def compute_target_curvatures_paper( mesh, L, constraint_vertices, constraint_curvatures, constraint_weight = 1. ):
    '''
    Returns an array of scalar curvature values, one for each vertex in the mesh.
    Parameters are
        'mesh', the input TriMesh,
        'L', the laplacian matrix corresponding to 'mesh',
        'constraint_vertices', the list of indices of constrained vertices,
        'constraint_curvatures', curvatures corresponding to the constrained vertices,
        (optional) 'constraint_weight', the scalar weight given to constraints (defaults to 1, the same as 'L')
    '''
    
    assert len( mesh.vs ) == L.shape[0]
    assert len( mesh.vs ) == L.shape[1]
    assert len( constraint_vertices ) == len( constraint_curvatures )
    
    assert False not in [ v >= 0 and v < len( mesh.vs ) for v in constraint_vertices ]
    
    assert constraint_weight >= 0.
    
    from scipy import sparse
    import scipy.sparse.linalg.dsolve as linsolve
    
    
    ## Clip curvatures to 0.
    print('negative constraint_curvatures:', constraint_curvatures[ constraint_curvatures < 0. ])
    if kClipNegativeCurvatureConstraints: constraint_curvatures[ constraint_curvatures < 0. ] = 0.
    
    N = len( mesh.vs )
    BV = len( constraint_vertices )
    M = N + BV
    A = assembly_matrix_t( (M,N) )
    b = zerosf( M )
    for i in range( len( constraint_vertices ) ):
        A[ N + i, constraint_vertices[i] ] = constraint_weight
        b[ N + i ] = constraint_curvatures[i] * constraint_weight
    
    ## Put the Laplacian matrix at the top
    A[:N,:N] = L
    rhs = A.T * b
    assert rhs.shape == (N,)
    system = A.T * A
    assert system.shape == (N,N)
    system = solve_matrix_t( system ) #, dtype = 'F' )
    # linsolve.use_solver( use = {'useUmfpack' : True} )
    factorized = linsolve.factorized( system )
    LMs = factorized( rhs )
    
    ## Clip curvatures to 0.
    print('negative LMs:', LMs[ LMs < 0. ])
    LMs[ LMs < 0. ] = 0.
    
    return LMs


def compute_target_curvatures_mass_laplacian( mesh, L, L2, areas, constraint_vertices, constraint_curvatures, constraint_weight = 1. ):
    '''
    Returns an array of scalar curvature values, one for each vertex in the mesh.
    Parameters are
        'mesh', the input TriMesh,
        'L', the laplacian matrix corresponding to 'mesh',
        'constraint_vertices', the list of indices of constrained vertices,
        'constraint_curvatures', curvatures corresponding to the constrained vertices,
        (optional) 'constraint_weight', the scalar weight given to constraints (defaults to 1, the same as 'L')
    '''
    
    assert len( mesh.vs ) == L.shape[0]
    assert len( mesh.vs ) == L.shape[1]
    assert L.shape == L2.shape
    assert len( constraint_vertices ) == len( constraint_curvatures )
    
    assert False not in [ v >= 0 and v < len( mesh.vs ) for v in constraint_vertices ]
    
    assert constraint_weight >= 0.
    
    from scipy import sparse
    import scipy.sparse.linalg.dsolve as linsolve
    
    
    ## Clip curvatures to 0.
    print('negative constraint_curvatures:', constraint_curvatures[ constraint_curvatures < 0. ])
    if kClipNegativeCurvatureConstraints: constraint_curvatures[ constraint_curvatures < 0. ] = 0.
    
    N = len( mesh.vs )
    BV = len( constraint_vertices )
    
    kReuseFactorization = False
    if kReuseFactorization:
        from scipy import sparse
        diagonal_matrix = sparse.identity(N).tocsr()
        #for i in xrange(N): diagonal_matrix[i,i] = 1./areas[i]
        diagonal_matrix.setdiag( areas )
        A = diagonal_matrix * L.tocsr()
    else:
        #A = L * L
        A = L2.copy()
    
    assert A.shape == (N,N)
    #b = zerosf( N )
    
    C = assembly_matrix_t( (BV, N) )
    C_rhs = zerosf( BV )
    for i in range( len( constraint_vertices ) ):
        C[ i, constraint_vertices[i] ] = 1.
        C_rhs[ i ] = constraint_curvatures[i]
    
    A = A.tocsr()
    C = C.tocsr()
    system = A + constraint_weight * (C.T*C)
    rhs = constraint_weight * asarray( C.T * C_rhs ) # + b
    assert rhs.shape == (N,)
    assert system.shape == (N,N)
    system = solve_matrix_t( system ) #, dtype = 'F' )
    
    system, rhs = improve_condition( system, rhs )
    
    # linsolve.use_solver( use = {'useUmfpack' : True} )
    #try:
    factorized = linsolve.factorized( system )
    if kReuseFactorization:
        LMs = factorized( diagonal_matrix.tocsr() * rhs )
        LMs = factorized( diagonal_matrix.tocsr() * LMs + constraint_weight * asarray( C.T * C_rhs ) )
    else:
        LMs = factorized( rhs )
    #except RuntimeError, val:
    #    print 'RuntimeError:', val
    #    return zeros( rhs.shape )
    
    ## Clip curvatures to 0.
    print('negative LMs:', LMs[ LMs < 0. ])
    LMs[ LMs < 0. ] = 0.
    
    return LMs


def gen_regular_grid( **kwargs ):
    '''
    Generates a regular grid of vertices on the xy plane.
    Parameters are xmin, xmax, ymin, ymax, xverts, and yverts.
    By default, the grid has 5 vertices in both the x and y directions,
    and goes from [0..1]x[0..1].
    Returns a mesh_gui.TriMesh.
    '''
    
    kwargs.setdefault( 'xmin', 0. )
    kwargs.setdefault( 'xmax', 1. )
    kwargs.setdefault( 'ymin', 0. )
    kwargs.setdefault( 'ymax', 1. )
    
    kwargs.setdefault( 'xverts', 5 )
    kwargs.setdefault( 'yverts', 5 )
    
    xmin = kwargs['xmin']
    xscale = kwargs['xmax'] - xmin
    ymin = kwargs['ymin']
    yscale = kwargs['ymax'] - ymin
    
    xverts = kwargs['xverts']
    yverts = kwargs['yverts']
    
    from mesh_gui import TriMesh
    mesh = TriMesh()
    for yi in range( yverts ):
        for xi in range( xverts ):
            mesh.vs.append( arrayf( ( xmin + xscale * xi/float(xverts-1), ymin + yscale * yi/float(yverts-1), 0. ) ) )
            
            if xi > 0 and yi > 0:
                ## alternating regular grid
                if len( mesh.vs ) % 2:
                    mesh.faces.append( (
                        (xi-1) + (yi-1)*xverts,
                        (xi-0) + (yi-1)*xverts,
                        (xi-0) + (yi-0)*xverts
                        ) )
                    mesh.faces.append( (
                        (xi-1) + (yi-1)*xverts,
                        (xi-0) + (yi-0)*xverts,
                        (xi-1) + (yi-0)*xverts
                        ) )
                else:
                    mesh.faces.append( (
                        (xi-1) + (yi-1)*xverts,
                        (xi-0) + (yi-1)*xverts,
                        (xi-1) + (yi-0)*xverts
                        ) )
                    mesh.faces.append( (
                        (xi-0) + (yi-1)*xverts,
                        (xi-0) + (yi-0)*xverts,
                        (xi-1) + (yi-0)*xverts
                        ) )
    
    mesh.topology_changed()
    return mesh

def gen_equilateral_grid( **kwargs ):
    '''
    Generates an equilateral grid of vertices on the xy plane.
    Parameters are xmin, xmax, ymin, ymax, xverts, and yverts.
    By default, the grid has 5 vertices in both the x and y directions,
    and goes from [0..1]x[0..1].
    Returns a mesh_gui.TriMesh.
    '''
    
    print('WARNING: This does not obey x,y min/max parameters and so sample the stroke all wrong.', file=sys.stderr)
    
    kwargs.setdefault( 'xmin', 0. )
    kwargs.setdefault( 'xmax', 1. )
    kwargs.setdefault( 'ymin', 0. )
    kwargs.setdefault( 'ymax', 1. )
    
    kwargs.setdefault( 'xverts', 5 )
    kwargs.setdefault( 'yverts', 5 )
    
    xmin = kwargs['xmin']
    xside = kwargs['xmax'] - xmin
    ymin = kwargs['ymin']
    yside = kwargs['ymax'] - ymin
    
    xverts = kwargs['xverts']
    yverts = kwargs['yverts']
    
    from generate_rectangular import generate_equilateral_points, generate_equilateral_triangles
    
    side_length = min( float(xside)/xverts, float(yside)/yverts )
    grid_points, initial_points = generate_equilateral_points(
        xmin, ymin, side_length, yverts-1, xverts-1,
        lambda x,y: (x,y,0.), lambda x,y: (x,y,0.), lambda x,y: (x,y,0.),
        0
        )
    triangles = generate_equilateral_triangles( yverts-1, xverts-1 )
    
    from mesh_gui import TriMesh
    mesh = TriMesh()
    mesh.vs = [ asarrayf( pt ) for pt in grid_points ]
    mesh.faces = [ (i,k,j) for i,j,k in triangles ]
    
    mesh.topology_changed()
    return mesh

def get_current_curvatures( mesh, L ):
    '''
    Returns a vector of length len( mesh.vs ) containing the scalar curvature of each mesh
    vertex.
    '''
    
    assert L.shape[0] == L.shape[1]
    assert L.shape[0] == len( mesh.vs )
    
    laplacians = L * asmatrix( mesh.vs )
    assert laplacians.shape == ( len( mesh.vs ), 3 )
    
    curvatures = ( asarray( laplacians ) * asarray( mesh.vertex_normals ) ).sum( 1 )
    #curvatures = sqrt( ( asarray( laplacians ) * asarray( laplacians ) ).sum( 1 ) )
    #curvatures = sqrt( ( asarray( laplacians ) * asarray( laplacians ) ).sum( 1 ) ) / 2.
    assert curvatures.shape[0] == len( mesh.vs )
    assert prod( curvatures.shape ) == len( mesh.vs )
    
    return curvatures.flatten()

def get_current_average_edge_lengths( mesh ):
    '''
    Returns a vector of length len( mesh.vs ) containing the average incident edge length of each
    mesh vertex.
    '''
    
    result = emptyf( len( mesh.vs ) )
    for vi in range( len( mesh.vs ) ):
        neighs = mesh.vertex_vertex_neighbors( vi )
        print('neighs:', neighs)
        result[vi] = sqrt( ( ( mesh.vs[vi] - mesh.vs[ neighs ] )**2 ).sum(1) ).sum() / len( neighs )
    return result

def gen_graph_laplacian_matrix( mesh, A ):
    '''
    Places the graph (unweighted) Laplacian operator matrix for 'mesh' into the upper-left of 'A'.
    Thus, 'A' must have at least len(mesh.vs) rows and columns.
    NOTE: The Laplacian operator matrix contains one row per vertex, not three.
    '''
    
    N = len( mesh.vs )
    if A is None: A = assembly_matrix_t( (N,N) )
    assert A.shape[0] >= N and A.shape[1] >= N
    
    for vi in range( N ):
        neighbor_indices = mesh.vertex_vertex_neighbors( vi )
        neigh_val = -1./len(neighbor_indices)
        A[ vi, vi ] = 1.
        
        A[ vi, tuple( neighbor_indices ) ] = neigh_val
        ## csc matrices don't support writable slices
        #for ni in neighbor_indices: A[ vi, ni ] = neigh_val
    
    return A, asarrayf( mesh.vertex_areas )

__graph_laplacian_matrix_cache = {}
def gen_weighted_graph_laplacian_matrix( mesh, A ):
    '''
    Places the graph (unweighted) Laplacian operator matrix for 'mesh' into the upper-left of 'A'.
    Thus, 'A' must have at least len(mesh.vs) rows and columns.
    NOTE: The Laplacian operator matrix contains one row per vertex, not three.
    '''
    
    #print 'hash( mesh ):', hash( mesh )
    
    N = len( mesh.vs )
    
    if mesh in __graph_laplacian_matrix_cache:
        print('cached!')
        A = __graph_laplacian_matrix_cache[ mesh ]
    else:
        A = assembly_matrix_t( (N,N) )
        assert A.shape[0] >= N and A.shape[1] >= N
        
        for vi in range( N ):
            w = 1.
            neighbor_indices = mesh.vertex_vertex_neighbors( vi )
            neigh_val = -w/len(neighbor_indices)
            A[ vi, vi ] = w
            
            ## Scipy 0.8 (or maybe 0.7+)'s lil matrices no longer allow setting a slice from a single scalar value.
            #A[ vi, tuple( neighbor_indices ) ] = neigh_val
            A[ vi, tuple( neighbor_indices ) ] = [neigh_val] * len( neighbor_indices )
            ## csc matrices don't support writable slices.
            #for ni in neighbor_indices: A[ vi, ni ] = neigh_val
        
        __graph_laplacian_matrix_cache[ mesh ] = A
    
    areas = asarrayf( mesh.vertex_areas )
    ## Diagonal matrix construction is funny.  See my notes in improve_condition().
    from scipy import sparse
    diagonal_matrix = sparse.identity(N).tocsr()
    #for i in xrange(N): diagonal_matrix[i,i] = 1./areas[i]
    diagonal_matrix.setdiag( 1./areas )
    A = diagonal_matrix * A.tocsr()
    
    '''
    ## The previous code
    Aold = assembly_matrix_t( (N,N) )
    for vi in xrange( N ):
        w = 1./areas[vi]
        neighbor_indices = mesh.vertex_vertex_neighbors( vi )
        neigh_val = -w/len(neighbor_indices)
        Aold[ vi, vi ] = w
        Aold[ vi, tuple( neighbor_indices ) ] = neigh_val
    
    print 'areas:', areas
    print '|A-Aold|^2:', ((A-Aold)**2).sum()
    A = Aold
    '''
    
    ## Scale A by the reciprocal of the average value of the diagonal, to that
    ## diagonals are, on average, close to 1.
    ## This has two effects:
    ## 1) The condition number of the system should be improved.
    ## 2) The magnitude of the energy using this matrix should be roughly normalized,
    ##    so that the weight of the least-squares constraint energy is meaningful.
    ## UPDATE: On second thought, scaling this matrix will take it *out* of scale
    ##         from other terms in the energy which are scale-dependent,
    ##         including position displacement.
    if kScaleLaplacianWeights:
        conditioner = float(N) / abs(A.diagonal()).sum()
        A *= conditioner
        areas *= conditioner
    
    return A, areas

def gen_cotangent_laplacian_matrix( mesh, A ):
    '''
    Places the Laplacian operator matrix for 'mesh' into the upper-left of 'A'.
    Thus, 'A' must have at least len(mesh.vs) rows and columns.
    NOTE: The Laplacian operator matrix contains one row per vertex, not three.
    '''
    
    import laplacian_editing
    
    N = len( mesh.vs )
    assert A is None or (A.shape[0] >= N and A.shape[1] >= N)
    
    areas = zerosf( N )
    
    ijs = []
    vals = []
    
    for vi in range( N ):
        neighbor_indices = mesh.vertex_vertex_neighbors( vi )
        has_boundary = mesh.vertex_is_boundary( vi )
        
        cot_alpha, cot_beta, area = laplacian_editing.cotangentWeights(
            mesh.vs[ vi ],
            [ mesh.vs[ ni ] for ni in neighbor_indices ],
            has_boundary
            )
        
        areas[vi] = area
        
        indices, values = laplacian_editing.Hn_coefficients(
            vi,
            neighbor_indices,
            has_boundary,
            cot_alpha, cot_beta, area
            )
        
        #print values
        
        #A[ vi, tuple(indices) ] = values
        ## csc matrices don't support writable slices
        #for col, val in izip( indices, values ): A[ vi, col ] = val
        
        ijs.extend( [ (vi,j) for j in indices ] )
        vals.extend( values )
    
    assert False not in [ isfinite( v ) for v in vals ]
    
    ijs = asarray(ijs,dtype=int32).T
    #print ijs
    #print vals
    A = solve_matrix_t( scipy.sparse.coo_matrix( (vals, ijs), shape = (N,N) ) )
    
    #print 'L.diagonal():'
    #print A.diagonal()
    #print 'L:'
    #print A
    #debugger()
    
    ## Scale A by the reciprocal of the average value of the diagonal, to that
    ## diagonals are, on average, close to 1.
    ## This has two effects:
    ## 1) The condition number of the system should be improved.
    ## 2) The magnitude of the energy using this matrix should be roughly normalized,
    ##    so that the weight of the least-squares constraint energy is meaningful.
    ## UPDATE: On second thought, scaling this matrix will take it *out* of scale
    ##         from other terms in the energy which are scale-dependent,
    ##         including position displacement.
    if kScaleLaplacianWeights:
        conditioner = float(N) / abs(A.diagonal()).sum()
        A *= conditioner
        areas *= conditioner
    
    #print 'L.diagonal():'
    #print A.diagonal()
    
    ## Actually, this return value is used to compensate for the unweighted laplacian matrix's lack
    ## of area.
    ## Since we take area into account in the coefficients above, return ones.
    return A, areas
    #return A, onesf( N )

def optimize_surface_linear( mesh, constraints, **kwargs ):
    '''
    Given a TriMesh 'mesh', Constraints 'constraints',
    minimizes linearized bending energy to find new mesh vertex positions.
    If 'graph_not_cotangent_weights' is False (the default is false), cotangent weights are used
    instead of graph laplacian weights.
    If 'vertex_normal_motion_only' is True (the default is false), vertices are only allowed to move
    in their normal directions.
    
    Returns 0 (to be compatible with optimize_surface2()).
    
    NOTE: Modifies 'mesh' in place.
    '''
    
    kwargs.setdefault( "graph_not_cotangent_weights", True )
    kwargs.setdefault( "vertex_normal_motion_only", False )
    kwargs.setdefault( "name", "" )
    
    graph_not_cotangent_weights = kwargs["graph_not_cotangent_weights"]
    vertex_normal_motion_only = kwargs["vertex_normal_motion_only"]
    name = kwargs["name"]
    
    from scipy import sparse
    
    ## Mesh vertices should be unique (must they be, or only positions of constraints?)
    assert len(set([ tuple(v) for v in mesh.vs ])) == len( mesh.vs )
    ## We want our mesh vertices to be a 2d array.
    mesh.vs = asarrayf( mesh.vs )
    
    #debugger()
    
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
    L, vert_areas = gen_laplacian( mesh, None )
    L2 = L*L
    
    import mymath
    xyz_hard_linear_constraints = list( constraints.mymath_xy_hard_and_linear_constraints() )
    xyz_hard_linear_constraints.append( constraints.mymath_z_hard_and_linear_constraints() )
    
    print('xyz_hard_linear_constraints:')
    pprint( xyz_hard_linear_constraints )
    #debugger()
    
    print('general_linear_constraints:')
    pprint( constraints.general_linear_constraints )
    
    ## We're going to assume that the system is coupled.
    ## So, put all degrees of freedom into a 3N system.
    
    A3 = assembly_matrix_t( ( 3*N, 3*N ), dtype = L2.dtype )
    A3[:N,:N] = L2
    A3[N:2*N,N:2*N] = L2
    A3[2*N:3*N,2*N:3*N] = L2
    
    b3 = zeros( 3*N, A3.dtype )
    
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
    
    solution3 = mymath.solve_system_with_constraints(
        A3, b3,
        hard_constraints,
        linear_constraints,
        'lsq',
        w_position
        )
    
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

def improve_condition( system, rhs ):
    '''
    Given a system matrix 'system' and right-hand side 'rhs' such that
    system * x = rhs, returns a new (system, rhs) with improved condition number
    that still solve for the same x.
    '''
    
    ## This code is slow.
    return system, rhs
    
    from scipy import sparse
    
    #from numpy.linalg import det
    #print 'before det( system ):', det( system.todense() )
    #if det( system.todense() ) < 1e-10: debugger()
    
    N = system.shape[0]
    assert system.shape[1] == N
    assert rhs.shape[0] == N
    
    ## We will multiply each row of the system so that the diagonal element becomes 1.
    conditioner = 1. / system.diagonal()
    
    ## Modify the matrix.
    ## This is wrong (maybe old scipy?):
    #system = sparse.spdiags( conditioner, 0, N, N ) * system
    ## Why does the following fail?  Some kind of scipy bug aggravated in optimize mode (python -O).
    #diagonal_matrix = solve_matrix_t( ( conditioner, arange(N).reshape((1,N)).repeat(2,0) ), (N, N) )
    diagonal_matrix = sparse.identity(N).tocsr()
    diagonal_matrix.setdiag( conditioner )
    #for i in xrange(N): diagonal_matrix[i,i] = conditioner[i]
    system = diagonal_matrix * system
    
    ## Modify the right-hand-side.
    if len( rhs.shape ) > 1: conditioner = conditioner.reshape( (N,1) )
    rhs = multiply( conditioner, rhs )
    
    #print 'after det( system ):', det( system.todense() )
    
    return system, rhs

def curvature_1d_fibermesh( v0, v1, v2 ):
    '''
    Inspired by the disassembled Java code of FiberMesh.
    '''
    
    v0 = asarrayf(v0)
    v1 = asarrayf(v1)
    v2 = asarrayf(v2)
    
    
    count = 2
    l = mag( v0 - v1 ) + mag( v1 - v2 )
    center = v0 + v2
    center *= 1./count
    l /= count
    vec = v1 - center
    
    ## consistent 90 degree rotation, CCW
    def rot90( xy ):
        return array( (-xy[1], xy[0]) )
    
    normal = dir( dir( rot90( v1 - v0 ) ) + dir( rot90( v2 - v1 ) ) )
    
    d = dot( vec, normal )
    w = 1. * l*l
    return abs( d/w )


def curvature_1d( v0, v1, v2 ):
    '''
    Returns the curvature at vertex v1 in the sequence v0, v1, v2.
    
    untested
    '''
    
    return curvature_1d_fibermesh( v0, v1, v2 )
    return mag( curve_laplacian_1d( v0, v1, v2 ) )
    
    v0 = asarrayf(v0)
    v1 = asarrayf(v1)
    v2 = asarrayf(v2)
    
    
    ## From Nealen's PhD thesis, page 22, which references Grinspun at al. [GDP+06].
    from helpers import angle2D
    return 2*sin( .5 * angle2D( v1 - v0, v2 - v1 ) )
    
def curve_laplacian_1d( v0, v1, v2 ):
    '''
    Returns the laplacian at vertex v1 in the sequence v0, v1, v2.
    
    untested
    '''
    
    v0 = asarrayf(v0)
    v1 = asarrayf(v1)
    v2 = asarrayf(v2)
    
    
    ## From Nealen's PhD thesis, page 23.
    a = mag( v2 - v1 )
    b = mag( v1 - v0 )
    w10 = a / (a+b)
    w12 = b / (a+b)
    
    L1w = v1 - w10 * v0 - w12 * v2
    return L1w

def save_stroke( stroke, filename ):
    out = open( filename, 'w' )
    print('stroke = [', file=out)
    for pt in stroke[:-1]: print('\t( %s, %s ),' % tuple(pt), file=out)
    print('\t( %s, %s )' % tuple( stroke[-1] ), file=out)
    print('\t]', file=out)
    print('', file=out)
    print('mesh = None', file=out)
    print('def main():', file=out)
    print('\tglobal mesh', file=out)
    print('\timport fibermesh', file=out)
    print('\timport sys', file=out)
    print('\tmesh = fibermesh.inflate( *( [stroke] + map( int, sys.argv[1:2] ) + sys.argv[2:] ) )', file=out)
    print('', file=out)
    print('if __name__ == "__main__": main()', file=out)
    del out

def test_display_mesh( mesh ):
    import sys
    try:
        #from OpenGL.GLUT import *
        from OpenGL.GLUT import glutInit, glutInitDisplayMode, GLUT_RGBA, GLUT_DOUBLE, GLUT_DEPTH, glutMainLoop
        #from OpenGL.GL import *
        #from OpenGL.GLU import *
    except:
        print('''
    ERROR: PyOpenGL not installed properly. 
            ''')
        sys.exit()
    
    from mesh_gui import TriMeshWindow
    
    glutInit( sys.argv )
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH )
    
    w1 = TriMeshWindow()
    w1.setWindowTitle( "test1" )
    w1.mesh = mesh
    
    glutMainLoop()

def test1():
    mesh = gen_regular_grid( xin = -1./3, xmax = .5, ymin = -.25, ymax = 1., xverts = 5, yverts = 7 )
    test_display_mesh( mesh )

def test2():
    #curve = [ (0,0), (1,0), (1,1), (0,1) ]
    
    N = 20
    #curve = (
    #    [ ( float(i)/N, 0. ) for i in xrange(N) ] +
    #    [ ( 1., float(i)/N ) for i in xrange(N) ] +
    #    [ ( 1. - float(i)/N, 1. ) for i in xrange(N) ] +
    #    [ ( 0., 1. - float(i)/N ) for i in xrange(N) ]
    #    )
    
    curve = [ ( cos( theta ), sin( theta ) ) for theta in linspace( 0., 2*pi, 4*N )[:-1] ]
    
    
    ## Smooth the curve
    curve = smooth_curve( curve, smooth = True )
    mesh = inflate( curve )
    #test_display_mesh( mesh )
    return mesh

def main():
    #test1()
    test2()

if __name__ == '__main__': main()
