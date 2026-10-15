from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *
from helpers import *
from vecutil import *
from surface_intersections import ray_mesh_intersections


try: from pydb import debugger
except ImportError:
    import pdb
    def debugger(): pdb.set_trace()

kOldZConstraint = False
kSilhouetteConstraints = True
kUseGenSilhouetteConstraint = None

class CurveBoundaryRegion( object ):
    def __init__( self, region_boundary ):
        ## A simple, closed 3d curve representing the boundary of the region.
        self.region_boundary = asarray( region_boundary )
        ## It should be a 3D curve (or empty).
        assert self.region_boundary.shape == ( len( self.region_boundary ), 3 ) or len( region_boundary ) == 0
        
        ## The inflated shape (triangle mesh).
        self.mesh = None
        
        ## Internal structure used during inflation of the mesh.
        self.llconstraints = None
        self.last_iterate_diff = 1e20
        self.iterate_threshold = 1e-4
        
        self.self_test()
    
    def get_mesh( self ): return self._mesh
    def set_mesh( self, mesh ):
        self._mesh = mesh
        self.last_iterate_diff = 1e20
    mesh = property( get_mesh, set_mesh )
    
    def get_llconstraints( self ): return self._llconstraints
    def set_llconstraints( self, llconstraints ):
        self._llconstraints = llconstraints
        self.last_iterate_diff = 1e20
    llconstraints = property( get_llconstraints, set_llconstraints )
    
    def self_test( self ):
        ## Make sure the region boundary has more than 2 vertices.
        ## UPDATE: We can't assert this, because sometimes the object is
        ##         constructed with an empty list.
        #assert len( self.region_boundary ) > 2
        assert asarray( self.region_boundary ).shape == ( len( self.region_boundary ), 3 ) or len( self.region_boundary ) == 0
    
    def iterate( self, N = 1, fast = True, use_cotangent = False, linear = False ):
        '''
        Returns True if something changed, False otherwise.
        (Will always return True if there is a mesh with constraints.)
        '''
        
        assert N >= 0
        if N == 0: return False
        
        if self.mesh and self.llconstraints:
            
            if linear:
                from curvesurface import optimize_surface_linear as optimize_surface
            else:
                #from curvesurface import optimize_surface2 as optimize_surface
                from curvesurface import optimize_surface_linear as optimize_surface
            
            #from fibermesh_cpp import optimize_surface3 as optimize_surface
            #self.last_iterate_diff = optimize_surface( self.mesh, self.llconstraints, N )
            
            print('self.last_iterate_diff:', self.last_iterate_diff)
            print('self.iterate_threshold:', self.iterate_threshold)
            if False: #self.last_iterate_diff < self.iterate_threshold:
                print("Already converged, not iterating")
                return False
            else:
                #self.last_iterate_diff = optimize_surface( self.mesh, self.llconstraints, 40, self.iterate_threshold, True )
                #self.last_iterate_diff = optimize_surface( self.mesh, self.llconstraints, N, self.iterate_threshold, False, True )
                #self.last_iterate_diff = optimize_surface( self.mesh, self.llconstraints, N, self.iterate_threshold, True, True )
                #self.last_iterate_diff = optimize_surface( self.mesh, self.llconstraints, N, self.iterate_threshold, True, False )
                
                self.last_iterate_diff = optimize_surface( self.mesh, self.llconstraints, max_iterations = N, absolute_convergence_threshold = self.iterate_threshold, graph_not_cotangent_weights = not use_cotangent, vertex_normal_motion_only = use_cotangent, stretching_weight = 1. )
            
            return self.last_iterate_diff >= self.iterate_threshold
        else:
            print("No mesh")
            return False
    
    def update_llconstraints( self ):
        print("WARNING: Who is calling this?")
        
        mesh, llconstraints = mesh_llconstraints_from_closed_curve_3D(
            self.region_boundary
            )
        self.llconstraints = llconstraints
    
    def ray_surface_intersection_ts( self, world_ray, ray_is_line = False ):
        '''
        Given a world-space ray (3d point, 3d vector direction), returns the t's along the ray
        corresponding to intersections with the shape's 3d surface.
        If 'ray_is_line' is False, only intersections in the positive direction along the ray
        will be returned.
        
        untested
        '''
        
        if self.mesh is None: raise RuntimeError("self.mesh is None!")
        
        assert len( world_ray[0] ) == 3
        assert len( world_ray[1] ) == 3
        world_ray = ( asarrayf( world_ray[0] ).ravel(), asarrayf( world_ray[1] ).ravel() )
        
        from surface_intersections import ray_mesh_intersections
        isects = ray_mesh_intersections( world_ray, self.mesh )
        ts = [ isect[0] for isect in isects ]
        if not ray_is_line: ts = [ t for t in ts if t > -kEps ]
        return ts
    
    def ray_surface_intersections( self, world_ray, ray_is_line = False ):
        '''
        Given a world-space ray (3d point, 3d vector direction), returns the 3d intersection points
        with the shape's 3d surface.
        If 'ray_is_line' is False, only intersections in the positive direction along the ray
        will be returned.
        
        untested
        '''
        
        ts = self.ray_surface_intersection_ts( world_ray, ray_is_line )
        return [ asarray( world_ray[0] ) + t * asarray( world_ray[1] ) for t in ts ]
    
    def offset_in_depth( self, z_offset, oblique ):
        '''
        Translate the object along the look vector (oblique[0], oblique[1], -1) by 'z_offset'.
        Inverts the oblique projection, so the object maintains the same (projected) appearance.
        '''
        
        self.mesh.vs = asarray( self.mesh.vs )
        self.mesh.vs[:,:] += z_offset * array( ( oblique[0], oblique[1], -1 ) )
        self.mesh.positions_changed()

def curve_3D_mesh_llconstraints_for_N_by_M_grid( N, M ):
    '''
    Returns a 3D curve, a mesh, and low-level constraints
    for an N-by-N grid.
    '''
    
    print("[[[Entered curve_3D_mesh_llconstraints_for_N_by_M_grid()]]]")
    
    name = ''
    
    assert int( N ) == N
    N = int( N )
    
    assert int( M ) == M
    M = int( M )
    
    from fibermesh import gen_regular_grid
    mesh = gen_regular_grid( xmin = -.5, xmax = .5, ymin = -.5, ymax = .5, xverts = N, yverts = M )
    mesh.write_OBJ( name + 'glued.obj' )
    mesh.vs = asarray( mesh.vs )
    
    xy_verts = []
    xy_pos = []
    
    z_verts = []
    z_vals = []
    
    ## 2
    ## mesh.boundary_vertices() is not necessarily in the same order as 'closed_curve_ND'.
    boundary_vertices = mesh.boundary_vertices()
    xy_verts = list( boundary_vertices )
    xy_pos = [ mesh.vs[ vi, :2 ] for vi in boundary_vertices ]
    z_verts = list( boundary_vertices )
    z_vals = [ mesh.vs[ vi, 2 ] for vi in boundary_vertices ]
    
    from fibermesh import Constraints
    constraints = Constraints()
    constraints.xy_constraint_vertices = xy_verts
    constraints.xy_constraint_positions = xy_pos
    constraints.z_constraint_vertices = z_verts
    constraints.z_constraint_values = z_vals
    
    ## These are the points on the curve, but they're not in sequential order.
    curve_3D = [ mesh.vs[ vi ] for vi in boundary_vertices ]
    ## Sort according to atan2 of the xy.
    from math import atan2
    curve_3D.sort( key = lambda xyz: atan2( xyz[1], xyz[0] ) )
    curve_3D = asarray( curve_3D )
    
    return ( curve_3D, mesh, constraints )

def mesh_llconstraints_from_closed_curve_3D( closed_curve_3D ):
    '''
    Given a sequence of 3D points representing a closed curve in space,
    returns the tuple ( mesh, low-level constraints ).
    These elements are suitable for passing to fibermesh.optimize_surface2().
    
    untested
    '''
    
    print("[[[Entered mesh_llconstraints_from_closed_curve_3D()]]]")
    
    name = ''
    
    #debugger()
    
    ## In this function, we are following steps 1-5 of inflate() in
    ## "annotations revisited data structures.txt":
    ## 1 create a disk (2d triangle-meshe) inside 'closed_curve_3D'.
    ## 2 constrain 'closed_curve_3D' in 3D.
    
    from fibermesh import gen_glued_mesh, attach_curve2d_to_mesh_vertices, Constraints
    
    ## 1
    closed_curve_2D = map_closed_curve_3D_to_2D( closed_curve_3D )
    #mesh = gen_mesh_from_closed_curve_2D_irregular( closed_curve_2D )
    mesh = gen_mesh_from_closed_curve_2D_grid( closed_curve_2D, 5 )
    mesh.write_OBJ( name + 'glued.obj' )
    mesh.vs = asarray( mesh.vs )
    
    xy_verts = []
    xy_pos = []
    
    z_verts = []
    z_vals = []
    
    ## 2
    ## mesh.boundary_vertices() is not necessarily in the same order as 'closed_curve_ND'.
    #boundary_vertices = mesh.boundary_vertices()
    xy_verts = list(range(len( closed_curve_3D )))
    xy_pos = list( asarray( closed_curve_3D )[:,:2] )
    z_verts = list( xy_verts )
    z_vals = list( asarray( closed_curve_3D )[:,2] )
    
    constraints = Constraints()
    constraints.xy_constraint_vertices = xy_verts
    constraints.xy_constraint_positions = xy_pos
    constraints.z_constraint_vertices = z_verts
    constraints.z_constraint_values = z_vals
    
    return ( mesh, constraints )

def gen_mesh_from_closed_curve_2D_sloppy( closed_curve_2D ):
    '''
    Given a sequence of 2D points 'closed_curve_2D', returns three values:
        a TriMesh object approximately triangulating it,
        a list 'boundary_vertices' of indices into the mesh's '.vs' property
        representing boundary vertices,
        and a sequence 'boundary_vi_to_curve_vi' mapping
        indices in 'boundary_vertices' to indices in 'closed_curve_2D'.
    
    NOTE: The returned mesh does not have 'closed_curve_2D' as its boundary,
          and may even have additional points.  Use the map 'boundary_vi_to_curve_vi'
          to determine which 'mesh' vertices indicated by 'boundary_vertices'
          should be constrained to which locations in 'closed_curve_2D'
          by the optimization.
    '''
    
    from fibermesh import gen_glued_mesh
    
    mesh, boundary_vertices, front_vertices, back_vertices = gen_glued_mesh( closed_curve_2D )
    mesh.remove_vertex_indices( back_vertices )
    assert set( mesh.boundary_vertices() ) == set( boundary_vertices )
    
    mesh.vs = asarray( mesh.vs )
    mesh.vs[:,2] = 0.
    mesh.positions_changed()
    
    boundary_vertices, boundary_vi_to_curve_vi = attach_curve2d_to_mesh_vertices( closed_curve_2D, mesh, mesh.boundary_vertices() )
    
    return mesh, boundary_vertices, boundary_vi_to_curve_vi

def gen_mesh_from_closed_curve_2D_and_interior_points( closed_curve_2D, interior_points ):
    '''
    Given a sequence of 2D points 'closed_curve_2D'
    representing a boundary curve and
    a set of 2D points 'interior_points' on the interior of the curve,
    returns a triangulated TriMesh with 'closed_curve_2D'
    as the first N vertices of the TriMesh's '.vs' field.
    '''
    
    ## Call Triangle.
    vs = list( closed_curve_2D ) + list( interior_points )
    from triangle import triangles_for_points
    triangles = triangles_for_points(
        vs,
        [ ( i, i+1 ) for i in range( len( closed_curve_2D ) - 1 ) ] + [ ( len( closed_curve_2D ) - 1, 0 ) ]
        )
    
    from mesh_gui import TriMesh
    mesh = TriMesh()
    mesh.vs = zeros( ( len( vs ), 3 ), dtype = float )
    mesh.vs[:,:2] = vs
    mesh.faces = triangles
    mesh.topology_changed()
    return mesh

def gen_mesh_from_closed_curve_2D_irregular( closed_curve_2D, num_interior_points = None ):
    '''
    Given a sequence of 2D points 'closed_curve_2D', returns a TriMesh
    homeomorphic to a disc whose boundary vertices are 'closed_curve_2D'
    and whose '.vs' field begins with 'closed_curve_2D'.
    NOTE: The vertices returned by TriMesh.boundary_vertices() are not
          necessarily in the same order as 'closed_curve_2D'.
    
    Optional parameter 'num_interior_points' specifies the number
    of interior points in the resulting mesh.
    '''
    
    if num_interior_points is None:
        #num_interior_points = len( closed_curve_2D ) * len( closed_curve_2D )
        num_interior_points = len( closed_curve_2D ) * 10
    
    interior_points = generate_N_points_inside_the_closed_curve_2D_random( closed_curve_2D, num_interior_points )
    
    return gen_mesh_from_closed_curve_2D_and_interior_points( closed_curve_2D, interior_points )

def gen_mesh_from_closed_curve_2D_grid( closed_curve_2D, num_interior_points = None ):
    '''
    Given a sequence of 2D points 'closed_curve_2D', returns a TriMesh
    homeomorphic to a disc whose boundary vertices are 'closed_curve_2D'
    and whose '.vs' field begins with 'closed_curve_2D'.
    NOTE: The vertices returned by TriMesh.boundary_vertices() are not
          necessarily in the same order as 'closed_curve_2D'.
    
    Optional parameter 'num_interior_points' specifies the minimum number
    of interior points in the resulting mesh.
    '''
    
    if num_interior_points is None:
        #num_interior_points = len( closed_curve_2D ) * len( closed_curve_2D )
        num_interior_points = len( closed_curve_2D ) * 10
    
    interior_points = generate_N_points_inside_the_closed_curve_2D_grid( closed_curve_2D, num_interior_points )
    
    return gen_mesh_from_closed_curve_2D_and_interior_points( closed_curve_2D, interior_points )

def generate_N_points_inside_the_closed_curve_2D_grid( closed_curve_2D, num_interior_points ):
    '''
    Given a sequence of 2D points 'closed_curve_2D',
    returns 'num_interior_points' points inside.
    
    NOTE: Returns at least N points.
    '''
    
    closed_curve_2D = asarray( closed_curve_2D )
    edges = append( closed_curve_2D, [closed_curve_2D[0]], axis = 0 )
    edges = concatenate( ( edges[:-1][:,newaxis,:], edges[1:][:,newaxis,:] ), axis = 1 )
    
    ### 1 Generate a grid and reject points outside the contour.
    ### 2 Count the remaining points, search for the parameter
    ###   that returns approximately (or at least) the right number of points.
    
    
    ### 1
    
    from helpers import get_contours_bbox
    from pnpoly_c import line_loop_contains_point as contour_contains
    ## bbox = ( minx, maxx, miny, maxy )
    bbox = get_contours_bbox( [closed_curve_2D] )
    
    ## hchelpers is ~/Work/TAU/hasm-code/experiments/src/helpers.py
    from hchelpers import points_min_and_max_distance_from_each_other
    min_dist = .25*points_min_and_max_distance_from_each_other( closed_curve_2D )[0]
    
    N = max( 1, ceil( sqrt( num_interior_points ) ) )
    while True:
        print('generate_N_points_inside_the_closed_curve_2D_grid(): Trying grid size N-by-N with N:', N)
        
        ## TODO Q: If I scale a grid to fit inside the bounding box, do I want
        ##         to take care to keep delta-x and delta-y uniform?
        ## A1: For now, keep it simple and make the number of
        ##     rows and columns the same.
        
        Xs, Ys = meshgrid( linspace( bbox[0], bbox[1], N+2 )[1:-1], linspace( bbox[2], bbox[3], N+2 )[1:-1] )
        grid_points = list( zip( Xs.ravel(), Ys.ravel() ) )
        if len( grid_points ) > num_interior_points:
            
            ## 2
            interior_points = [ pt for pt in grid_points if contour_contains( closed_curve_2D, pt ) ]
            if len( interior_points ) > num_interior_points:
                
                print('passed interior point test; are they far from the line loop?')
                #interior_points = [ pt for pt in interior_points if closest_distsqr_and_point_and_edge_index_on_line_loop_to_point( closed_curve_2D, pt )[0] > min_dist*min_dist ]
                
                from edge_distances import distancesSqr_to_edges
                interior_points = asarray( interior_points )
                min_distsqr_to_edges = distancesSqr_to_edges( interior_points.T, edges ).min( -2 )
                interior_points = interior_points[ min_distsqr_to_edges > min_dist*min_dist ]
                
                if len( interior_points ) > num_interior_points:
                    break
        
        ## It's N-squared, so increase by sqrt(1.5)
        N = int(round( sqrt(1.5) * N ))
    
    return interior_points

def generate_N_points_inside_the_closed_curve_2D_random( closed_curve_2D, num_interior_points ):
    '''
    Given a sequence of 2D points 'closed_curve_2D',
    returns 'num_interior_points' points inside.
    
    fails to return in a reasonable amount of time
    '''
    
    closed_curve_2D = asarray( closed_curve_2D )
    
    ## We expect a 2D curve, so the shape should be N-by-2:
    assert closed_curve_2D.shape == ( len( closed_curve_2D ), 2 )
    ## We require more than 2 points in the curve.
    assert closed_curve_2D.shape[0] > 2
    
    ## Generate at least 1 interior point.
    assert num_interior_points > 0
    
    from helpers import get_contours_bbox
    from pnpoly_c import line_loop_contains_point as contour_contains
    ## bbox = ( minx, maxx, miny, maxy )
    bbox = get_contours_bbox( [closed_curve_2D] )
    def sample_point_inside_bbox():
        x, y = random.uniform( bbox[0], bbox[1] ), random.uniform( bbox[2], bbox[3] )
        return x, y
    def sample_point_inside_closed_curve_2D():
        while True:
            x, y = sample_point_inside_bbox()
            if contour_contains( closed_curve_2D, ( x,y ) ): break
        
        return ( x, y )
    
    ## hchelpers is ~/Work/TAU/hasm-code/experiments/src/helpers.py
    from hchelpers import points_min_and_max_distance_from_each_other, sample_N_points_more_than_d_from_each_other, all_points_are_more_than_d_from_p
    
    min_dist = .25*points_min_and_max_distance_from_each_other( closed_curve_2D )[0]
    
    def sample_point_inside_and_far_from_closed_curve_2D():
        while True:
            x, y = sample_point_inside_closed_curve_2D()
            if all_points_are_more_than_d_from_p( closed_curve_2D, min_dist, ( x, y ) ): break
        
        return ( x, y )
    
    interior_points = sample_N_points_more_than_d_from_each_other(
        num_interior_points,
        min_dist,
        sample_point_inside_and_far_from_closed_curve_2D
        )
    
    return interior_points

def map_closed_curve_3D_to_2D( closed_curve_3D ):
    closed_curve_3D = asarray( closed_curve_3D )
    
    ## We expect a 3D curve, so the shape should be N-by-3:
    assert closed_curve_3D.shape == ( len( closed_curve_3D ), 3 )
    ## We require more than 2 points in the curve.
    N = closed_curve_3D.shape[0]
    assert N > 2
    
    ts = linspace( 0, 2*pi, N+1 )[:-1]
    closed_curve_2D = empty( ( closed_curve_3D.shape[0], 2 ), dtype = float )
    closed_curve_2D[:,0] = cos( ts )
    closed_curve_2D[:,1] = sin( ts )
    return closed_curve_2D

def test_contour_interior():
    pts = [ ( -1,-1 ), ( 1, -1 ), ( 1, 1 ), ( -1, 1 ) ]
    
    ## This isn't very good, because 4 random points may be self-intersecting
    ## when viewed as a polyline loop.
    #import random
    #pts = [ ( random.uniform( -1, 1 ), random.uniform( -1, 1 ) ) for i in xrange(4) ]
    
    print('pts:', pts)
    #interior_pts = generate_N_points_inside_the_closed_curve_2D_random( pts, 10 )
    #print 'generate_N_points_inside_the_closed_curve_2D_random( pts, 10 ):'
    interior_pts = generate_N_points_inside_the_closed_curve_2D_grid( pts, 10 )
    #interior_pts = generate_N_points_inside_the_closed_curve_2D_grid( pts, 25 )
    print('generate_N_points_inside_the_closed_curve_2D_grid( pts, 10 ):')
    print(interior_pts)
    
    from mesh_gui import TriMeshWindow
    class DotWindow( TriMeshWindow ):
        def displayFunc( self ):
            TriMeshWindow.displayFunc( self )
            
            from OpenGL.GL import glPointSize, glDisable, GL_LIGHTING, glColor3f, glColor4f, glBegin, GL_POINTS, glEnd, glVertex3f, glVertex2f
            glDisable( GL_LIGHTING )
            glPointSize( 10 )
            glBegin( GL_POINTS )
            for dots, color in self.dots_and_color:
                glColor3f( *color ) if len( color ) == 3 else glColor4f( *color )
                for pt in dots:
                    glVertex3f( *pt ) if len( pt ) == 3 else glVertex2f( *pt )
            glEnd()
    
    from OpenGL.GLUT import glutInit, glutInitDisplayMode, GLUT_RGBA, GLUT_DOUBLE, GLUT_DEPTH, glutMainLoop
    import sys
    glutInit( sys.argv )
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH )
    w = DotWindow()
    w.dots_and_color = [ ( pts, ( 1, 0, 0 ) ), ( interior_pts, ( 0, 1, 0 ) ) ]
    glutMainLoop()

def main():
    test_contour_interior()

if __name__ == '__main__': main()
