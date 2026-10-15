from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *
from helpers import *
from vecutil import *
from geometry import line_line_intersection
from surface_intersections import ray_mesh_intersections


try: from pydb import debugger
except ImportError:
    import pdb
    def debugger(): pdb.set_trace()

kNewCrossSectionConstraints = True
kCrossSectionConstraints1ConstrainsEdgeNormals = True
## This is a mod value, so 1 means no subsetting, and larger numbers are the number to skip - 1.
kCrossSectionConstraints1ConstrainsEdgeNormalsSubset = 2
kCrossSectionConstraints2ConstrainsIntersectionPoints = False ## These are broken!  Why do I even want them?
kCrossSectionConstraints2ConstrainsEdgeRotation = False
kCrossSectionConstraints2ConstrainsEdgeRotation2 = True
kCrossSectionConstraints2ConstrainsCullSome = False
kOldZConstraint = False
kSilhouetteConstraints = True
kUseGenSilhouetteConstraint = None

class RotundRegion( object ):
    def __init__( self, region_boundary, silhouette_edges, sharp_cut_edges, spine = None ):
        ## A simple, closed 2d curve representing the boundary of the region.
        self.region_boundary = region_boundary
        
        ## The index of tagged edges in region_boundary.
        ## Each index i corresponds to ( region_boundary[i], region_boundary[i+1] ),
        ## where region_boundary[ len(region_boundary) ] loops back to region_boundary[0].
        self.silhouette_edges = silhouette_edges
        self.sharp_cut_edges = sharp_cut_edges
        
        ## A list of marked-up cross sections.
        ## Each entry is: (
        ##                  a 2d anchor point (this is the point picking the cross section on the untilted shape),
        ##                  the destination 2d point at the closest point on the tilted cross section
        ##                )
        ## If there is a third element, then it is a a user-drawn (non-circular) cross section (a simple, closed, CCW 2d curve).
        ## (This can be parameterized by ordering the points in ths list so that the beginning is (1,0) of the unit circle.)
        ## NOTE: There are several ways of mapping anchor points onto cross sections:
        ##       1. Using the spine.
        ##       2. Using a local spine (formed by the vector between the anchor source and destination).
        ##       3. Finding the shortest line segment between two boundary points through the anchor point.
        ##       4. Finding the geodesic shortest path through the point, then straightening it (or not).
        self.cross_sections = []
        
        ## The spine path (a simple, open 2d curve inside of the region).
        self.spine = None
        
        ## The inflated shape (triangle mesh).
        self.mesh = None
        ## The last inflated shape before cross sections were added.
        self.mesh_no_cross_sections = None
        ## The lifted spine (list of 3d vertices).
        self.spine3d = None
        ## A list of line loops corresponding to 3d cross sections constraints.
        self.cross_sections3d = None
        
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
        ## Make sure edge indices are in range.
        assert False not in [ ei >= 0 and ei < len( self.region_boundary ) for ei in self.silhouette_edges ]
        assert False not in [ ei >= 0 and ei < len( self.region_boundary ) for ei in self.sharp_cut_edges ]
        ## Make sure there are no duplicate edges.
        assert len( set( self.silhouette_edges ) ) == len( self.silhouette_edges )
        assert len( set( self.sharp_cut_edges ) ) == len( self.sharp_cut_edges )
        ## Make sure silhouette_edges and sharp_cut_edges are disjoint.
        assert len( set( self.sharp_cut_edges ).intersection( set( self.silhouette_edges ) ) ) == 0
        ## Make sure all edges in region_boundary are marked.
        ## TODO Q: Do I need this?  Should I remove this restriction?
        ## A1: Yes.  I added support for it.
        # assert False not in [ ei in self.silhouette_edges or ei in self.sharp_cut_edges for ei in xrange(len( self.region_boundary )) ]
    
    def iterate( self, N = 1, fast = True, use_cotangent = False, linear = False ):
        '''
        Returns True if something changed, False otherwise.
        (Will always return True if there is a mesh with constraints.)
        '''
        
        assert N >= 0
        if N == 0: return False
        
        if self.mesh and self.llconstraints:
            
            if linear:
                if fast:
                    from fibermesh_cpp import optimize_surface_linear as optimize_surface
                else:
                    from fibermesh import optimize_surface_linear as optimize_surface
            else:
                if fast:
                    from fibermesh_cpp import optimize_surface2 as optimize_surface
                else:
                    from fibermesh import optimize_surface2 as optimize_surface
            
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
                
                self.last_iterate_diff = optimize_surface( self.mesh, self.llconstraints, max_iterations = N, absolute_convergence_threshold = self.iterate_threshold, graph_not_cotangent_weights = not use_cotangent, vertex_normal_motion_only = use_cotangent )
                if self.mesh_no_cross_sections is None:
                    from mesh_gui import TriMesh
                    self.mesh_no_cross_sections = TriMesh()
                    self.mesh_no_cross_sections.faces = array( self.mesh.faces )
                    self.mesh_no_cross_sections.vs = array( self.mesh.vs )
                    self.mesh_no_cross_sections.positions_changed()
                elif len( self.cross_sections ) == 0:
                    self.mesh_no_cross_sections.vs[:] = self.mesh.vs
                    self.mesh_no_cross_sections.positions_changed()
            
            return self.last_iterate_diff >= self.iterate_threshold
        else:
            print("No mesh")
            return False
    
    def update_llconstraints( self, oblique = None ):
        if oblique is None:
            mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2(
                self
                )
        else:
            mesh, llconstraints, xyz_cross_sections = mesh_llconstraints_from_constraints2_oblique(
                self, oblique
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


class Connection( object ):
    def __init__( self, connection_curve, region1, region1_attaches_front, region2, region2_attaches_front ):
        ## A user-drawn (simple, open, 2d) curve specifying the connection (intersection curve) between two regions.
        self.connection_curve = connection_curve
        ## An internal, generated 3d intersection curve between the two regions' inflated mesh representations.
        self.connection_curve_precise = None
        
        ## A pair of references to RotundRegions.  These are the regions connected by this curve.
        self.regions = ( region1, region2 )
        ## A pair of booleans specifying which sides of the inflated regions are attached along
        ## connection_curve.  True means the front, False means the back.
        ## NOTE: In the case of a region having more than two faces (front, back, front, back, ...)
        ##       under the connection curve, we can store an integer offset into the sorted list of
        ##       intersection t's along the ray under the point in the look direction.
        self.regions_attaches_front = ( region1_attaches_front, region2_attaches_front )
    
    def ReplaceConnectionCurve( self, curve ):
        self.connection_curve = curve
        self.connection_curve_precise = None


def mesh_llconstraints_from_constraints( silhouette_curve, cross_sections ):
    '''
    Given a 2d curve (list of 2d points) for the silhouette,
    and a list of 2d curves for cross sections,
    returns a ( mesh, low-level constraints, list of cross sections in 3d ).
    The first two elements are suitable for passing to fibermesh.optimize_surface2(),
    and the third element is a list of 3d cross sections curves corresponding to each 2d cross
    section in 'cross_sections'.
    
    This assumes that the curves have already been smoothed and self-intersections removed.
    '''
    
    name = ''
    
    from fibermesh import gen_glued_mesh, attach_curve2d_to_mesh_vertices
    
    mesh, boundary_vertices, front_vertices, back_vertices = gen_glued_mesh( silhouette_curve )
    mesh.write_OBJ( name + 'glued.obj' )
    
    ## Get the silhouette (xy) constraint
    boundary_vertices, boundary_vi_to_curve_vi = attach_curve2d_to_mesh_vertices( silhouette_curve, mesh, boundary_vertices )
    xy_verts = list( boundary_vertices )
    xy_pos = [ silhouette_curve[ si ] for si in boundary_vi_to_curve_vi ]
    
    
    ## Get the relative z constraints
    ## 1 Find straight lines (from the beginning of the cross section to the end).
    ## 2 Intersect the straight lines with the mesh.
    ## 3 Assign (relative) z values to the mesh vertices near the intersection.
    ## 4 Also add new xy constraints for the nearby mesh vertices.
    
    
    ## 1
    straight_cross_sections = [ straighten( cross_section ) for cross_section in cross_sections ]
    
    
    ## 2
    xsections_intersected_edges = [ intersect_curve_with_mesh_edges( straight_cross_section, mesh ) for straight_cross_section in straight_cross_sections ]
    xsection_nearby_front_vertices = []
    for cross_section_edges in xsections_intersected_edges:
        xsection_nearby_front_vertices.append( [] )
        for i,j in cross_section_edges:
            ## i and j must be in front_vertices and must not have already been used by an older
            ## cross section.
            if i in front_vertices and True not in [ i in vs for vs in xsection_nearby_front_vertices[:-1] ]:
                xsection_nearby_front_vertices[-1].append( i )
            if j in front_vertices and True not in [ j in vs for vs in xsection_nearby_front_vertices[:-1] ]:
                xsection_nearby_front_vertices[-1].append( j )
        xsection_nearby_front_vertices[-1] = list(set( xsection_nearby_front_vertices[-1] ))
    assert len( cross_sections ) == len( straight_cross_sections )
    assert len( straight_cross_sections ) == len( xsection_nearby_front_vertices )
    cross_section_verts_curv_indices = [
        attach_curve2d_to_mesh_vertices( straight_cross_section, mesh, nearby_front_vertices )
        for straight_cross_section, nearby_front_vertices in zip( straight_cross_sections, xsection_nearby_front_vertices )
        ]
    #debugger()
    
    
    ## 3
    zrel_groups_vert = []
    zrel_groups_vals = []
    for ((verts, curve_indices), cross_section) in zip( cross_section_verts_curv_indices, cross_sections ):
        assert len( verts ) == len( curve_indices )
        zrel_groups_vert.append( [] )
        zrel_groups_vals.append( [] )
        for mvi, cvi in zip( verts, curve_indices ):
            zval = magic_zval( cross_section[cvi], cross_section )
            if zval > 0.:
                zrel_groups_vert[-1].append( mvi )
                zrel_groups_vals[-1].append( zval )
        
        #print '## Keep only two.  Why is this still so harsh?'
        #del zrel_groups_vert[-1][2:]
        #del zrel_groups_vals[-1][2:]
    
    ## This is the same as above, but doesn't cull if the zval < 0.
    #zrel_groups_vert = [ verts for (verts, curve_indices) in cross_section_verts_curv_indices ]
    #zrel_groups_vals = [
    #    [ magic_zval( cross_section[i], cross_section ) for i in curve_indices ]
    #    for ((verts, curve_indices), cross_section) in izip( cross_section_verts_curv_indices, cross_sections )
    #    ]
    
    
    ## 4
    for (verts, curve_indices), straight_cross_section in zip( cross_section_verts_curv_indices, straight_cross_sections ):
        for mvi, cvi in zip( verts, curve_indices ):
            if mvi not in xy_verts:
                xy_verts.append( mvi )
                #xy_pos.append( straight_cross_section[ cvi ] )
                xy_pos.append( mesh.vs[ mvi ][:2] )
    
    ## We need at least one (arbitrary) z value constraint
    z_verts = [boundary_vertices[0]]
    z_vals = [0.]
    ## Prefer an endpoint in zrel_groups if we have one
    zrel_z_vert = [ verts[0] for verts in zrel_groups_vert if len( verts ) > 0 ][:1]
    #print 'zrel_z_vert:', zrel_z_vert
    if len( zrel_z_vert ) > 0: z_verts = [ zrel_z_vert[0] ]
    
    ## Actually, keep the silhouette orthogonal
    #z_verts = [ mvi for mvi in boundary_vertices if True not in [ mvi in verts for verts in zrel_groups_vert ] ]
    #z_vals = [0.] * len(z_verts)
    #print 'boundary_vertices:', boundary_vertices
    #print 'zverts:', zverts
    
    ### Convert zrel constraints into z constraints
    #for verts, relvals in izip( zrel_groups_vert, zrel_groups_vals ):
    #    for mvi, zrel in izip( verts, relvals ):
    #        if mvi not in z_verts:
    #            z_verts.append( mvi )
    #            z_vals.append( zrel )
    #zrel_groups_vert = []
    #zrel_groups_vals = []
    
    constraints = Struct()
    constraints.xy_constraint_vertices = xy_verts
    constraints.xy_constraint_positions = xy_pos
    constraints.z_constraint_vertices = z_verts
    constraints.z_constraint_values = z_vals
    constraints.z_relative_constrained_groups_vertices = zrel_groups_vert
    constraints.z_relative_constrained_groups_values = zrel_groups_vals
    
    ## This should produce the same results as regular, orthogonal FiberMesh
    #constraints.z_constraint_vertices = xy_verts
    #constraints.z_constraint_values = [0.] * len(xy_verts)
    #constraints.z_relative_constrained_groups_vertices = []
    #constraints.z_relative_constrained_groups_values = []
    
    return (
        mesh, constraints, [
            [ (sx,sy, magic_zval( (cx,cy), cross_section )) for (cx,cy), (sx,sy) in zip( cross_section, straight_cross_section ) ]
            for cross_section, straight_cross_section in zip( cross_sections, straight_cross_sections )
            ]
        )

def magic_zval( curve_space_pt, curve ):
    '''
    Given a 2d point in curve space, return its height from the line between the curve's endpoints.
    The parameter 'curve' is a list of 2d points representing an open path (like the
    front-facing part of a cross section).
    '''
    
    ## 1 Translate, rotate, and (not) uniformly scale 'curve' so its first and last points
    ##   are at (0,0) and (mag(curve[0] - curve[-1],0)
    ## 2 Transform 'curve_space_pt' the same way.
    ## 3 The result is a magic factor times the y value of the transformed 'curve_space_pt'.
    
    xend = mag(curve[0] - curve[-1] )
    
    xform_pts = solve_for_fixed_points( [curve_space_pt], (curve[0], curve[-1]), ((0,0), (xend,0)) )
    height = xform_pts[0][1]
    
    ## Draw a picture of a sphere.  Measure the maximum displacement of the cross section.
    ## That displacement should be the same as its radius.
    ## This is our factor.
    ## I drew a sphere with diameter: 625 - 372 = 253 ~= 250
    ## The cross section I drew had maximum y displacement of 384 - 334 = 50
    ## So, the factor is 5.
    #magic_constant = 5.
    ## 5 seems a high.
    magic_constant = 2.
    ## For debugging, use 1.
    #magic_constant = 1.
    
    return magic_constant * height

def mesh_llconstraints_from_constraints2( region ):
    '''
    Given a RotundRegion 'region',
    returns the tuple ( mesh, low-level constraints, list of 3d curves representing cross sections ).
    The first two elements are suitable for passing to fibermesh.optimize_surface2(),
    and the third element is a list of 3d cross sections curves corresponding to each 2d cross
    section in 'region.cross_sections'.
    
    This assumes that the curves have already been smoothed and self-intersections removed.
    
    tested
    '''
    
    name = ''
    
    ## In this function, we are following steps 1-5 of inflate() in
    ## "annotations revisited data structures.txt":
    ## 1 glue two disks (2d triangle-meshes) together along 'region_boundary',
    ##   slightly offset from each other in z.
    ## 2 constrain 'region_boundary' in xy.
    ## 3 constrain normals of edges in 'silhouette_edges' to have z = 0.
    ## 4 for edges in 'sharp_cut_edges', duplicate the boundary vertices and disconnect the back and
    ##   front vertices from each other.  leave the front vertices constrained in xy, but
    ##   do not constrain the xy or z of the back vertices.
    ## 
    ## 5 for each cross_section in 'cross_sections', let 'cs' be the mesh vertices (front and back)
    ##   near cross_section.  (See the NOTE near the definition of the cross_sections in RotundRegion.)
    ##   Determine the desired tilt by rotating the cross section circle about the line between the
    ##   cross section's two points on the mesh's silhouette such that the anchor point
    ##   (cross_section[0]) is as close as possible to its target point (cross_section[1]).
    ##   On each vertex in 'cs', place an xy constraint at the corresponding image-plane location
    ##   on the tilted cross-section.  Also place a relative z constraint with the first vertex
    ##   in 'cs'.  (This is an arbitrary choice, but the relative constraints must be relative to
    ##   some vertex in the cross section's z).
    
    from fibermesh import gen_glued_mesh, attach_curve2d_to_mesh_vertices, Constraints
    
    ## 1
    mesh, boundary_vertices, front_vertices, back_vertices = gen_glued_mesh( region.region_boundary )
    ## These sets should be disjoint.
    assert len( set( boundary_vertices ).intersection( set( front_vertices ) ) ) == 0
    assert len( set( boundary_vertices ).intersection( set( back_vertices ) ) ) == 0
    assert len( set( front_vertices ).intersection( set( back_vertices ) ) ) == 0
    ## If we already have a mesh, re-use it:
    if region.mesh_no_cross_sections: mesh = region.mesh_no_cross_sections
    mesh.write_OBJ( name + 'glued.obj' )
    mesh.vs = asarray( mesh.vs )
    
    xy_verts = []
    xy_pos = []
    
    zrel_groups_vert = []
    zrel_groups_vals = []
    
    xy_vertex_weighted_sum_indices = []
    xy_vertex_weighted_sum_weights = []
    xy_vertex_weighted_sum_rhss = []
    z_vertex_weighted_sum_indices = []
    z_vertex_weighted_sum_weights = []
    z_vertex_weighted_sum_rhss = []
    
    general_linear_constraints = []
    
    ## 2
    boundary_vertices, boundary_vi_to_curve_vi = attach_curve2d_to_mesh_vertices( region.region_boundary, mesh, boundary_vertices )
    #xy_verts = list( boundary_vertices )
    #xy_pos = [ region.region_boundary[ si ] for si in boundary_vi_to_curve_vi ]
    ## Only constrain boundary vertices that are endpoints of silhouette or sharp-cut edges.
    constrain_vertices = set()
    constrain_vertices.update( region.silhouette_edges )
    constrain_vertices.update( [ ( ei + 1 ) % len( region.region_boundary ) for ei in region.silhouette_edges ] )
    constrain_vertices.update( region.sharp_cut_edges )
    constrain_vertices.update( [ ( ei + 1 ) % len( region.region_boundary ) for ei in region.sharp_cut_edges ] )
    for bvi, cvi in zip( boundary_vertices, boundary_vi_to_curve_vi ):
        if cvi in constrain_vertices:
            xy_verts.append( bvi )
            xy_pos.append( region.region_boundary[ cvi ] )
            
            ## This is what the _oblique version of this function constrains; it is a least-squares
            ## constraint instead of a hard constraint, and this showed that the constraint
            ## weight wasn't large enough.
            '''
            for coord in (0,1):
                general_linear_constraints.append( (
                    ## The list of ( vertex index, coord ) tuples.
                    [ (bvi,coord) ],
                    ## The list of coefficients for each term.
                    ## We add the oblique parameter * z because we want xy.screen = xy.world + ab*z.
                    [ 1. ],
                    ## The right hand side value.
                    region.region_boundary[ cvi ][coord]
                    ) )
            '''
    
    ## 3
    if kSilhouetteConstraints: add_silhouette_constraints(
        region,
        mesh, boundary_vertices, front_vertices, back_vertices,
        boundary_vi_to_curve_vi,
        general_linear_constraints
        )
    ## TODO: Step 4.
    
    ## 5
    if not kNewCrossSectionConstraints:
        xyz_cross_sections = add_cross_section_constraints(
            region,
            mesh, boundary_vertices, front_vertices, back_vertices,
            xy_verts, xy_pos,
            zrel_groups_vert, zrel_groups_vals,
            general_linear_constraints
            )
    else:
        xyz_cross_sections = add_cross_section_constraints2(
            region,
            mesh, boundary_vertices, front_vertices, back_vertices,
            xy_verts, xy_pos,
            zrel_groups_vert, zrel_groups_vals,
            xy_vertex_weighted_sum_indices, xy_vertex_weighted_sum_weights, xy_vertex_weighted_sum_rhss,
            z_vertex_weighted_sum_indices, z_vertex_weighted_sum_weights, z_vertex_weighted_sum_rhss
            )
    
    ## We need at least one (arbitrary) z value constraint.
    ## TODO Q: Does this cause a propagation problem because our system is nonlinear?
    if kOldZConstraint:
        z_verts = [boundary_vertices[0]]
        z_vals = [0.]
        ## Prefer an endpoint in zrel_groups if we have one.
        for verts in zrel_groups_vert:
            if len( verts ) > 0:
                zverts = [ verts[0] ]
                #print 'zrel_z_vert:', zverts[0]
                break
    else:
        z_verts = []
        z_vals = []
        
        ## Constrain all points to have average z = 0.
        #z_vertex_weighted_sum_indices.append( range(len(mesh.vs)) )
        #z_vertex_weighted_sum_weights.append( onesf( len(mesh.vs) ) )
        z_vertex_weighted_sum_indices.append( boundary_vertices )
        z_vertex_weighted_sum_weights.append( onesf( len(boundary_vertices) ) )
        z_vertex_weighted_sum_rhss.append( 0. )
    
    ## Actually, keep the silhouette orthogonal
    #z_verts = [ mvi for mvi in boundary_vertices if True not in [ mvi in verts for verts in zrel_groups_vert ] ]
    #z_vals = [0.] * len(z_verts)
    #print 'boundary_vertices:', boundary_vertices
    #print 'zverts:', zverts
    
    ### Convert zrel constraints into z constraints
    #for verts, relvals in izip( zrel_groups_vert, zrel_groups_vals ):
    #    for mvi, zrel in izip( verts, relvals ):
    #        if mvi not in z_verts:
    #            z_verts.append( mvi )
    #            z_vals.append( zrel )
    #zrel_groups_vert = []
    #zrel_groups_vals = []
    
    ## Constrain xy of all vertices, leaving a double-sided height field.
    '''
    for i in xrange(len( mesh.vs )):
        #if i not in xy_verts and i not in boundary_vertices:
        if i not in xy_verts and not contour_contains( region.region_boundary, mesh.vs[i][:2] ):
            pass
            #xy_verts.append( i )
            #xy_pos.append( mesh.vs[i] )
    '''
    
    constraints = Constraints()
    constraints.xy_constraint_vertices = xy_verts
    constraints.xy_constraint_positions = xy_pos
    constraints.z_constraint_vertices = z_verts
    constraints.z_constraint_values = z_vals
    constraints.z_relative_constrained_groups_vertices = zrel_groups_vert
    constraints.z_relative_constrained_groups_values = zrel_groups_vals
    constraints.xy_vertex_weighted_sum_indices = xy_vertex_weighted_sum_indices
    constraints.xy_vertex_weighted_sum_weights = xy_vertex_weighted_sum_weights
    constraints.xy_vertex_weighted_sum_rhss = xy_vertex_weighted_sum_rhss
    constraints.z_vertex_weighted_sum_indices = z_vertex_weighted_sum_indices
    constraints.z_vertex_weighted_sum_weights = z_vertex_weighted_sum_weights
    constraints.z_vertex_weighted_sum_rhss = z_vertex_weighted_sum_rhss
    constraints.general_linear_constraints = general_linear_constraints
    
    ## This should produce the same results as regular, orthogonal FiberMesh
    #constraints.z_constraint_vertices = xy_verts
    #constraints.z_constraint_values = [0.] * len(xy_verts)
    #constraints.z_relative_constrained_groups_vertices = []
    #constraints.z_relative_constrained_groups_values = []
    
    return ( mesh, constraints, xyz_cross_sections )

def mesh_llconstraints_from_constraints2_oblique( region, oblique ):
    '''
    Given a RotundRegion 'region' and a pair of numbers 'oblique' representing
    the x and y shift according to z ( x + oblique[0]*z, y + oblique[1]*z ),
    returns the tuple ( mesh, low-level constraints, list of 3d curves representing cross sections ).
    The first two elements are suitable for passing to fibermesh.optimize_surface2(),
    and the third element is a list of 3d cross sections curves corresponding to each 2d cross
    section in 'region.cross_sections'.
    
    This assumes that the curves have already been smoothed and self-intersections removed.
    
    tested
    '''
    
    print("[[[Entered mesh_llconstraints_from_constraints2_oblique()]]]")
    
    name = ''
    
    assert len( oblique ) == 2
    
    ## In this function, we are following steps 1-5 of inflate() in
    ## "annotations revisited data structures.txt":
    ## 1 glue two disks (2d triangle-meshes) together along 'region_boundary',
    ##   slightly offset from each other in z.
    ## 2 constrain 'region_boundary' in xy.
    ## 3 constrain normals of edges in 'silhouette_edges' to have z = 0.
    ## 4 for edges in 'sharp_cut_edges', duplicate the boundary vertices and disconnect the back and
    ##   front vertices from each other.  leave the front vertices constrained in xy, but
    ##   do not constrain the xy or z of the back vertices.
    ## 
    ## 5 for each cross_section in 'cross_sections', let 'cs' be the mesh vertices (front and back)
    ##   near cross_section.  (See the NOTE near the definition of the cross_sections in RotundRegion.)
    ##   Determine the desired tilt by rotating the cross section circle about the line between the
    ##   cross section's two points on the mesh's silhouette such that the anchor point
    ##   (cross_section[0]) is as close as possible to its target point (cross_section[1]).
    ##   On each vertex in 'cs', place an xy constraint at the corresponding image-plane location
    ##   on the tilted cross-section.  Also place a relative z constraint with the first vertex
    ##   in 'cs'.  (This is an arbitrary choice, but the relative constraints must be relative to
    ##   some vertex in the cross section's z).
    
    from fibermesh import gen_glued_mesh, attach_curve2d_to_mesh_vertices, Constraints
    
    ## 1
    mesh, boundary_vertices, front_vertices, back_vertices = gen_glued_mesh( region.region_boundary )
    ## These sets should be disjoint.
    assert len( set( boundary_vertices ).intersection( set( front_vertices ) ) ) == 0
    assert len( set( boundary_vertices ).intersection( set( back_vertices ) ) ) == 0
    assert len( set( front_vertices ).intersection( set( back_vertices ) ) ) == 0
    ## If we already have a mesh, re-use it:
    if region.mesh_no_cross_sections: mesh = region.mesh_no_cross_sections
    mesh.write_OBJ( name + 'glued.obj' )
    mesh.vs = asarray( mesh.vs )
    
    xy_verts = []
    xy_pos = []
    
    z_verts = []
    z_vals = []
    
    zrel_groups_vert = []
    zrel_groups_vals = []
    
    xy_vertex_weighted_sum_indices = []
    xy_vertex_weighted_sum_weights = []
    xy_vertex_weighted_sum_rhss = []
    z_vertex_weighted_sum_indices = []
    z_vertex_weighted_sum_weights = []
    z_vertex_weighted_sum_rhss = []
    
    general_linear_constraints = []
    
    #if not hasattr( region, '_decs' ):
    region._decs = []
    
    ## 2
    boundary_vertices, boundary_vi_to_curve_vi = attach_curve2d_to_mesh_vertices( region.region_boundary, mesh, boundary_vertices )
    #xy_verts = list( boundary_vertices )
    #xy_pos = [ region.region_boundary[ si ] for si in boundary_vi_to_curve_vi ]
    ## Only constrain boundary vertices that are endpoints of silhouette or sharp-cut edges.
    constrain_vertices = set()
    constrain_vertices.update( region.silhouette_edges )
    constrain_vertices.update( [ ( ei + 1 ) % len( region.region_boundary ) for ei in region.silhouette_edges ] )
    constrain_vertices.update( region.sharp_cut_edges )
    constrain_vertices.update( [ ( ei + 1 ) % len( region.region_boundary ) for ei in region.sharp_cut_edges ] )
    for bvi, cvi in zip( boundary_vertices, boundary_vi_to_curve_vi ):
        if cvi in constrain_vertices:
            for coord in (0,1):
                general_linear_constraints.append( (
                    ## The list of ( vertex index, coord ) tuples.
                    [ (bvi,coord), (bvi,2) ],
                    ## The list of coefficients for each term.
                    ## We add the oblique parameter * z because we want xy.screen = xy.world + ab*z.
                    [ 1., oblique[coord] ],
                    ## The right hand side value.
                    region.region_boundary[ cvi ][coord]
                    ) )
            
            ( ( #region._decs.append( (
                'line', (0,0,0),
                append( region.region_boundary[ cvi ], 0 ) + 1*array( (oblique[0], oblique[1], -1) ),
                append( region.region_boundary[ cvi ], 0 ) - 1*array( (oblique[0], oblique[1], -1) )
                ) )
    
    ## 3
    if kSilhouetteConstraints: add_silhouette_constraints(
        region,
        mesh, boundary_vertices, front_vertices, back_vertices,
        boundary_vi_to_curve_vi,
        general_linear_constraints,
        oblique
        )
    ## TODO: Step 4.
    
    ## 5
    xyz_cross_sections = add_cross_section_constraints2_oblique(
        region,
        mesh, boundary_vertices, front_vertices, back_vertices,
        xy_verts, xy_pos,
        zrel_groups_vert, zrel_groups_vals,
        xy_vertex_weighted_sum_indices, xy_vertex_weighted_sum_weights, xy_vertex_weighted_sum_rhss,
        z_vertex_weighted_sum_indices, z_vertex_weighted_sum_weights, z_vertex_weighted_sum_rhss,
        general_linear_constraints,
        oblique
        )
    
    ## We need at least one z value constraint.
    ## TODO Q: Does this cause a propagation problem because our system is nonlinear?
    kZConstraintType = 2
    if kZConstraintType == 1:
        z_verts = [boundary_vertices[0]]
        z_vals = [0.]
        ## Prefer an endpoint in zrel_groups if we have one.
        for verts in zrel_groups_vert:
            if len( verts ) > 0:
                zverts = [ verts[0] ]
                #print 'zrel_z_vert:', zverts[0]
                break
    elif kZConstraintType == 2:
        ## Find the boundary vertex with xy gradient closest to oblique direction.
        if mag2( array( oblique ) ) < kEps:
            z_verts = [boundary_vertices[0]]
            z_vals = [0.]
        else:
            oblique_dir = dir( array( oblique ) )
            
            silhouette_vertices = set()
            silhouette_vertices.update( region.silhouette_edges )
            silhouette_vertices.update( [ ( ei + 1 ) % len( region.region_boundary ) for ei in region.silhouette_edges ] )
            
            boundary_vertices, boundary_vi_to_curve_vi = attach_curve2d_to_mesh_vertices( region.region_boundary, mesh, boundary_vertices )
            max_cosangle = -1.
            max_cosangle_i = None
            for i in range(len( boundary_vi_to_curve_vi )):
                if boundary_vi_to_curve_vi[i] in silhouette_vertices:
                    rbN = len( region.region_boundary )
                    bvtangent = dir(
                        dir(
                            array( region.region_boundary[ ( boundary_vi_to_curve_vi[i] + 1 ) % rbN ] ) -
                            array( region.region_boundary[ boundary_vi_to_curve_vi[i] ] )
                        ) + dir(
                            array( region.region_boundary[ boundary_vi_to_curve_vi[i] ] ) -
                            array( region.region_boundary[ ( boundary_vi_to_curve_vi[i] - 1 + rbN ) % rbN ] )
                        )
                        )
                    cosangle = abs( dot( bvtangent, oblique_dir ) )
                    
                    if cosangle > max_cosangle:
                        max_cosangle = cosangle
                        max_cosangle_i = i
            
            assert max_cosangle_i is not None
            z_verts = [boundary_vertices[max_cosangle_i]]
            z_vals = [0.]
            print('constraint z = 0 for vertex', z_verts[0])
    elif kZConstraintType == 3:
        ## Q: Does it work if we don't constrain any z on its own?
        ## A: No.  It explodes immediately.
        pass
    else:
        z_verts = []
        z_vals = []
        
        ## Constrain all points to have average z = 0.
        #z_vertex_weighted_sum_indices.append( range(len(mesh.vs)) )
        #z_vertex_weighted_sum_weights.append( onesf( len(mesh.vs) ) )
        z_vertex_weighted_sum_indices.append( boundary_vertices )
        z_vertex_weighted_sum_weights.append( onesf( len(boundary_vertices) ) )
        z_vertex_weighted_sum_rhss.append( 0. )
    
    ## Actually, keep the silhouette orthogonal
    #z_verts = [ mvi for mvi in boundary_vertices if True not in [ mvi in verts for verts in zrel_groups_vert ] ]
    #z_vals = [0.] * len(z_verts)
    #print 'boundary_vertices:', boundary_vertices
    #print 'zverts:', zverts
    
    ### Convert zrel constraints into z constraints
    #for verts, relvals in izip( zrel_groups_vert, zrel_groups_vals ):
    #    for mvi, zrel in izip( verts, relvals ):
    #        if mvi not in z_verts:
    #            z_verts.append( mvi )
    #            z_vals.append( zrel )
    #zrel_groups_vert = []
    #zrel_groups_vals = []
    
    ## Constrain xy of all vertices, leaving a double-sided height field.
    '''
    for i in xrange(len( mesh.vs )):
        #if i not in xy_verts and i not in boundary_vertices:
        if i not in xy_verts and not contour_contains( region.region_boundary, mesh.vs[i][:2] ):
            pass
            #xy_verts.append( i )
            #xy_pos.append( mesh.vs[i] )
    '''
    
    constraints = Constraints()
    constraints.xy_constraint_vertices = xy_verts
    constraints.xy_constraint_positions = xy_pos
    constraints.z_constraint_vertices = z_verts
    constraints.z_constraint_values = z_vals
    constraints.z_relative_constrained_groups_vertices = zrel_groups_vert
    constraints.z_relative_constrained_groups_values = zrel_groups_vals
    constraints.xy_vertex_weighted_sum_indices = xy_vertex_weighted_sum_indices
    constraints.xy_vertex_weighted_sum_weights = xy_vertex_weighted_sum_weights
    constraints.xy_vertex_weighted_sum_rhss = xy_vertex_weighted_sum_rhss
    constraints.z_vertex_weighted_sum_indices = z_vertex_weighted_sum_indices
    constraints.z_vertex_weighted_sum_weights = z_vertex_weighted_sum_weights
    constraints.z_vertex_weighted_sum_rhss = z_vertex_weighted_sum_rhss
    constraints.general_linear_constraints = general_linear_constraints
    
    ## This should produce the same results as regular, orthogonal FiberMesh
    #constraints.z_constraint_vertices = xy_verts
    #constraints.z_constraint_values = [0.] * len(xy_verts)
    #constraints.z_relative_constrained_groups_vertices = []
    #constraints.z_relative_constrained_groups_values = []
    
    return ( mesh, constraints, xyz_cross_sections )

def add_cross_section_constraints(
    region,
    mesh, boundary_vertices, front_vertices, back_vertices,
    xy_verts, xy_pos,
    zrel_groups_vert, zrel_groups_vals,
    linear_constraints
    ):
    '''
    Implements Step 5 from inflate() in "annotations revisited data structures.txt" in the
    following manner:
    a We'll find a straight line for each cross section by intersecting with the mesh the line
      passing through the anchor point and perpendicular to the line from it to its destination
      location.
    b Compute the tilted circle implied by the cross section and the anchor point's destination.
    c For intersected mesh edges, assign the edges' vertices to lie along the tilted circle
      using xy and relative z constraints.
    d Constrain edge normals to have the circles' normals.
    
    Returns a list of line loops corresponding to the xyz cross sections.
    
    tested
    '''
    
    if not kCrossSectionConstraints1ConstrainsEdgeNormals:
        ## We skip these constraints by renaming output parameter 'linear_constraints'.
        linear_constraints = []
    
    xyz_cross_sections = []
    for cross_section in region.cross_sections:
        ## a
        anchor = asarrayf( cross_section[0] )
        ## Find the direction perpendicular to the spine at the anchor point.
        anchor_to_dest = asarrayf( cross_section[1] ) - anchor
        anchor_to_dest_perp = anchor_to_dest.copy()
        anchor_to_dest_perp[0], anchor_to_dest_perp[1] = anchor_to_dest_perp[1], -anchor_to_dest_perp[0]
        if mag2( anchor_to_dest_perp ) < kEps:
            print("WARNING: No spine and anchor didn't move.  Using <0,1> as the spine_perp.")
            anchor_to_dest_perp[:] = (0,1)
        boundary1, cs_front_vertices, boundary2, cs_back_vertices = get_cross_section_vertices(
            anchor, anchor_to_dest_perp,
            mesh, boundary_vertices, front_vertices, back_vertices
            )
        expect( len( cs_front_vertices ) > 0 )
        expect( len( cs_back_vertices ) > 0 )
        
        ## b
        ## Give the two boundary vertices constant z.  Note that this z floats,
        ## so optimization may change it.
        ## To make my asserts easier, and because the cross sections float,
        ## center cross sections on the z = 0. plane.
        mesh.vs[ boundary1 ][2] = 0.
        mesh.vs[ boundary2 ][2] = mesh.vs[ boundary1 ][2]
        ## Get the cross section's tilt.
        #center, plane_perp = get_tilted_circle_plane_for_cross_section(
        radians = get_rotation_about_boundary1to2axis_for_cross_section(
            mesh.vs[ boundary1 ],
            mesh.vs[ boundary2 ],
            anchor_to_dest
            )
        cross_section_radius = .5 * len( mesh.vs[ boundary2 ] - mesh.vs[ boundary1 ] )
        boundary1to2axis = ( mesh.vs[ boundary2 ] - mesh.vs[ boundary1 ] ) / (2*cross_section_radius)
        ## The two boundary vertices have the same z, so the vector between them should have z = 0.
        assert abs( boundary1to2axis[2] ) < kEps
        center = .5 * ( mesh.vs[ boundary2 ] + mesh.vs[ boundary1 ] )
        
        ## c
        u = rotate( (0,0,1), boundary1to2axis, radians )
        ## 'u' should never rotate towards the back of the model.
        assert u[2] >= 0.
        ## Make 'u' the same length as boundary1to2axis:
        u = mag( boundary1to2axis ) * dir( u )
        zrel_groups_vert.append( [] )
        zrel_groups_vals.append( [] )
        
        ## For circular cross sections, the two boundary vertices have the same z.
        ## NOTE: Cross sections must touch silhouette boundaries.
        zrel_groups_vert[-1].append( boundary1 )
        zrel_groups_vals[-1].append( center[2] )
        zrel_groups_vert[-1].append( boundary2 )
        zrel_groups_vals[-1].append( center[2] )
        
        ## NOTE: All boundary vertices have already been fixed in xy.
        
        xyz_cross_sections.append( [] )
        ## Make a copy of mesh.vs[ boundary1 and 2 ] because otherwise they follow the mesh while
        ## the rest of the points in xyz_cross_sections don't.
        xyz_cross_sections[-1].append( list( mesh.vs[ boundary1 ] ) )
        
        last_vertex_index = boundary1
        last_vertex_pos = mesh.vs[ boundary1 ]
        N = len( cs_front_vertices )
        for i, fvi in enumerate( cs_front_vertices ):
            theta = pi * float(i + 1) / ( N+1 )
            assert theta > 0 and theta < pi
            pos = center + cross_section_radius * ( -cos( theta ) * boundary1to2axis + sin( theta ) * u )
            assert pos[2] >= 0. or debugger()
            ## Don't move these vertices in xy.
            xy_verts.append( fvi )
            xy_pos.append( pos[:2] )
            ## In z, move them rigidly with the other vertices in the cross section.
            zrel_groups_vert[-1].append( fvi )
            zrel_groups_vals[-1].append( pos[2] )
            
            ## d
            if (i % kCrossSectionConstraints1ConstrainsEdgeNormalsSubset) == 0:
                linear_constraints.extend(
                    gen_normal_constraint( mesh, dir( pos - center ), last_vertex_index, fvi, pos - last_vertex_pos )
                    )
            last_vertex_index = fvi
            last_vertex_pos = pos
            
            xyz_cross_sections[-1].append( pos )
        
        ## d
        if (N % kCrossSectionConstraints1ConstrainsEdgeNormalsSubset) == 0:
            linear_constraints.extend(
                gen_normal_constraint( mesh, dir( pos - center ), last_vertex_index, boundary2, mesh.vs[ boundary2 ] - last_vertex_pos )
                )
        ## Make a copy of mesh.vs[ boundary1 and 2 ] because otherwise they follow the mesh while
        ## the rest of the points in xyz_cross_sections don't.
        xyz_cross_sections[-1].append( list( mesh.vs[ boundary2 ] ) )
        
        last_vertex_index = boundary2
        last_vertex_pos = mesh.vs[ boundary2 ]
        N = len( cs_back_vertices )
        for i, fvi in enumerate( cs_back_vertices ):
            theta = pi + pi * float(i + 1) / ( N+1 )
            assert theta > pi and theta < 2*pi
            pos = center + cross_section_radius * ( -cos( theta ) * boundary1to2axis + sin( theta ) * u )
            assert pos[2] <= 0. or debugger()
            ## Don't move these vertices in xy.
            xy_verts.append( fvi )
            xy_pos.append( pos[:2] )
            ## In z, move them rigidly with the other vertices in the cross section.
            zrel_groups_vert[-1].append( fvi )
            zrel_groups_vals[-1].append( pos[2] )
            
            ## d
            if ((i + len( cs_front_vertices ) + 1) % kCrossSectionConstraints1ConstrainsEdgeNormalsSubset) == 0:
                linear_constraints.extend(
                gen_normal_constraint( mesh, dir( pos - center ), last_vertex_index, fvi, pos - last_vertex_pos )
                )
            last_vertex_index = fvi
            last_vertex_pos = pos
            
            xyz_cross_sections[-1].append( pos )
        
        ## d
        if ((len( cs_front_vertices ) + 1 + N + 1) % kCrossSectionConstraints1ConstrainsEdgeNormalsSubset) == 0:
            linear_constraints.extend(
                gen_normal_constraint( mesh, dir( pos - center ), last_vertex_index, boundary1, mesh.vs[ boundary1 ] - last_vertex_pos )
                )
    
    return xyz_cross_sections


def add_cross_section_constraints2(
    region,
    mesh, boundary_vertices, front_vertices, back_vertices,
    xy_verts, xy_pos,
    zrel_groups_vert, zrel_groups_vals,
    xy_vertex_weighted_sum_indices, xy_vertex_weighted_sum_weights, xy_vertex_weighted_sum_rhss,
    z_vertex_weighted_sum_indices, z_vertex_weighted_sum_weights, z_vertex_weighted_sum_rhss
    ):
    '''
    Implements Step 5 from inflate() in "annotations revisited data structures.txt" in the
    following manner:
    a We'll find a straight line for each cross section by intersecting with the mesh the line
      passing through the anchor point and perpendicular to the line from it to its destination
      location.
    b Compute the tilted circle implied by the cross section and the anchor point's destination.
    c For intersected mesh edges, assign the point of intersection (as a weighted sum of the edges'
      vertices) to lie along the tilted circle using xy and relative z constraints.
      Also constrain the edge to be perpendicular to the tilted circle's normal
      (use the cross section's rotation? only if the edges haven't been tilted yet.).
    
    Returns a list of line loops corresponding to the xyz cross sections.
    
    used
    '''
    
    constrained_edges = set()
    
    xyz_cross_sections = []
    for cross_section in region.cross_sections:
        xyz_cross_sections.append( [] )
        
        ## a
        anchor = asarrayf( cross_section[0] )
        ## Find the direction perpendicular to the spine at the anchor point.
        anchor_to_dest = asarrayf( cross_section[1] ) - anchor
        anchor_to_dest_perp = anchor_to_dest.copy()
        anchor_to_dest_perp[0], anchor_to_dest_perp[1] = anchor_to_dest_perp[1], -anchor_to_dest_perp[0]
        if mag2( anchor_to_dest_perp ) < kEps:
            print("WARNING: No spine and anchor didn't move.  Using <0,1> as the spine_perp.")
            anchor_to_dest_perp[:] = (0,1)
        boundary1, cs_front_edges, boundary2, cs_back_edges = get_cross_section_edges(
            anchor, anchor_to_dest_perp,
            mesh, boundary_vertices, front_vertices, back_vertices
            )
        expect( len( cs_front_edges ) > 0 )
        expect( len( cs_back_edges ) > 0 )
        
        '''
        print 'boundary1:', boundary1
        print 'boundary2:', boundary2
        print 'cs_front_edges:', cs_front_edges
        print 'cs_back_edges:', cs_back_edges
        '''
        
        cross_section_all_edges = set()
        cross_section_all_edges.add( frozenset( boundary1[0] ) )
        cross_section_all_edges.add( frozenset( boundary2[0] ) )
        cross_section_all_edges.update( [ frozenset( edge_t[0] ) for edge_t in cs_front_edges ] )
        cross_section_all_edges.update( [ frozenset( edge_t[0] ) for edge_t in cs_back_edges ] )
        if len( constrained_edges.intersection( cross_section_all_edges ) ) > 0:
            print('cross section intersects another cross section; skipping:', cross_section)
            continue
        
        constrained_edges.update( cross_section_all_edges )
        
        ## Remove the first and last edges, which connect to boundary vertices.
        ## I don't know if this is necessary, but my constraints are seeming too harsh,
        ## so perhaps this will help.
        ## TODO: Revisit this.
        cs_front_edges = cs_front_edges[1:-1]
        cs_back_edges = cs_back_edges[1:-1]
        
        region._decs = []
        
        ## b
        boundary1pos = lerp(
            mesh.vs[ boundary1[0][0] ], mesh.vs[ boundary1[0][1] ],
            boundary1[1]
            )
        boundary2pos = lerp(
            mesh.vs[ boundary2[0][0] ], mesh.vs[ boundary2[0][1] ],
            boundary2[1]
            )
        ## To make my asserts easier, and because the cross sections float,
        ## center cross sections on the z = 0. plane.
        boundary1pos[2] = 0.
        boundary2pos[2] = boundary1pos[2]
        ## Get the cross section's tilt.
        #center, plane_perp = get_tilted_circle_plane_for_cross_section(
        radians = get_rotation_about_boundary1to2axis_for_cross_section(
            boundary1pos,
            boundary2pos,
            anchor_to_dest
            )
        cross_section_radius = .5 * len( boundary2pos - boundary1pos )
        boundary1to2axis = ( boundary2pos - boundary1pos ) / (2*cross_section_radius)
        ## The two boundary vertices are constrained to have the same z,
        ## so the vector between them should have z = 0.
        assert abs( boundary1to2axis[2] ) < kEps
        center = .5 * ( boundary2pos + boundary1pos )
        
        ## c
        u = rotate( (0,0,1), boundary1to2axis, radians )
        ## 'u' should never rotate towards the back of the model.
        assert u[2] >= 0.
        ## Make 'u' the same length as boundary1to2axis:
        u = mag( boundary1to2axis ) * dir( u )
        
        ## We will now walk through all intersected edges and
        ##   1 constrain the edge intersection point in xy and relative z.
        ##   2 constrain the edge vector to a rotated version of itself, where the rotation is the
        ##      minimal rotation such that the dot product with the circle normal is 0.
        
        ## NOTES:
        ##   1 This causes some mesh faces to self-intersect.  Perhaps when I'm constraining edges
        ##     that run nearly parallel to the cross section line?  No, it still happens when I
        ##     skip nearly parallel edges.
        ##   2 This works and could be the solution on its own.  However, it is very dependent
        ##     on how inflated the mesh is.  If I don't allow the surface to inflate all the way,
        ##     the edges rotate and are constrained in what looks like too small a circle.
        ##     Perhaps because I am spacing things out uniformly, instead of projecting
        ##     intersections onto the closest points of the circle.
        ##   Perhaps the best thing to do is use the geodesic shortest path curve as in
        ##   add_cross_section_constraints() and add constraints to the vertices opposite the edges
        ##   so that the normal points in the direction I want.  This is the same as my proposed
        ##   silhouette-enforcing constraints on the boundary curve.
        
        if kCrossSectionConstraints2ConstrainsIntersectionPoints:
            ## 1
            ## No need to constrain the boundary1 and boundary2 edges in xy, because their endpoints
            ## are already constrained.
            ## Still constrain their z values, however.
            ## NOTE: We have no relative z_vertex_weighted_sum constraints, so we'll have to do it
            ##       manually.  Use boundary1 as the reference point.  This means we only constrain
            ##       boundary2's z value relatively then.
            
            ## This is what an absolute constraint on boundary1 would look like:
            # z_vertex_weighted_sum_indices.append( boundary1[0] )
            # z_vertex_weighted_sum_weights.append( ( boundary1[1], 1. - boundary1[1] ) )
            assert boundary1[1] > -kEps and boundary1[1] < 1. + kEps
            # z_vertex_weighted_sum_rhss.append( 0. )
            
            ## We'll make a similar constraint on boundary2 and subtract it from boundary1's constraint.
            z_vertex_weighted_sum_indices.append( list( boundary2[0] ) + list( boundary1[0] ) )
            z_vertex_weighted_sum_weights.append( (
                boundary2[1], 1. - boundary2[1], - boundary1[1], -( 1. - boundary1[1] ),
                ) )
            assert boundary2[1] > -kEps and boundary2[1] < 1. + kEps
            z_vertex_weighted_sum_rhss.append( 0. )
        
        xyz_cross_sections[-1].append( boundary1pos )
        
        ## Now add xy and relative z constraints on the intersections in front_edges, which lead
        ## from boundary1 to boundary2.
        N = len( cs_front_edges )
        for i, edge_t in enumerate( cs_front_edges ):
            if dot( dir( boundary1to2axis[:2] ), dir( ( mesh.vs[ edge_t[0][0] ] - mesh.vs[ edge_t[0][1] ] )[:2] ) ) > .5:
                print('skipping parallel edge')
                continue
            
            region._decs.append( (
                'line',
                ( 1, 0, 0 ),
                mesh.vs[ edge_t[0][0] ],
                mesh.vs[ edge_t[0][1] ]
                ) )
            
            theta = pi * float(i + 1) / ( N+1 )
            assert theta > 0 and theta < pi
            pos = center + cross_section_radius * ( -cos( theta ) * boundary1to2axis + sin( theta ) * u )
            assert pos[2] >= 0. or debugger()
            ## Set the intersection position's xy.
            if kCrossSectionConstraints2ConstrainsIntersectionPoints:
                xy_vertex_weighted_sum_indices.append( edge_t[0] )
                xy_vertex_weighted_sum_weights.append( ( edge_t[1], 1. - edge_t[1] ) )
                xy_vertex_weighted_sum_rhss.append( pos[:2] )
                assert edge_t[1] > -kEps and edge_t[1] < 1. + kEps
                
                ## In z, move them rigidly with the other vertices in the cross section.
                z_vertex_weighted_sum_indices.append( list( edge_t[0] ) + list( boundary1[0] ) )
                z_vertex_weighted_sum_weights.append( (
                    edge_t[1], 1. - edge_t[1], - boundary1[1], -( 1. - boundary1[1] ),
                    ) )
                z_vertex_weighted_sum_rhss.append( pos[2] )
            
            if kCrossSectionConstraints2ConstrainsEdgeRotation:
                ## 2
                ## Finally, rotate the edge so that it's on the plane perpendicular to the normal
                ## vector dir(pos-center).
                ## Using the intersection point of the edge as our origin, project the edge
                ## onto the plane passing through the origin with normal dir(pos-center).
                ## Rescale the edge so its length is preserved, and make that the constraint.
                edge0 = mesh.vs[ edge_t[0][0] ]
                edge1 = mesh.vs[ edge_t[0][1] ]
                edgepos = lerp(
                    edge0, edge1,
                    edge_t[1]
                    )
                normal = dir( pos - center )
                edge0proj = project_point_to_plane( edge0, ( edgepos, normal ) )
                edge1proj = project_point_to_plane( edge1, ( edgepos, normal ) )
                edge0proj = mag( edge0 - edgepos ) * dir( edge0proj - edgepos ) + edgepos
                edge1proj = mag( edge1 - edgepos ) * dir( edge1proj - edgepos ) + edgepos
                ## Constrain the edge to this vector.
                xy_vertex_weighted_sum_indices.append( edge_t[0] )
                xy_vertex_weighted_sum_weights.append( [1, -1] )
                xy_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[:2] )
                z_vertex_weighted_sum_indices.append( edge_t[0] )
                z_vertex_weighted_sum_weights.append( [1, -1] )
                z_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[2] )
                
                region._decs.append( (
                    'line',
                    ( .5, .5, 0 ),
                    edge0proj,
                    edge1proj
                    ) )
            
            if kCrossSectionConstraints2ConstrainsEdgeRotation2:
                ## 2
                ## Finally, rotate the edge so that it's on the plane perpendicular to the normal
                ## vector dir(pos-center).
                ## Using the intersection point of the edge as our origin, project the edge
                ## onto the plane passing through the origin with normal dir(pos-center).
                ## Rescale the edge so its length is preserved, and make that the constraint.
                edge0 = mesh.vs[ edge_t[0][0] ]
                edge1 = mesh.vs[ edge_t[0][1] ]
                edge0proj = center + rotate( edge0 - center, boundary1to2axis, radians )
                edge1proj = center + rotate( edge1 - center, boundary1to2axis, radians )
                ## Constrain the edge to this vector.
                xy_vertex_weighted_sum_indices.append( edge_t[0] )
                xy_vertex_weighted_sum_weights.append( [1, -1] )
                xy_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[:2] )
                z_vertex_weighted_sum_indices.append( edge_t[0] )
                z_vertex_weighted_sum_weights.append( [1, -1] )
                z_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[2] )
                
                region._decs.append( (
                    'line',
                    ( .5, .5, 0 ),
                    edge0proj,
                    edge1proj
                    ) )
            
            xyz_cross_sections[-1].append( pos )
        
        xyz_cross_sections[-1].append( boundary2pos )
        
        ## Now add xy and relative z constraints on the intersections in back_edges, which lead
        ## from boundary2 to boundary1.
        N = len( cs_back_edges )
        for i, edge_t in enumerate( cs_back_edges ):
            if dot( dir( boundary1to2axis[:2] ), dir( ( mesh.vs[ edge_t[0][0] ] - mesh.vs[ edge_t[0][1] ] )[:2] ) ) > .5:
                print('skipping parallel edge')
                continue
            
            region._decs.append( (
                'line',
                ( 1, 0, 0 ),
                mesh.vs[ edge_t[0][0] ],
                mesh.vs[ edge_t[0][1] ]
                ) )
            
            theta = pi + pi * float(i + 1) / ( N+1 )
            assert theta > pi and theta < 2*pi
            pos = center + cross_section_radius * ( -cos( theta ) * boundary1to2axis + sin( theta ) * u )
            assert pos[2] <= 0. or debugger()
            ## Set the intersection position's xy.
            if kCrossSectionConstraints2ConstrainsIntersectionPoints:
                xy_vertex_weighted_sum_indices.append( edge_t[0] )
                xy_vertex_weighted_sum_weights.append( ( edge_t[1], 1. - edge_t[1] ) )
                xy_vertex_weighted_sum_rhss.append( pos[:2] )
                assert edge_t[1] > -kEps and edge_t[1] < 1. + kEps
                
                ## In z, move them rigidly with the other vertices in the cross section.
                z_vertex_weighted_sum_indices.append( list( edge_t[0] ) + list( boundary1[0] ) )
                z_vertex_weighted_sum_weights.append( (
                    edge_t[1], 1. - edge_t[1], - boundary1[1], -( 1. - boundary1[1] ),
                    ) )
                z_vertex_weighted_sum_rhss.append( pos[2] )
            
            if kCrossSectionConstraints2ConstrainsEdgeRotation:
                ## 2
                ## Finally, rotate the edge so that it's on the plane perpendicular to the normal
                ## vector dir(pos-center).
                ## Using the intersection point of the edge as our origin, project the edge
                ## onto the plane passing through the origin with normal dir(pos-center).
                ## Rescale the edge so its length is preserved, and make that the constraint.
                edge0 = mesh.vs[ edge_t[0][0] ]
                edge1 = mesh.vs[ edge_t[0][1] ]
                edgepos = lerp(
                    edge0, edge1,
                    edge_t[1]
                    )
                normal = dir( pos - center )
                edge0proj = project_point_to_plane( edge0, ( edgepos, normal ) )
                edge1proj = project_point_to_plane( edge1, ( edgepos, normal ) )
                edge0proj = mag( edge0 - edgepos ) * dir( edge0proj - edgepos ) + edgepos
                edge1proj = mag( edge1 - edgepos ) * dir( edge1proj - edgepos ) + edgepos
                ## Constrain the edge to this vector.
                xy_vertex_weighted_sum_indices.append( edge_t[0] )
                xy_vertex_weighted_sum_weights.append( [1, -1] )
                xy_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[:2] )
                z_vertex_weighted_sum_indices.append( edge_t[0] )
                z_vertex_weighted_sum_weights.append( [1, -1] )
                z_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[2] )
                
                region._decs.append( (
                    'line',
                    ( .5, .5, 0 ),
                    edge0proj,
                    edge1proj
                    ) )
            
            if kCrossSectionConstraints2ConstrainsEdgeRotation2:
                ## 2
                ## Finally, rotate the edge so that it's on the plane perpendicular to the normal
                ## vector dir(pos-center).
                ## Using the intersection point of the edge as our origin, project the edge
                ## onto the plane passing through the origin with normal dir(pos-center).
                ## Rescale the edge so its length is preserved, and make that the constraint.
                edge0 = mesh.vs[ edge_t[0][0] ]
                edge1 = mesh.vs[ edge_t[0][1] ]
                edge0proj = center + rotate( edge0 - center, boundary1to2axis, radians )
                edge1proj = center + rotate( edge1 - center, boundary1to2axis, radians )
                ## Constrain the edge to this vector.
                xy_vertex_weighted_sum_indices.append( edge_t[0] )
                xy_vertex_weighted_sum_weights.append( [1, -1] )
                xy_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[:2] )
                z_vertex_weighted_sum_indices.append( edge_t[0] )
                z_vertex_weighted_sum_weights.append( [1, -1] )
                z_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[2] )
                
                region._decs.append( (
                    'line',
                    ( .5, .5, 0 ),
                    edge0proj,
                    edge1proj
                    ) )
            
            xyz_cross_sections[-1].append( pos )
    
    
    if kCrossSectionConstraints2ConstrainsCullSome:
        ## Cull half of my constraints.
        xy_vertex_weighted_sum_indices[:] = xy_vertex_weighted_sum_indices[::2]
        xy_vertex_weighted_sum_weights[:] = xy_vertex_weighted_sum_weights[::2]
        xy_vertex_weighted_sum_rhss[:] = xy_vertex_weighted_sum_rhss[::2]
        
        z_vertex_weighted_sum_indices[:] = z_vertex_weighted_sum_indices[::2]
        z_vertex_weighted_sum_weights[:] = z_vertex_weighted_sum_weights[::2]
        z_vertex_weighted_sum_rhss[:] = z_vertex_weighted_sum_rhss[::2]
        '''
        ## Simplify by converting to normal xy and z constraints.
        ## This doesn't go crazy.
        for (indices, weights), xyrhs in izip( izip( xy_vertex_weighted_sum_indices, xy_vertex_weighted_sum_weights ), xy_vertex_weighted_sum_rhss ):
            if indices[0] in xy_verts: continue
            xy_verts.append( indices[0] )
            xy_pos.append( xyrhs )
        
        del xy_vertex_weighted_sum_indices[:]
        del xy_vertex_weighted_sum_weights[:]
        del xy_vertex_weighted_sum_rhss[:]
        '''
        
        ## Simplify by making the constraints not relative.
        xy_vertex_weighted_sum_indices[:] = list(map( list, xy_vertex_weighted_sum_indices ))
        xy_vertex_weighted_sum_weights[:] = list(map( list, xy_vertex_weighted_sum_weights ))
        for (indices, weights), xyrhs in zip( zip( xy_vertex_weighted_sum_indices, xy_vertex_weighted_sum_weights ), xy_vertex_weighted_sum_rhss ):
            '''
            ## Make them single vertex constraints.
            del indices[1:]
            del weights[1:]
            weights[0] = 1.
            '''
            del indices[2:]
            del weights[2:]
            assert abs( sum( weights ) - 1. ) < kEps
        
        for (indices, weights), zrhs in zip( zip( z_vertex_weighted_sum_indices, z_vertex_weighted_sum_weights ), z_vertex_weighted_sum_rhss ):
            pass
            #z_verts.append( indices[0] )
            #z_pos.append( zrhs )
        
        del z_vertex_weighted_sum_indices[:]
        del z_vertex_weighted_sum_weights[:]
        del z_vertex_weighted_sum_rhss[:]
    
    return xyz_cross_sections

def decorate_ellipse( decs, el, color ):
    NS = 40
    decs.append( [
        'line',
        color,
        el[0] + cos( 0. ) * el[1][0] + sin( 0. ) * el[1][1]
        ] )
    for i in range(NS):
        theta = 2*pi* float(i+1) / NS
        decs[-1].append(
            el[0] + cos( theta ) * el[1][0] + sin( theta ) * el[1][1]
            )

def add_cross_section_constraints2_oblique( *args, **kwargs ):
    #return add_cross_section_constraints2_obliqueD( *args, **kwargs )
    return add_cross_section_constraints2_obliqueA( *args, **kwargs )

def add_cross_section_constraints2_obliqueD(
    region,
    mesh, boundary_vertices, front_vertices, back_vertices,
    xy_verts, xy_pos,
    zrel_groups_vert, zrel_groups_vals,
    xy_vertex_weighted_sum_indices, xy_vertex_weighted_sum_weights, xy_vertex_weighted_sum_rhss,
    z_vertex_weighted_sum_indices, z_vertex_weighted_sum_weights, z_vertex_weighted_sum_rhss,
    general_linear_constraints,
    oblique
    ):
    '''
    Implements Step 5 from inflate() in "annotations revisited data structures.txt" in the
    following manner:
    a We'll find a straight line for each cross section by intersecting with the mesh the line
      passing through the anchor point and perpendicular to the line from it to its destination
      location.
    b Compute the rotation axis and radians implied by the anchor point and the anchor point's destination.
    c For intersected mesh edges, constrain the edge to be the edge rotated by the axis and radians.
    
    Returns a list of line loops corresponding to the xyz cross sections.
    
    used
    '''
    
    if not hasattr( region, '_decs' ): region._decs = []
    
    
    constrained_edges = set()
    
    xyz_cross_sections = []
    for cross_section in region.cross_sections:
        xyz_cross_sections.append( [] )
        
        ## a
        anchor = asarrayf( cross_section[0] )
        ## Find the direction perpendicular to the spine at the anchor point.
        anchor_to_dest = asarrayf( cross_section[1] ) - anchor
        anchor_to_dest_perp = anchor_to_dest.copy()
        anchor_to_dest_perp[0], anchor_to_dest_perp[1] = anchor_to_dest_perp[1], -anchor_to_dest_perp[0]
        if mag2( anchor_to_dest_perp ) < kEps:
            print("WARNING: No spine and anchor didn't move.  Using <0,1> as the spine_perp.")
            anchor_to_dest_perp[:] = (0,1)
        
        ## TODO Q: What is a cross section under oblique projection?
        ##         Candidates:
        ##             A Use slices of the oblique skewed mesh that project obliquely to the same xy.
        ##               Pros: Similar to / Degenerates to the orthographic case.
        ##               Cons: The slice of the mesh is unintuitive.  The slices are skewed
        ##                     opposite the oblique skew (lines which are in the -skew direction
        ##                     project to points), which is strange for e.g. the front face of a
        ##                     smoothed cuboid with xyz edge lengths <1,1,3>.
        ##             B Use slices of the oblique skewed mesh that project orthographically to the same xy.
        ##               This is the same as option A.
        ##             C Find the point on the mesh under the xy pos (invert the oblique projection).
        ##               Then use the slice of the mesh that intersects the point and has the
        ##               same skew as the oblique projection.
        ##               Pros: Similar to / Degenerates to the orthographic case.  The slice of the
        ##                     mesh will be the same when clicking on a part that would be visible
        ##                     under orthographic projection.
        ##               Cons: The slice of the mesh will be the non-existent when clicking on a
        ##                     part that would not be visible under orthographic projection.
        ##             D Find the point on the mesh under the xy pos (invert the oblique projection).
        ##               Then use the normal of the mesh and the destination point projected
        ##               onto the tangent plane to find a cross section slice of the mesh.
        ##               Pros: Intuitive slice, if the existing 3d structure is intuitive.
        ##               Cons: The slices are not the same slices as in the orthographic case.
        ##                     Q: Can I remove the component of the normals that projects into the
        ##                        viewing plane?
        ##                     A: No.  In the ortho case this will always give me <0,0,1>.
        ##                        In the oblique case this will always give me <a,b,-1>.
        get_cross_section_edgesD_extra = {}
        boundary1, cs_front_edges, boundary2, cs_back_edges = get_cross_section_edgesD(
            anchor, anchor_to_dest_perp,
            mesh, boundary_vertices, front_vertices, back_vertices,
            oblique,
            get_cross_section_edgesD_extra
            )
        expect( len( cs_front_edges ) > 0 )
        expect( len( cs_back_edges ) > 0 )
        
        kNormalDecorationMag = .2
        region._decs.append( (
                'line',
                ( 1, 0, 0 ),
                get_cross_section_edgesD_extra['tangent_plane'][0],
                get_cross_section_edgesD_extra['tangent_plane'][0] + kNormalDecorationMag * get_cross_section_edgesD_extra['tangent_plane'][1]
                ) )
        
        region._decs.append( (
                'line',
                ( .5, 0, .5 ),
                get_cross_section_edgesD_extra['tangent_plane'][0],
                get_cross_section_edgesD_extra['tangent_plane'][0] + kNormalDecorationMag * get_cross_section_edgesD_extra['spine_perp_proj']
                ) )
        
        '''
        print 'boundary1:', boundary1
        print 'boundary2:', boundary2
        print 'cs_front_edges:', cs_front_edges
        print 'cs_back_edges:', cs_back_edges
        '''
        
        cross_section_all_edges = set()
        cross_section_all_edges.add( frozenset( boundary1[0] ) )
        cross_section_all_edges.add( frozenset( boundary2[0] ) )
        cross_section_all_edges.update( [ frozenset( edge_t[0] ) for edge_t in cs_front_edges ] )
        cross_section_all_edges.update( [ frozenset( edge_t[0] ) for edge_t in cs_back_edges ] )
        if len( constrained_edges.intersection( cross_section_all_edges ) ) > 0:
            print('cross section intersects another cross section; skipping:', cross_section)
            continue
        
        constrained_edges.update( cross_section_all_edges )
        
        ## b
        ## Before rotation:
        xyz_cross_sections[-1] = [
            lerp( mesh.vs[ edge_t[0][0] ], mesh.vs[ edge_t[0][1] ], edge_t[1] )
            for edge_t in [boundary1] + cs_front_edges + [boundary2] + cs_back_edges
            ]
        
        total_length = 0.
        center = zerosf(3)
        for i in range(len( xyz_cross_sections[-1] )):
            cspts = xyz_cross_sections[-1]
            
            pti = cspts[i]
            ptim1 = cspts[ ( i + len( cspts ) - 1 ) % len( cspts ) ]
            ptip1 = cspts[ ( i + 1 ) % len( cspts ) ]
            
            length = .5 * ( len( pti - ptim1 ) + len( pti - ptip1 ) )
            total_length += length
            center += length * pti
        center /= total_length
        axis = dir( get_cross_section_edgesD_extra['spine_perp_proj'] )
        start = get_cross_section_edgesD_extra['tangent_plane'][0]
        ## To find the radians, we want to find the closest point on the circle traced by 'start'
        ## around 'axis' to the oblique projection ray under 'anchor + anchor_to_dest'.
        anchor_dest3d = append( anchor + anchor_to_dest, 0 )
        anchor_dest3d_dir = arrayf( ( oblique[0], oblique[1], -1 ) )
        ellipse = (
            center,
            [
                start - center,
                mag( start - center ) * dir( cross( axis, start - center ) )
            ]
            )
        cos_theta, sin_theta = closest_point_coeffs_on_tilted_ellipse2d_to_line(
            ellipse, ( anchor_dest3d, anchor_dest3d + anchor_dest3d_dir )
            )
        region._decs.append( (
            'point',
            ( 0, 0, 0 ),
            ellipse2d_pt
            ) )
        region._decs.append( (
            'point',
            ( .7, .2, .7 ),
            ellipse2d[0] + cos_theta * ellipse2d[1][0] + sin_theta * ellipse2d[1][1]
            ) )
        region._decs.append( (
            'point',
            ( .2, .2, .7 ),
            ellipse[0] + cos_theta * ellipse[1][0] + sin_theta * ellipse[1][1]
            ) )
        from math import atan2
        ## TODO Q: Hmm, for an orthogonal projection the above code finds the points rotated
        ##         >pi/2 radians.
        ## A1: I could clip it, but that would just clip to pi/2 all the time.
        ## A2: Because 0 radians is properly set up to mean no rotation, I can just use the
        ##     absolute value of cos_theta, and that will always pick the point on the front half.
        ## atan2() with |cos| is the same as asin()
        #radians = atan2( sin_theta, abs( cos_theta ) )
        radians = asin( sin_theta )
        print("cos( theta ):", cos_theta, "; sin( theta ):", sin_theta, "; radians:", radians)
        region._decs.append( (
            'point',
            ( 0, 0, .4 ),
            ellipse[0] + abs( cos_theta ) * ellipse[1][0] + sin_theta * ellipse[1][1]
            ) )
        #debugger()
        
        ## After rotation:
        xyz_cross_sections[-1] = [
            center + rotate( pt - center, axis, radians )
            for pt in xyz_cross_sections[-1]
            ]
        
        ## We will now walk through all intersected edges and
        ##   1 constrain the edge vector to a rotated version of itself, where the rotation is the
        ##      minimal rotation such that the dot product with the circle normal is 0.
        
        ## Remove the first and last edges, which connect to boundary vertices.
        ## I don't know if this is necessary, but my constraints are seeming too harsh,
        ## so perhaps this will help.
        ## TODO: Revisit this.
        cs_front_edges = cs_front_edges[1:-1]
        cs_back_edges = cs_back_edges[1:-1]
        
        ## Add xy and relative z constraints on the intersections in front/back_edges.
        for i, edge_t in enumerate( cs_front_edges + cs_back_edges ):
            ## 1
            ## Finally, rotate the edge so that it's on the plane perpendicular to the normal
            ## vector dir(pos-center).
            ## Using the intersection point of the edge as our origin, project the edge
            ## onto the plane passing through the origin with normal dir(pos-center).
            ## Rescale the edge so its length is preserved, and make that the constraint.
            edge0 = mesh.vs[ edge_t[0][0] ]
            edge1 = mesh.vs[ edge_t[0][1] ]
            edge0rot = center + rotate( edge0 - center, axis, radians )
            edge1rot = center + rotate( edge1 - center, axis, radians )
            ## Constrain the edge to this vector.
            xy_vertex_weighted_sum_indices.append( edge_t[0] )
            xy_vertex_weighted_sum_weights.append( [1, -1] )
            xy_vertex_weighted_sum_rhss.append( ( edge0rot - edge1rot )[:2] )
            z_vertex_weighted_sum_indices.append( edge_t[0] )
            z_vertex_weighted_sum_weights.append( [1, -1] )
            z_vertex_weighted_sum_rhss.append( ( edge0rot - edge1rot )[2] )
            
            region._decs.append( (
                'line',
                ( 1, 0, 0 ),
                edge0,
                edge1
                ) )
            
            region._decs.append( (
                'line',
                ( .5, .5, 0 ),
                edge0rot,
                edge1rot
                ) )
    
    return xyz_cross_sections

def add_cross_section_constraints2_obliqueA(
    region,
    mesh, boundary_vertices, front_vertices, back_vertices,
    xy_verts, xy_pos,
    zrel_groups_vert, zrel_groups_vals,
    xy_vertex_weighted_sum_indices, xy_vertex_weighted_sum_weights, xy_vertex_weighted_sum_rhss,
    z_vertex_weighted_sum_indices, z_vertex_weighted_sum_weights, z_vertex_weighted_sum_rhss,
    general_linear_constraints,
    oblique
    ):
    '''
    Implements Step 5 from inflate() in "annotations revisited data structures.txt" in the
    following manner:
    a We'll find a straight line for each cross section by intersecting with the mesh the line
      passing through the anchor point and perpendicular to the line from it to its destination
      location.
    b Compute the sphere and geodesics implied by the anchor point and the anchor point's destination.
    c For intersected mesh edges, constrain the edge to intersect the sphere and have the appropriate normal.
    
    Returns a list of line loops corresponding to the xyz cross sections.
    
    used
    '''
    
    if not hasattr( region, '_decs' ): region._decs = []
    
    
    constrained_edges = set()
    
    xyz_cross_sections = []
    for cross_section in region.cross_sections:
        xyz_cross_sections.append( [] )
        
        ## a
        anchor = asarrayf( cross_section[0] )
        ## Find the direction perpendicular to the spine at the anchor point.
        anchor_to_dest = asarrayf( cross_section[1] ) - anchor
        anchor_to_dest_perp = anchor_to_dest.copy()
        anchor_to_dest_perp[0], anchor_to_dest_perp[1] = anchor_to_dest_perp[1], -anchor_to_dest_perp[0]
        if mag2( anchor_to_dest_perp ) < kEps:
            print("WARNING: No spine and anchor didn't move.  Using <0,1> as the spine_perp.")
            anchor_to_dest_perp[:] = (0,1)
        
        ## TODO Q: What is a cross section under oblique projection?
        ##         Candidates:
        ##             A Use slices of the oblique skewed mesh that project obliquely to the same xy.
        ##               Pros: Similar to / Degenerates to the orthographic case.
        ##               Cons: The slice of the mesh is unintuitive.  The slices are skewed
        ##                     opposite the oblique skew (lines which are in the -skew direction
        ##                     project to points), which is strange for e.g. the front face of a
        ##                     smoothed cuboid with xyz edge lengths <1,1,3>.
        ##             B Use slices of the oblique skewed mesh that project orthographically to the same xy.
        ##               This is the same as option A.
        ##             C Find the point on the mesh under the xy pos (invert the oblique projection).
        ##               Then use the slice of the mesh that intersects the point and has the
        ##               same skew as the oblique projection.
        ##               Pros: Similar to / Degenerates to the orthographic case.  The slice of the
        ##                     mesh will be the same when clicking on a part that would be visible
        ##                     under orthographic projection.
        ##               Cons: The slice of the mesh will be the non-existent when clicking on a
        ##                     part that would not be visible under orthographic projection.
        ##             D Find the point on the mesh under the xy pos (invert the oblique projection).
        ##               Then use the normal of the mesh and the destination point projected
        ##               onto the tangent plane to find a cross section slice of the mesh.
        ##               Pros: Intuitive slice, if the existing 3d structure is intuitive.
        ##               Cons: The slices are not the same slices as in the orthographic case.
        ##                     Q: Can I remove the component of the normals that projects into the
        ##                        viewing plane?
        ##                     A: No.  In the ortho case this will always give me <0,0,1>.
        ##                        In the oblique case this will always give me <a,b,-1>.
        boundary1, cs_front_edges, boundary2, cs_back_edges = get_cross_section_edgesA(
            anchor, anchor_to_dest_perp,
            mesh, boundary_vertices, front_vertices, back_vertices,
            oblique
            )
        expect( len( cs_front_edges ) > 0 )
        expect( len( cs_back_edges ) > 0 )
        
        
        '''
        print 'boundary1:', boundary1
        print 'boundary2:', boundary2
        print 'cs_front_edges:', cs_front_edges
        print 'cs_back_edges:', cs_back_edges
        '''
        
        cross_section_all_edges = set()
        cross_section_all_edges.add( frozenset( boundary1[0] ) )
        cross_section_all_edges.add( frozenset( boundary2[0] ) )
        cross_section_all_edges.update( [ frozenset( edge_t[0] ) for edge_t in cs_front_edges ] )
        cross_section_all_edges.update( [ frozenset( edge_t[0] ) for edge_t in cs_back_edges ] )
        if len( constrained_edges.intersection( cross_section_all_edges ) ) > 0:
            print('cross section intersects another cross section; skipping:', cross_section)
            continue
        
        constrained_edges.update( cross_section_all_edges )
        
        ## b
        ## Before constraining:
        xyz_cross_sections[-1] = [
            lerp( mesh.vs[ edge_t[0][0] ], mesh.vs[ edge_t[0][1] ], edge_t[1] )
            for edge_t in [boundary1] + cs_front_edges + [boundary2] + cs_back_edges
            ]
        
        boundary1pos_before = lerp( mesh.vs[ boundary1[0][0] ], mesh.vs[ boundary1[0][1] ], boundary1[1] )
        boundary2pos_before = lerp( mesh.vs[ boundary2[0][0] ], mesh.vs[ boundary2[0][1] ], boundary2[1] )
        
        oblique_dir = arrayf( ( oblique[0], oblique[1], -1 ) )
        ## These can be shifted by adding s*oblique_dir to each for any value of s.
        t_a_minus_t_b = tangent_circle_between_rays_t_a_minus_t_b(
            boundary1pos_before,
            boundary2pos_before,
            oblique_dir
            )
        ## We'll shift boundary1pos_before and boundary2pos_before by equal amounts.
        boundary1pos = boundary1pos_before + .5 * t_a_minus_t_b * oblique_dir
        boundary2pos = boundary2pos_before - .5 * t_a_minus_t_b * oblique_dir
        center = .5 * ( boundary1pos + boundary2pos )
        radius = .5 * mag( boundary1pos - boundary2pos )
        
        ## We want the geodesic which passes through boundary1pos, boundary2pos, and
        ## the ray through the anchor's destination point.
        ## 1 Find the point on the sphere under (by projection) the anchor's destination point.
        anchor_dest_sphere_pt = closest_point_on_sphere_to_ray(
            ( center, radius ),
            ## Pull the ray's point back in the view direction so that the closest point is the one closer to the viewer.
            ( append( anchor + anchor_to_dest, 0. ) - ( center[2] + radius + 1 )*oblique_dir, oblique_dir ),
            ray_is_line = True
            )[0]
        ## 2 Find the other axis of the geodesic, by taking cross products.  (The first axis
        ##   is boundary1pos - center.)
        axis1 = dir( boundary1pos - center )
        axis2dest = dir( cross( cross( axis1, anchor_dest_sphere_pt - center ), axis1 ) )
        ## The normal of the plane this circle lies on:
        axis12dest_normal = dir( cross( axis1, axis2dest ) )
        
        ## We also want the geodesic which passes through boundary1pos, boundary2pos, and
        ## the ray through the anchor's source point.
        ## Find the point on the sphere under (by projection) the anchor's source point.
        anchor_source_sphere_pt = closest_point_on_sphere_to_ray(
            ( center, radius ),
            ## Pull the ray's point back in the view direction so that the closest point is the one closer to the viewer.
            ( append( anchor, 0. ) - ( center[2] + radius + 1 )*oblique_dir, oblique_dir ),
            ray_is_line = True
            )[0]
        ## Find the other axis of the geodesic, by taking cross products.  (The first axis is the
        ## same as axis1 above.)
        axis2source = dir( cross( cross( axis1, anchor_source_sphere_pt - center ), axis1 ) )
        ## The normal of the plane this circle lies on:
        axis12source_normal = dir( cross( axis1, axis2source ) )
        
        ## After placement on the geodesic:
        xyz_cross_sections[-1] = [
            ## Simple projection doesn't work, especially if the circle's plane is perpendicular to the picked slice of edges.
            center + radius * dir( dot( pt - center, axis1 ) * axis1 + dot( pt - center, axis2source ) * axis2dest )
            for pt in xyz_cross_sections[-1]
            ]
        
        ## We will now walk through all intersected edges and
        ##   1 constrain the edge vector to a rotated version of itself, where the rotation is the
        ##      minimal rotation such that the dot product with the circle normal is 0.
        
        ## Remove the first and last edges, which connect to boundary vertices.
        ## I don't know if this is necessary, but my constraints are seeming too harsh,
        ## so perhaps this will help.
        ## TODO Q: Revisit this.  Is this still necessary?
        ## A1: If I remove the first and last edge, the cylindrical constraints collapse,
        ##     and I need an additional distance constraint 'constrain_edge_to_cylinder_center'.
        #cs_front_edges = cs_front_edges[1:-1]
        #cs_back_edges = cs_back_edges[1:-1]
        
        ## Add xy and relative z constraints on the intersections in front/back_edges.
        constrain_edge_lerpval = False
        constrain_edge_to_cylinder = True
        constrain_edge_to_cylinder_center = False
        for edge_t in ( cs_front_edges + cs_back_edges ):
            edge0 = mesh.vs[ edge_t[0][0] ]
            edge1 = mesh.vs[ edge_t[0][1] ]
            
            region._decs.append( (
                'line',
                ( 1, 0, 0 ),
                edge0,
                edge1
                ) )
            
            if constrain_edge_lerpval:
                ## Constrain the edge so that it intersects the anchor destination geodesic circle.
                ## TODO Q: This doesn't work; it constrains edges to self-intersect.  How to fix it?
                pt = lerp( edge0, edge1, edge_t[1] )
                pt_new = center + radius * dir( dot( pt - center, axis1 ) * axis1 + dot( pt - center, axis2source ) * axis2dest )
                
                xy_vertex_weighted_sum_indices.append( edge_t[0] )
                xy_vertex_weighted_sum_weights.append( ( edge_t[1], 1. - edge_t[1] ) )
                xy_vertex_weighted_sum_rhss.append( pt_new[:2] )
                ## TODO Q: Shouldn't these be relative?
                #z_vertex_weighted_sum_indices.append( edge_t[0] )
                #z_vertex_weighted_sum_weights.append( ( edge_t[1], 1. - edge_t[1] ) )
                #z_vertex_weighted_sum_rhss.append( pt_new[2] )
            
            if constrain_edge_to_cylinder:
                ## Constrain the edge so that it intersects the cylinder through the anchor
                ## destination geodesic circle.
                
                edge0proj = (
                    center +
                    radius * dir(
                        dot( edge0 - center, axis1 ) * axis1 + dot( edge0 - center, axis2source ) * axis2dest
                        ) +
                    dot( edge0 - center, axis12source_normal ) * axis12dest_normal
                    )
                edge1proj = (
                    center +
                    radius * dir(
                        dot( edge1 - center, axis1 ) * axis1 + dot( edge1 - center, axis2source ) * axis2dest
                        ) +
                    dot( edge1 - center, axis12source_normal ) * axis12dest_normal
                    )
                
                ## Constrain the edge to this vector.
                xy_vertex_weighted_sum_indices.append( edge_t[0] )
                xy_vertex_weighted_sum_weights.append( [1, -1] )
                xy_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[:2] )
                z_vertex_weighted_sum_indices.append( edge_t[0] )
                z_vertex_weighted_sum_weights.append( [1, -1] )
                z_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[2] )
                
                region._decs.append( (
                    'line',
                    ( .5, .5, 0 ),
                    edge0proj,
                    edge1proj
                    ) )
            
            if constrain_edge_to_cylinder_center:
                ## Also constrain the edge to lie radius away from the cylinder's axis.
                ## TODO Q: Should I do this only for the two edges in the middle of cs_front/back_edges?
                if edge_t not in ( cs_front_edges[ len(cs_front_edges)//2 ], cs_back_edges[ len(cs_back_edges)//2 ] ): continue
                #print "constrain_edge_to_cylinder_center operating on edge", edge_t
                
                edge01 = lerp( edge0, edge1, edge_t[1] )
                edge01proj = center + radius * dir( dot( edge01 - center, axis1 ) * axis1 + dot( edge01 - center, axis2source ) * axis2dest )
                
                ## This expression corresponds to
                ## lerp( edge_t ) - .5 * ( lerp( boundary1 ) + lerp( boundary2 ) )
                expr = (
                    edge_t[1], 1. - edge_t[1],
                    -.5 * boundary1[1], -.5 * ( 1. - boundary1[1] ),
                    -.5 * boundary2[1], -.5 * ( 1. - boundary2[1] )
                    )
                verts = list( edge_t[0] ) + list( boundary1[0] ) + list( boundary2[0] )
                rhs = edge01proj - center
                
                xy_vertex_weighted_sum_indices.append( verts )
                xy_vertex_weighted_sum_weights.append( expr )
                xy_vertex_weighted_sum_rhss.append( rhs[:2] )
                z_vertex_weighted_sum_indices.append( verts )
                z_vertex_weighted_sum_weights.append( expr )
                z_vertex_weighted_sum_rhss.append( rhs[2] )
    
    return xyz_cross_sections

def gen_normal_constraint( mesh, n_xyz, v2, v1, v2_to_v1 ):
    '''
    Implements 'step d' in add_cross_section_constraints().
    Returns a list of general_linear_constraints to add to the fibermesh.Constraints object.
    '''
    
    ## For a description of this function, see "note on normal/silhouette constraints" at the bottom
    ## of this file.
    
    print('gen_normal_constraint:', v2, v1)
    assert v2 != v1
    
    linear_constraints = []
    
    ## Gather vertices involved.
    v2v1he = mesh.halfedges[ mesh.directed_edge2he_index( ( v2, v1 ) ) ]
    v4 = mesh.halfedges[ v2v1he.next_he ].to_vertex
    v3 = mesh.halfedges[ mesh.halfedges[ v2v1he.opposite_he ].next_he ].to_vertex
    assert v3 != v2
    assert v4 != v2
    assert v4 != v3
    #if len( set( [ v1, v2, v3, v4 ] ) ) != len( [ v1, v2, v3, v4 ] ): return []
    #v2v1_dir = dir( mesh.vs[ v1 ] - mesh.vs[ v2 ] )
    v2v1_dir = dir( v2_to_v1 )
    s4 = 1. / mag( cross( v2v1_dir, mesh.vs[ v4 ] - mesh.vs[ v2 ] ) )
    s3 = 1. / mag( cross( v2v1_dir, mesh.vs[ v2 ] - mesh.vs[ v3 ] ) )
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [ (v4,2), (v2,2), (v4,1), (v2,1), (v2,2), (v3,2), (v2,1), (v3,1) ],
        ## The list of coefficients for each term.
        [
            s4 * v2v1_dir[1], -s4 * v2v1_dir[1], -s4 * v2v1_dir[2], s4 * v2v1_dir[2],
            s3 * v2v1_dir[1], -s3 * v2v1_dir[1], -s3 * v2v1_dir[2], s3 * v2v1_dir[2]
        ],
        ## The right hand side value.
        n_xyz[0]
        ) )
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [ (v4,0), (v2,0), (v4,2), (v2,2), (v2,0), (v3,0), (v2,2), (v3,2) ],
        ## The list of coefficients for each term.
        [
            s4 * v2v1_dir[2], -s4 * v2v1_dir[2], -s4 * v2v1_dir[0], s4 * v2v1_dir[0],
            s3 * v2v1_dir[2], -s3 * v2v1_dir[2], -s3 * v2v1_dir[0], s3 * v2v1_dir[0]
        ],
        ## The right hand side value.
        n_xyz[1]
        ) )
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [ (v4,1), (v2,1), (v4,0), (v2,0), (v2,1), (v3,1), (v2,0), (v3,0) ],
        ## The list of coefficients for each term.
        [
            s4 * v2v1_dir[0], -s4 * v2v1_dir[0], -s4 * v2v1_dir[1], s4 * v2v1_dir[1],
            s3 * v2v1_dir[0], -s3 * v2v1_dir[0], -s3 * v2v1_dir[1], s3 * v2v1_dir[1]
        ],
        ## The right hand side value.
        n_xyz[2]
        ) )
    
    ## Normalize the whole thing.
    current_normal = asarray( [
        sum( [
            mesh.vs[ vi, coord ] * coeff
            for ((vi, coord), coeff) in zip( linear_constraints[c][0], linear_constraints[c][1] )
            ] )
        for c in (0,1,2)
        ] )
    current_normal_mag = mag( current_normal )
    linear_constraints = [
        ( verts, coeffs, rhs * current_normal_mag )
        for verts, coeffs, rhs in linear_constraints
        ]
    
    return linear_constraints

def gen_normal_constraint2( mesh, n_xyz, v2, v1, v2_to_v1 ):
    '''
    Implements 'step d' in add_cross_section_constraints().
    Returns a list of general_linear_constraints to add to the fibermesh.Constraints object.
    '''
    
    ## For a description of this function, see "note on normal/silhouette constraints" at the bottom
    ## of this file.
    
    print('gen_normal_constraint:', v2, v1)
    assert v2 != v1
    
    linear_constraints = []
    
    ## Gather vertices involved.
    v2v1he = mesh.halfedges[ mesh.directed_edge2he_index( ( v2, v1 ) ) ]
    v4 = mesh.halfedges[ v2v1he.next_he ].to_vertex
    v3 = mesh.halfedges[ mesh.halfedges[ v2v1he.opposite_he ].next_he ].to_vertex
    assert v3 != v2
    assert v4 != v2
    assert v4 != v3
    #if len( set( [ v1, v2, v3, v4 ] ) ) != len( [ v1, v2, v3, v4 ] ): return []
    #v2v1_dir = dir( mesh.vs[ v1 ] - mesh.vs[ v2 ] )
    v2v1_dir = dir( v2_to_v1 )
    
    def nx_gen( mesh ):
        s4 = 1. / mag( cross( v2v1_dir, mesh.vs[ v4 ] - mesh.vs[ v2 ] ) )
        s3 = 1. / mag( cross( v2v1_dir, mesh.vs[ v2 ] - mesh.vs[ v3 ] ) )
        return [
            s4 * v2v1_dir[1], -s4 * v2v1_dir[1], -s4 * v2v1_dir[2], s4 * v2v1_dir[2],
            s3 * v2v1_dir[1], -s3 * v2v1_dir[1], -s3 * v2v1_dir[2], s3 * v2v1_dir[2]
        ]
    def ny_gen( mesh ):
        s4 = 1. / mag( cross( v2v1_dir, mesh.vs[ v4 ] - mesh.vs[ v2 ] ) )
        s3 = 1. / mag( cross( v2v1_dir, mesh.vs[ v2 ] - mesh.vs[ v3 ] ) )
        return [
            s4 * v2v1_dir[2], -s4 * v2v1_dir[2], -s4 * v2v1_dir[0], s4 * v2v1_dir[0],
            s3 * v2v1_dir[2], -s3 * v2v1_dir[2], -s3 * v2v1_dir[0], s3 * v2v1_dir[0]
        ]
    def nz_gen( mesh ):
        s4 = 1. / mag( cross( v2v1_dir, mesh.vs[ v4 ] - mesh.vs[ v2 ] ) )
        s3 = 1. / mag( cross( v2v1_dir, mesh.vs[ v2 ] - mesh.vs[ v3 ] ) )
        return [
            s4 * v2v1_dir[0], -s4 * v2v1_dir[0], -s4 * v2v1_dir[1], s4 * v2v1_dir[1],
            s3 * v2v1_dir[0], -s3 * v2v1_dir[0], -s3 * v2v1_dir[1], s3 * v2v1_dir[1]
        ]
    
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [ (v4,2), (v2,2), (v4,1), (v2,1), (v2,2), (v3,2), (v2,1), (v3,1) ],
        ## The list of coefficients for each term.
        nx_gen,
        ## The right hand side value.
        n_xyz[0]
        ) )
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [ (v4,0), (v2,0), (v4,2), (v2,2), (v2,0), (v3,0), (v2,2), (v3,2) ],
        ## The list of coefficients for each term.
        ny_gen,
        ## The right hand side value.
        n_xyz[1]
        ) )
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [ (v4,1), (v2,1), (v4,0), (v2,0), (v2,1), (v3,1), (v2,0), (v3,0) ],
        ## The list of coefficients for each term.
        nz_gen,
        ## The right hand side value.
        n_xyz[2]
        ) )
    
    print("TODO: normalize the whole thing")
    
    return linear_constraints

def get_cross_section_vertices( anchor, spine_perp, mesh, boundary_vertices, front_vertices, back_vertices ):
    '''
    Given a 2d point 'anchor', a 2d vector perpendicular to the spine at 'anchor',
    a TriMesh 'mesh', and a partition of the vertices in 'mesh' into three lists containing
    indices of its vertices along the silhouette ('boundary_vertices') in front ('front_vertices'),
    and in back ('back_vertices'),
    returns a tuple of four items (the names here are suggestions):
        boundary1: the index of the vertex in mesh corresponding to the first boundary vertex
          in the direction 'spine_perp' from 'anchor',
        front_vertices: the list of mesh vertices along the path from 'boundary1' to 'boundary2'
          on the front of the mesh,
        boundary2: the index of the vertex in mesh corresponding to the first boundary vertex
          in the direction '-spine_perp' from 'anchor',
        back_vertices: the list of mesh vertices along the path from 'boundary2' to 'boundary1'
          on the back of the mesh
    
    tested
    '''
    
    #debugger()
    
    assert len( anchor ) == 2
    assert len( spine_perp ) == 2
    anchor = asarray( anchor )
    spine_perp = asarrayf( spine_perp )
    assert mag2( spine_perp ) > kEps
    ## Normalize spine_perp
    spine_perp = dir( spine_perp )
    
    ## In reality, I don't want to move too many mesh vertices, so snap to the boundary points.
    ## Specifically, shoot the rays ( anchor, +/- spine_perp ) and find the closest two boundary
    ## vertices.  Then set anchor to the midpoint of them, and update spine_perp accordingly.
    
    ## Let's cast the ray from anchor.
    ## TODO: Factor this into intersect_ray_with_mesh_edges_2d( ray, mesh, candidate_vertices ).
    ## First find the bounding sphere diameter.
    ## We'll use the diagonal of the axis-aligned bounding rectangle.
    mesh.vs = asarray( mesh.vs )
    maxxy = mesh.vs[:,:2].max(0)
    minxy = mesh.vs[:,:2].min(0)
    bounding_sphere_diameter = sqrt( ( ( maxxy - minxy )**2 ).sum() )
    
    eps = kEps
    ## Find boundary1, the boundary in the +spine_perp direction.
    line_segment = ( anchor, anchor + bounding_sphere_diameter * spine_perp )
    boundary1 = []
    for edge in mesh.edges:
        if edge[0] in boundary_vertices and edge[1] in boundary_vertices:
            t1, t2 = line_line_intersection( line_segment, ( mesh.vs[ edge[0] ], mesh.vs[ edge[1] ] ) )
            if(
                t1 > -eps and
                t1 < 1 + eps and
                t2 > -eps and
                t2 < 1 + eps
                ):
                boundary1.append( ( t1, t2, edge ) )
    
    assert len( boundary1 ) > 0
    boundary1.sort()
    ## Take the closest edge.
    boundary1 = boundary1[0]
    ## Take the closest vertex.
    boundary1 = boundary1[2][ 0 if boundary1[1] < .5 else 1 ]
    
    ## Find boundary2, the boundary in the -spine_perp direction.
    line_segment = ( anchor, anchor - bounding_sphere_diameter * spine_perp )
    boundary2 = []
    for edge in mesh.edges:
        if edge[0] in boundary_vertices and edge[1] in boundary_vertices:
            t1, t2 = line_line_intersection( line_segment, ( mesh.vs[ edge[0] ], mesh.vs[ edge[1] ] ) )
            if(
                t1 > -eps and
                t1 < 1 + eps and
                t2 > -eps and
                t2 < 1 + eps
                ):
                boundary2.append( ( t1, t2, edge ) )
    
    assert len( boundary2 ) > 0
    boundary2.sort()
    ## Take the closest edge.
    boundary2 = boundary2[0]
    ## Take the closest vertex.
    boundary2 = boundary2[2][ 0 if boundary2[1] < .5 else 1 ]
    
    ## Update anchor.
    anchor = .5 * ( mesh.vs[ boundary1 ] + mesh.vs[ boundary2 ] )
    ## We don't need spine_perp anymore, and it's out of date.
    del spine_perp
    #spine_perp = roll( mesh.vs[ boundary2 ] - mesh.vs[ boundary1 ], 1 ) * array((1,-1))
    
    metric = lambda u, v: mag( mesh.vs[u][:2] - mesh.vs[v][:2] )
    front_vertices = shortest_path_on_mesh( mesh, boundary1, boundary2, restrict_vertices = [boundary1] + list(front_vertices) + [boundary2], metric = metric )
    back_vertices = shortest_path_on_mesh( mesh, boundary2, boundary1, restrict_vertices = [boundary1] + list(back_vertices) + [boundary2], metric = metric )
    ## We're returning front and back vertices separately from the two boundary (path start and end)
    ## vertices, so remove the start and end from the paths.
    assert front_vertices[0] == boundary1
    assert front_vertices[-1] == boundary2
    assert back_vertices[0] == boundary2
    assert back_vertices[-1] == boundary1
    front_vertices = front_vertices[1:-1]
    back_vertices = back_vertices[1:-1]
    
    return boundary1, front_vertices, boundary2, back_vertices


def get_cross_section_edges( anchor, spine_perp, mesh, boundary_vertices, front_vertices, back_vertices ):
    '''
    Given a 2d point 'anchor', a 2d vector perpendicular to the spine at 'anchor',
    a TriMesh 'mesh', and a partition of the vertices in 'mesh' into three lists containing
    indices of its vertices along the silhouette ('boundary_vertices') in front ('front_vertices'),
    and in back ('back_vertices'),
    returns a tuple of four items (the names here are suggestions):
        boundary1: an edge (pair of vertices) and a t value such that
          lerp( mesh.vs[ edge[0] ], mesh.vs[ edge[1] ], t ) corresponds to the first
          boundary edge intersected in the direction 'spine_perp' from 'anchor',
        front_edges: the ordered list of intersected edges and t's (the same format as boundary1)
          along the line from 'boundary1' to 'boundary2' on the front of the mesh,
        boundary2: the first boundary edge intersected in the direction 'spine_perp' from 'anchor'
          (in the same format as boundary1),
        back_edges: the ordered list of intersected edges and t's (the same format as boundary1)
          along the line from 'boundary2' to 'boundary1' on the back of the mesh
    NOTE: The first and last edge in both 'front_edges' and 'back_edges' will have one vertex
          in 'boundary_vertices' and one in 'front/back_vertices'.
    
    used
    '''
    
    #debugger()
    
    assert len( anchor ) == 2
    assert len( spine_perp ) == 2
    anchor = asarray( anchor )
    spine_perp = asarrayf( spine_perp )
    assert mag2( spine_perp ) > kEps
    ## Normalize spine_perp
    spine_perp = dir( spine_perp )
    
    ## In reality, I don't want to move too many mesh vertices, so snap to the boundary points.
    ## Specifically, shoot the rays ( anchor, +/- spine_perp ) and find the closest two boundary
    ## vertices.  Then set anchor to the midpoint of them, and update spine_perp accordingly.
    
    ## Cast the ray from anchor in the +spine_perp direction.
    eps = kEps
    line_segment = ( anchor, anchor + spine_perp )
    boundary_intersections = []
    front_intersections = []
    back_intersections = []
    for edge in mesh.edges:
        t1, t2 = line_line_intersection( line_segment, ( mesh.vs[ edge[0] ], mesh.vs[ edge[1] ] ) )
        ## We only care if the intersection is within the mesh edge line segment.
        ## t1 only tells us if the edge is in the positive or negative 'spine_perp' direction.
        if(
            t2 > -eps and
            t2 < 1 + eps
            ):
            if edge[0] in boundary_vertices and edge[1] in boundary_vertices:
                boundary_intersections.append( ( t1, t2, edge ) )
            elif edge[0] in front_vertices or edge[1] in front_vertices:
                assert edge[0] in front_vertices or edge[0] in boundary_vertices
                assert edge[1] in front_vertices or edge[1] in boundary_vertices
                front_intersections.append( ( t1, t2, edge ) )
            elif edge[0] in back_vertices or edge[1] in back_vertices:
                assert edge[0] in back_vertices or edge[0] in boundary_vertices
                assert edge[1] in back_vertices or edge[1] in boundary_vertices
                back_intersections.append( ( t1, t2, edge ) )
            else:
                raise RuntimeError("this can't happen!")
    
    assert len( boundary_intersections ) > 0
    assert len( front_intersections ) > 0
    assert len( back_intersections ) > 0
    
    ## Find the first boundary intersections in the positive and negative t1 directions.
    boundary1 = min( [ isect for isect in boundary_intersections if isect[0] > 0. ] )
    boundary2 = max( [ isect for isect in boundary_intersections if isect[0] < 0. ] )
    
    ## Find all front intersections between the boundary intersections and order them.
    front_edges = [
        isect
        for isect in front_intersections
        if isect[0] > boundary2[0] and isect[0] < boundary1[0]
        ]
    front_edges.sort()
    ## Reverse front_edges, because it leads from boundary1 to boundary2.
    front_edges.reverse()
    
    ## Find all back intersections between the boundary intersections and order them.
    back_edges = [
        isect
        for isect in back_intersections
        if isect[0] > boundary2[0] and isect[0] < boundary1[0]
        ]
    back_edges.sort()
    
    ## Finally, swizzle the return values into the format we said we would return:
    ## ( edge vertex 1, edge vertex 2 ), t along the edge
    boundary1 = ( boundary1[2], boundary1[1] )
    boundary2 = ( boundary2[2], boundary2[1] )
    front_edges = [ ( isect[2], isect[1] ) for isect in front_edges ]
    back_edges = [ ( isect[2], isect[1] ) for isect in back_edges ]
    
    return boundary1, front_edges, boundary2, back_edges


def get_cross_section_edgesA( anchor, spine_perp, mesh, boundary_vertices, front_vertices, back_vertices, oblique ):
    '''
    Given a 2d point 'anchor', a 2d vector perpendicular to the spine at 'anchor',
    a TriMesh 'mesh', and a partition of the vertices in 'mesh' into three lists containing
    indices of its vertices along the silhouette ('boundary_vertices') in front ('front_vertices'),
    and in back ('back_vertices'),
    returns a tuple of four items (the names here are suggestions):
        boundary1: an edge (pair of vertices) and a t value such that
          lerp( mesh.vs[ edge[0] ], mesh.vs[ edge[1] ], t ) corresponds to the first
          boundary edge intersected in the direction 'spine_perp' from 'anchor',
        front_edges: the ordered list of intersected edges and t's (the same format as boundary1)
          along the line from 'boundary1' to 'boundary2' on the front of the mesh,
        boundary2: the first boundary edge intersected in the direction 'spine_perp' from 'anchor'
          (in the same format as boundary1),
        back_edges: the ordered list of intersected edges and t's (the same format as boundary1)
          along the line from 'boundary2' to 'boundary1' on the back of the mesh
    NOTE: The first and last edge in both 'front_edges' and 'back_edges' will have one vertex
          in 'boundary_vertices' and one in 'front/back_vertices'.
    
    ## This function is so named because it implements option D in add_cross_section_constraints2_obliqueD():
    ##             A Use slices of the oblique skewed mesh that project obliquely to the same xy.
    ##               Pros: Similar to / Degenerates to the orthographic case.
    ##               Cons: The slice of the mesh is unintuitive.  The slices are skewed
    ##                     opposite the oblique skew (lines which are in the -skew direction
    ##                     project to points), which is strange for e.g. the front face of a
    ##                     smoothed cuboid with xyz edge lengths <1,1,3>.
    
    untested
    '''
    
    #debugger()
    
    assert len( anchor ) == 2
    assert len( spine_perp ) == 2
    anchor = asarray( anchor )
    spine_perp = asarrayf( spine_perp )
    assert mag2( spine_perp ) > kEps
    ## Normalize spine_perp
    spine_perp = dir( spine_perp )
    
    assert len( oblique ) == 2
    
    
    ## 1 Construct the skewed mesh.
    ## 2 Find cross section edges the usual way.
    
    
    ## 1
    ## Since this is a linear transform, any lerp values we compute on
    ## the skewed mesh are also correct for the unskewed mesh.
    from mesh_gui import TriMesh
    mesh_skew = TriMesh()
    mesh_skew.faces = array( mesh.faces )
    skew = identity(3)
    skew[0,2] = oblique[0]
    skew[1,2] = oblique[1]
    mesh_skew.vs = dot( skew, asarray( mesh.vs ).T ).T
    mesh_skew.topology_changed()
    
    ## 2
    return get_cross_section_edges( anchor, spine_perp, mesh_skew, boundary_vertices, front_vertices, back_vertices )


def get_cross_section_edgesD( anchor, spine_perp, mesh, boundary_vertices, front_vertices, back_vertices, oblique, extra = {} ):
    '''
    Given a 2d point 'anchor', a 2d vector perpendicular to the spine at 'anchor',
    a TriMesh 'mesh', and a partition of the vertices in 'mesh' into three lists containing
    indices of its vertices along the silhouette ('boundary_vertices') in front ('front_vertices'),
    and in back ('back_vertices'),
    returns a tuple of four items (the names here are suggestions):
        boundary1: an edge (pair of vertices) and a t value such that
          lerp( mesh.vs[ edge[0] ], mesh.vs[ edge[1] ], t ) corresponds to the first
          boundary edge intersected in the direction 'spine_perp' from 'anchor',
        front_edges: the ordered list of intersected edges and t's (the same format as boundary1)
          along the line from 'boundary1' to 'boundary2' on the front of the mesh,
        boundary2: the first boundary edge intersected in the direction 'spine_perp' from 'anchor'
          (in the same format as boundary1),
        back_edges: the ordered list of intersected edges and t's (the same format as boundary1)
          along the line from 'boundary2' to 'boundary1' on the back of the mesh
    NOTE: The first and last edge in both 'front_edges' and 'back_edges' will have one vertex
          in 'boundary_vertices' and one in 'front/back_vertices'.
    
    ## This function is so named because it implements option D in add_cross_section_constraints2_obliqueD():
    ##             D Find the point on the mesh under the xy pos (invert the oblique projection).
    ##               Then use the normal of the mesh and the destination point projected
    ##               onto the tangent plane to find a cross section slice of the mesh.
    ##               Pros: Intuitive slice, if the existing 3d structure is intuitive.
    ##               Cons: The slices are not the same slices as in the orthographic case.
    ##                     Q: Can I remove the component of the normals that projects into the
    ##                        viewing plane?
    ##                     A: No.  In the ortho case this will always give me <0,0,1>.
    ##                        In the oblique case this will always give me <a,b,-1>.
    
    tested
    '''
    
    #debugger()
    
    assert len( anchor ) == 2
    assert len( spine_perp ) == 2
    anchor = asarray( anchor )
    spine_perp = asarrayf( spine_perp )
    assert mag2( spine_perp ) > kEps
    ## Normalize spine_perp
    spine_perp = dir( spine_perp )
    
    assert len( oblique ) == 2
    
    
    ## 1 Construct the skewed mesh.
    ## 2 Intersect anchor with the skewed mesh to get the intersection point and face 'isect_pt' and 'isect_face'.
    ## 3 Find 'isect_normal', the unskewed mesh's normal at 'isect_pt'.
    ## 4 Project 'spine_perp' onto the tangent plane at 'isect_pt' to get 'spine_perp_proj'.
    ## 5 Rotate the un-skewed mesh and 'spine_perp_proj' so that 'isect_normal' points in the <0,0,1> direction.
    ## 6 Find cross section edges the usual way, using the rotated 'isect_pt' and 'spine_perp_proj'.
    
    
    ## 1
    ## Since this is a linear transform, any lerp values we compute on
    ## the skewed mesh are also correct for the unskewed mesh.
    from mesh_gui import TriMesh
    mesh_skew = TriMesh()
    mesh_skew.faces = array( mesh.faces )
    skew = identity(3)
    skew[0,2] = oblique[0]
    skew[1,2] = oblique[1]
    mesh_skew.vs = dot( skew, asarray( mesh.vs ).T ).T
    mesh_skew.topology_changed()
    
    ## 2
    isects = ray_mesh_intersections( ( append( anchor, 0. ), (0,0,-1) ), mesh_skew )
    expect( len( isects ) > 0 )
    if len( isects ) == 0: return []
    isect = min( isects )
    isect_pt = append( anchor, 0. ) + isect[0] * array((0.,0.,-1.))
    isect_face = isect[1]
    
    ## 3
    isect_normal = mesh.face_normals[ isect_face ]
    
    ## 4
    ## I'll do this by intersecting the oblique projection ray with the tangent plane.
    ## (I could also do this by projecting onto the skewed tangent plane and then unskewing.)
    ray_pt = append( array( anchor ) + array( spine_perp ), 0 )
    assert len( ray_pt ) == 3
    ## Shoot the ray in the negative screen-space z direction.
    ## (Invert the (oblique) projection matrix to get the actual ray.)
    ray_dir = ( oblique[0], oblique[1], -1 )
    spine_perp_proj = ray_plane_intersection( (ray_pt, ray_dir), (isect_pt, isect_normal), True )
    assert spine_perp_proj is not None
    spine_perp_proj = spine_perp_proj - isect_pt
    
    ## 5
    ## Re-use (rename) the mesh_skew structure.
    mesh_rot = mesh_skew
    del mesh_skew
    rot_mat = asarray( get_rotation_matrix( *get_axis_angle_rotation( isect_normal, (0,0,1) ) ) )
    mesh_rot.vs = dot( rot_mat, asarray( mesh.vs ).T ).T
    mesh_rot.positions_changed()
    
    isect_pt_rot = dot( rot_mat, isect_pt )
    spine_perp_proj_rot = dot( rot_mat, spine_perp_proj )
    
    extra['tangent_plane'] = ( isect_pt, isect_normal.copy() )
    extra['spine_perp_proj'] = spine_perp_proj
    
    return get_cross_section_edges( isect_pt_rot[:2], spine_perp_proj_rot[:2], mesh_rot, boundary_vertices, front_vertices, back_vertices )


def get_rotation_about_boundary1to2axis_for_cross_section( boundary1, boundary2, anchor_to_dest ):
    '''
    Given a cross section diameter defined as the line between 3d points 'boundary1' and 'boundary2',
    returns the rotation radians necessary to offset the point directly between them (and lifted
    in z) by 2d point 'anchor_to_dest' given the rotation axis (boundary2 - boundary1).
    If 'anchor_to_dest' is not orthogonal to (boundary2 - boundary1), it is projected onto that
    plane.
    
    tested
    '''
    
    #debugger()
    
    assert len( boundary1 ) == 3
    assert len( boundary2 ) == 3
    assert len( anchor_to_dest ) == 2
    
    boundary1 = asarrayf( boundary1 )
    boundary2 = asarrayf( boundary2 )
    anchor_to_dest = append( asarrayf( anchor_to_dest ), 0. )
    
    b1to2center = .5 * ( boundary1 + boundary2 )
    b1to2 = boundary2 - boundary1
    cross_section_radius = .5 * mag( b1to2 )
    
    ## Project 'anchor_to_dest' onto the plane perpendicular to the axis and passing through
    ## the center.  We project onto this plane because rotations around the axis must stay
    ## on this plane.  (They must also lie in a circle on this plane.)
    anchor_to_dest = project_point_to_plane(
        b1to2center + anchor_to_dest,
        ( b1to2center, b1to2 )
        ) - b1to2center
    
    ## Now find the right lift in z so that the cross section passes through the
    ## anchor_dest point in xy.
    ## We want the radians of rotation about the axis b1to2.  All the action is taking
    ## place on a plane perpendicular to the axis.  Consider the 2d graph of this plane where the
    ## vector we're rotating points in the y direction and the vector orthogonal to the axis
    ## points along the x axis.  Let x be the magnitude of the destination vector (from the origin).
    ## Then pi/2 - acos( x / radius ) is the magnitude of rotation about axis.
    ## We can determine the correct sign by taking the cross product of the start and end rotated
    ## vectors and testing if it points in the same direction as the axis (dot product > 0).
    ## Clip x to 1 in case the user drags the cross section too far.
    radians = pi/2 - acos( min( mag( anchor_to_dest[:2] ) / cross_section_radius, 1. ) )
    if dot( cross( b1to2, array( (0,0,1) ) ), anchor_to_dest ) < 0.:
        radians = -radians
    
    return radians

def add_silhouette_constraints(
    region,
    mesh, boundary_vertices, front_vertices, back_vertices,
    boundary_vi_to_curve_vi,
    linear_constraints,
    oblique = None
    ):
    assert oblique is None or len( oblique ) == 2
    
    ## The mesh's boundary_vertices and region.region_boundary aren't the same thing; at best,
    ## the mesh's boundary_vertices are subsampled from region.region_boundary.
    ## While in general it may be possible for the mesh's boundary vertices to have edge-less gaps,
    ## where two vertices are opposite a perpendicular edge between them, the way our mesh is
    ## constructed this is impossible.
    ## So, what I need to find is the set of edges in boundary_vertices with both endpoints'
    ## corresponding curve_vi's in the set of vertices in region.silhouette_edges.
    ## These are the edges I want to add silhouette constraints to.
    
    region_boundary_silhouette_vertices = set()
    region_boundary_silhouette_vertices.update( region.silhouette_edges )
    region_boundary_silhouette_vertices.update( [ ( ei + 1 ) % len( region.region_boundary ) for ei in region.silhouette_edges ] )
    
    silhouette_boundary_edges = []
    for i in range(len( boundary_vertices )):
        if boundary_vi_to_curve_vi[i] not in region_boundary_silhouette_vertices: continue
        
        vi = boundary_vertices[i]
        vi_neighs = mesh.vertex_vertex_neighbors( vi )
        for j in range( i+1, len( boundary_vertices ) ):
            if boundary_vi_to_curve_vi[j] not in region_boundary_silhouette_vertices: continue
            
            vj = boundary_vertices[j]
            if vj in vi_neighs:
                silhouette_boundary_edges.append( frozenset( ( i, j ) ) )
    assert len( silhouette_boundary_edges ) == len( set( silhouette_boundary_edges ) )
    
    for i,j in silhouette_boundary_edges:
        v1, v2 = boundary_vertices[i], boundary_vertices[j]
        v1pos = asarray( region.region_boundary[ boundary_vi_to_curve_vi[i] ] )
        v2pos = asarray( region.region_boundary[ boundary_vi_to_curve_vi[j] ] )
        v2_to_v1 = v1pos - v2pos
        if oblique is None:
            #linear_constraints.extend( gen_silhouette_constraint2( mesh, v2, v1, v2_to_v1 ) )
            linear_constraints.extend( kUseGenSilhouetteConstraint( mesh, v2, v1, v2_to_v1 ) )
        else:
            linear_constraints.extend( gen_silhouette_constraint_oblique( mesh, v2, v1, v2_to_v1, oblique ) )

def gen_silhouette_constraint( mesh, v2, v1, v2_to_v1 ):
    '''
    Implements 'Step 3' in mesh_llconstraints_from_constraints2().
    Returns a list of general_linear_constraints to add to the fibermesh.Constraints object.
    '''
    
    ## For a description of this function, see "note on normal/silhouette constraints" at the bottom
    ## of this file.
    
    print('gen_silhouette_constraint:', v2, v1)
    assert v2 != v1
    
    linear_constraints = []
    
    ## Gather vertices involved.
    v2v1he = mesh.halfedges[ mesh.directed_edge2he_index( ( v2, v1 ) ) ]
    v4 = mesh.halfedges[ v2v1he.next_he ].to_vertex
    v3 = mesh.halfedges[ mesh.halfedges[ v2v1he.opposite_he ].next_he ].to_vertex
    assert v3 != v2
    assert v4 != v2
    assert v4 != v3
    #if len( set( [ v1, v2, v3, v4 ] ) ) != len( [ v1, v2, v3, v4 ] ): return []
    #v2v1_dir = dir( mesh.vs[ v1 ] - mesh.vs[ v2 ] )
    v2v1_dir = dir( v2_to_v1 )
    #s4 = 1. / mag( cross( v2v1_dir, mesh.vs[ v4 ] - mesh.vs[ v2 ] ) )
    #s3 = 1. / mag( cross( v2v1_dir, mesh.vs[ v2 ] - mesh.vs[ v3 ] ) )
    s4 = 1.
    s3 = 1.
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [ (v4,1), (v2,1), (v4,0), (v2,0), (v2,1), (v3,1), (v2,0), (v3,0) ],
        ## The list of coefficients for each term.
        [
            s4 * v2v1_dir[0], -s4 * v2v1_dir[0], -s4 * v2v1_dir[1], s4 * v2v1_dir[1],
            s3 * v2v1_dir[0], -s3 * v2v1_dir[0], -s3 * v2v1_dir[1], s3 * v2v1_dir[1]
        ],
        ## The right hand side value.
        0.
        ) )
    
    return linear_constraints

def gen_silhouette_constraint2( mesh, v2, v1, v2_to_v1 ):
    '''
    Implements 'Step 3' in mesh_llconstraints_from_constraints2().
    Returns a list of general_linear_constraints to add to the fibermesh.Constraints object.
    '''
    
    ## This silhouette constraint is simple: the two opposite vertices of the face flaps of
    ## a silhouette edge must have equal x and y components.
    ## NOTE: This won't work if there is a tilted cross section.
    
    print('gen_silhouette_constraint2:', v2, v1)
    assert v2 != v1
    
    linear_constraints = []
    
    ## Gather vertices involved.
    v2v1he = mesh.halfedges[ mesh.directed_edge2he_index( ( v2, v1 ) ) ]
    v4 = mesh.halfedges[ v2v1he.next_he ].to_vertex
    v3 = mesh.halfedges[ mesh.halfedges[ v2v1he.opposite_he ].next_he ].to_vertex
    assert v3 != v2
    assert v4 != v2
    assert v4 != v3
    ## v3.x == v4.x
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [ (v4,0), (v3,0) ],
        ## The list of coefficients for each term.
        [ 1., -1. ],
        ## The right hand side value.
        0.
        ) )
    ## v3.y == v4.y
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [ (v4,1), (v3,1) ],
        ## The list of coefficients for each term.
        [ 1., -1. ],
        ## The right hand side value.
        0.
        ) )
    ## I don't need for this to be true.  I added it to help debug gen_silhouette_constraint().
    '''
    ## v3.z == -v4.z
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [ (v4,2), (v3,2) ],
        ## The list of coefficients for each term.
        [ 1., 1. ],
        ## The right hand side value.
        0.
        ) )
    '''
    
    return linear_constraints
kUseGenSilhouetteConstraint = gen_silhouette_constraint2

def gen_silhouette_constraint_oblique( mesh, v2, v1, v2_to_v1skew, oblique ):
    '''
    Implements 'Step 3' in mesh_llconstraints_from_constraints2().
    Returns a list of general_linear_constraints to add to the fibermesh.Constraints object.
    
    'oblique' is a sequence of 2 real numbers, a and b, which skew x and y into x + a*z and y + b*z.
    'v2_to_v1skew' is a sequence of 2 real numbers, the evaluation of the expressions
    ( v1.x + a*v1.z - v2.x - a*v2.z ) and ( v1.y + b*v1.z - v2.y - b*v2.z ) with constraint
    constants substituted in.
    '''
    
    ## For a description of this function, see "note on normal/silhouette constraints" and the note
    ## on "oblique silhouettes" at the bottom of this file.
    
    print('gen_silhouette_constraint_oblique:', v2, v1)
    assert v2 != v1
    
    linear_constraints = []
    
    ## Gather vertices involved.
    v2v1he = mesh.halfedges[ mesh.directed_edge2he_index( ( v2, v1 ) ) ]
    v4 = mesh.halfedges[ v2v1he.next_he ].to_vertex
    v3 = mesh.halfedges[ mesh.halfedges[ v2v1he.opposite_he ].next_he ].to_vertex
    assert v3 != v2
    assert v4 != v2
    assert v4 != v3
    #s4 = 1. / mag( cross( v2v1_dir, mesh.vs[ v4 ] - mesh.vs[ v2 ] ) )
    #s3 = 1. / mag( cross( v2v1_dir, mesh.vs[ v2 ] - mesh.vs[ v3 ] ) )
    s4 = 1.
    s3 = 1.
    linear_constraints.append( (
        ## The list of ( vertex index, coord ) tuples.
        [
            (v4,1), (v4,2), (v2,1), (v2,2),
            (v4,0), (v4,2), (v2,0), (v2,2),
            (v2,1), (v2,2), (v3,1), (v3,2),
            (v2,0), (v2,2), (v3,0), (v3,2)
        ],
        ## The list of coefficients for each term.
        [
            s4 * v2_to_v1skew[0], s4 * oblique[1] * v2_to_v1skew[0], -s4 * v2_to_v1skew[0], -s4 * oblique[1] * v2_to_v1skew[0],
            -s4 * v2_to_v1skew[1], -s4 * oblique[0] * v2_to_v1skew[1], s4 * v2_to_v1skew[1], s4 * oblique[0] * v2_to_v1skew[1],
            s3 * v2_to_v1skew[0], s3 * oblique[1] * v2_to_v1skew[0], -s3 * v2_to_v1skew[0], -s3 * oblique[1] * v2_to_v1skew[0],
            -s3 * v2_to_v1skew[1], -s3 * oblique[0] * v2_to_v1skew[1], s3 * v2_to_v1skew[1], s3 * oblique[0] * v2_to_v1skew[1]
        ],
        ## The right hand side value.
        0.
        ) )
    
    return linear_constraints

def tangent_circle_between_rays( a, b, d ):
    '''
    Similar to tangent_circle_between_rays_t_a_minus_t_b().
    
    This function returns a tuple of two points, r_a( t_a ) and r_b( t_b ),
    the solution where t_b = 0.
    '''
    a = asarrayf( a )
    b = asarrayf( b )
    d = asarrayf( d )
    
    t_a_minus_t_b = tangent_circle_between_rays_t_a_minus_t_b( a, b, d )
    return a + t_a_minus_t_b * d, b

def tangent_circle_between_rays_t_a_minus_t_b( a, b, d ):
    '''
    Given two points 'a' and 'b' and a direction 'd', let us define two rays
        r_a(t_a) = a + t_a * d
    and
        r_b(t_b) = b + t_b * d
    
    We want to find a circle (or sphere) such that the rays are tangent to it and the points they
    touch are diametrically opposite each other.
    
    Following the derivation under the pages headed "Tilting circles" in my black notebook:
        I'm looking for a center c and radius r.
        If r_a and r_b are diametrically opposed, then c = 1/2 * ( r_a(t_a) + r_b(t_b) ) and
        r = | r_a(t_a) - r_b(t_b) |.
        Let us use the implicit sphere 1/2 * | p - c |^2 = 1/2 * r^2 for a given point p;
        the gradient with respect to p gives us the normal on the surface: p - c.
        If r_a and r_b are tangent to the circle (sphere), then
          dot( ( r_a - c ), d ) = 0 = dot( ( r_b - c ), d ) = dot( r_b - 1/2 * ( r_a + r_b ), d )
          = dot( 1/2 r_b - 1/2 r_a, d ) = 0 = dot( r_b - r_a, d )
        <=>
          dot( r_b, d ) = dot( r_a, d ) = dot( a + t_a * d, d ) = dot( b + t_b * d, d )
          = dot( a, d ) + t_a * dot( d, d ) = dot( b, d ) + t_b * dot( d, d )
        <=>
          dot( a - b, d ) = -( t_a - t_b ) * dot( d, d )
        <=>
          t_a - t_b = - dot( a - b, d ) / dot( d, d )
      
    This gives us the solution, up to translation along d.
    
    This function returns t_a - t_b.
    '''
    
    assert len( a ) == len( b )
    assert len( a ) == len( d )
    
    a = asarrayf( a )
    b = asarrayf( b )
    d = asarrayf( d )
    
    t_a_minus_t_b = - dot( a - b, d ) / dot( d, d )
    return t_a_minus_t_b

def straighten( curve, N = 100 ):
    '''
    The parameter 'curve' is a list of 2d points representing an open path (like the
    front-facing part of a cross section).
    Returns a straight line sampled 'N' times (default: 100) between the curve's endpoints.
    '''
    
    curve = asarray( curve )
    
    return [ curve[0] + float(i)/(N-1)*( curve[-1] - curve[0] ) for i in range(N) ]


def intersect_curve_with_mesh_edges( curve, mesh ):
    '''
    The parameter 'curve' is a list of 2d points representing an open path (like the
    front-facing part of a cross section).
    Returns a list of edges in TriMesh 'mesh' which intersect the curve.
    '''
    
    edges = []
    for edge in mesh.edges:
        for i in range( len( curve ) - 1 ):
            if line_segments_cross( ( mesh.vs[edge[0]], mesh.vs[edge[1]] ), ( curve[i], curve[i+1] ) ):
                edges.append( edge )
    
    return edges

def place_in_depth_according_to_connection( connection, frozen_regions, oblique = [0.,0.] ):
    '''
    Given a Connection 'connection',
    translates the connection.regions in z so that connection.regions[0] and connection.regions[1]
    intersect along connection.connection_curve.
    The set 'frozen_regions' contains the set of all regions whose z depth is already frozen in z.
    One of connection.regions must be in 'frozen_regions' and one not.
    
    This assumes that the curves in 'connection' have already been smoothed and self-intersections
    removed.
    
    tested
    '''
    
    #print 'place_in_depth_according_to_connection( oblique = %s )' % (oblique,)
    
    ## In this function, we are following place_in_depth_according_to_connection() from
    ## "annotations revisited data structures.txt":
    ## NOTE: This routine should be called on all connections in arbitrary order.  We must also
    ##       take care to move each region in z at most once.
    ##       We also assume that the graph of connections is a tree (or forest).  We can enforce this
    ##       by checking when the user tries to add a new connection.
    ## 
    ## 1 initially (if connection_curve_precise is not set),
    ##   move the two regions in z until their mesh intersection curve looks
    ##   as similar as possible to 'connection_curve'.
    ##   set 'connection_curve_precise' to the precise intersection curve.
    ## 
    ## NOTE: when drawing, the z buffer should naturally clip the two regions along this connection
    ##       curve.
    ## 
    ## 2 if the connection curve is smooth, solve for some surface offsets near the
    ##   'connection_curve_precise' projected onto both regions that smoothes the intersection.
    ## TODO Q: Given 'connection_curve_precise', should we join two connected regions' systems when
    ##         solving, or should we create a small system involving the inflated meshes' vertices
    ##         in an area around 'connection_curve_precise' and solve for offsets?
    
    assert len( connection.regions ) == 2
    assert ( connection.regions[0] in frozen_regions ) != ( connection.regions[1] in frozen_regions )
    
    free = [
        ( region, region_attaches_front )
        for region, region_attaches_front in zip( connection.regions, connection.regions_attaches_front )
        if region not in frozen_regions
        ][0]
    fixed = [
        ( region, region_attaches_front )
        for region, region_attaches_front in zip( connection.regions, connection.regions_attaches_front )
        if region in frozen_regions
        ][0]
    
    ## 1
    ## For every point along the connection curve, we want to minimize the difference in z
    ## between the two regions at the point.
    ## So our energy is
    ##   E = \sum_{ v \in connection_curve } ( v.z_free + z_{free offset} - v.z_fixed )^2
    ## To minimize this, we take the derivative of E with respect to our one free variable,
    ## z_{free offset}, and set it equal to 0:
    ##   dE/dz_{free offset} = 0 = \sum 2 * ( v.z_free + z_{free offset} - v.z_fixed )
    ##                       = \sum 2 * ( v.z_free - v.z_fixed ) + \sum 2 * z_{free offset}
    ## <=>
    ##   z_{free offset} = (-1/N) * \sum ( v.z_free - v.z_fixed )
    ## where N = len( connection_curve )
    
    ## To account for non-uniform sampling, we should multiply each vertex's contribution
    ## by half its incident edge lengths:
    ##   E = \sum_{ v \in connection_curve } ( v.z_free + z_{free offset} - v.z_fixed )^2 * .5 * ( len( v.edge_in ) + len( v.edge_out ) )
    ## To minimize this, we take the derivative of E with respect to our one free variable,
    ## z_{free offset}, and set it equal to 0:
    ##   dE/dz_{free offset} = 0 = \sum 2 * ( v.z_free + z_{free offset} - v.z_fixed ) * .5 * ( len( v.edge_in ) + len( v.edge_out ) )
    ##                       = \sum ( v.z_free - v.z_fixed ) * ( len( v.edge_in ) + len( v.edge_out ) ) + \sum z_{free offset} * ( len( v.edge_in ) + len( v.edge_out ) )
    ## <=>
    ##   z_{free offset} = (-1/( \sum len( v.edge_in ) + len( v.edge_out ) )) * \sum ( v.z_free - v.z_fixed ) * ( len( v.edge_in ) + len( v.edge_out ) )
    
    ## Data for each vertex along the connection curve.
    v_z_free = []
    v_z_fixed = []
    v_edge_in_len = []
    v_edge_out_len = []
    ## Precompute all edge lengths.
    edge_lengths = line_strip_seg_lengths( connection.connection_curve )
    ## Precompute offsets into the intersection array.
    free_isect_off = 0 if free[1] else 1
    fixed_isect_off = 0 if fixed[1] else 1
    
    ## Shoot the ray in the negative screen-space z direction.
    ## (Invert the (oblique) projection matrix to get the actual ray.)
    ray_dir = asarray( ( oblique[0], oblique[1], -1. ) )
    #ray_z_height = 1 + max( free[0].mesh.vs[:,2].max(), fixed[0].mesh.vs[:,2].max() )
    ray_z_height = 0.
    
    ## Iterate over connection curve points.
    for cpi, curve_pt in enumerate( connection.connection_curve ):
        ray_pt = append( curve_pt, ray_z_height )
        
        free_isects = free[0].ray_surface_intersection_ts( ( ray_pt, ray_dir ), ray_is_line = True )
        free_isects.sort()
        ## Skip the point if it doesn't intersect with one of the regions.
        if len( free_isects ) == 0: continue
        ## Because we have a closed surface, there should be an even number of intersections.
        expect( ( len( free_isects ) % 2 ) == 0 )
        ## Unless it folds back on itself, there should be exactly two intersections.
        ## TODO: If the surface folds back on itself, we have to choose which panel
        ##       to intersect with.  'free_isect_off' should choose this, but obeying
        ##       it slavishly can lead to a discontinuity: if the surface folds back on
        ##       itself (like a U shape), but the second leg of the U only covers part of
        ##       the first leg, then the front-most (0-th) panel will jump from the back
        ##       to the front leg.  Perhaps we should maintain continuity and choose the
        ##       'free_isect_off' panel that the connection curve mostly lies on, and
        ##       then use continuity of the intersection depth (and whether the number
        ##       of intersections changes) to choose which of the rest of the points'
        ##       depths to use (or ignore everything with a discontinuous jump).
        if not ( len( free_isects ) == 2 ): print("WARNING: Connection curve over a fold-over!")
        
        fixed_isects = fixed[0].ray_surface_intersection_ts( ( ray_pt, ray_dir ), ray_is_line = True )
        fixed_isects.sort()
        ## Skip the point if it doesn't intersect with one of the regions.
        if len( fixed_isects ) == 0: continue
        ## Because we have a closed surface, there should be an even number of intersections.
        expect( ( len( fixed_isects ) % 2 ) == 0 )
        ## Unless it folds back on itself, there should be exactly two intersections.
        ## TODO: See the long TODO above near expect( len( free_isects ) == 2 ).
        if not ( len( fixed_isects ) == 2 ): print("WARNING: Connection curve over a fold-over!")
        
        v_z_free.append( ray_pt[2] + free_isects[free_isect_off] * ray_dir[2] )
        v_z_fixed.append( ray_pt[2] + fixed_isects[fixed_isect_off] * ray_dir[2] )
        if cpi == 0:
            v_edge_in_len.append( 0. )
        else:
            v_edge_in_len.append( edge_lengths[cpi - 1] )
        
        if cpi == len( edge_lengths ):
            v_edge_out_len.append( 0. )
        else:
            v_edge_out_len.append( edge_lengths[cpi] )
    
    ## All 4 lists should have the same length.
    assert len( v_z_free ) == len( v_z_fixed )
    assert len( v_edge_in_len ) == len( v_edge_out_len )
    assert len( v_z_free ) == len( v_edge_in_len )
    
    ## In the case of a connection curve consisting of a single point, simply set the only
    ## in and out edge lengths to 1.  The lengths are only for normalization, so with a single
    ## point it doesn't matter what they are so long as they are non-zero.
    if len( edge_lengths ) == 0:
        assert 1 == len( v_edge_in_len )
        assert 1 == len( v_edge_out_len )
        
        v_edge_in_len[0] = 1.
        v_edge_out_len[0] = 1.
    
    ## Choose z_offset so that the connection_curve is close to the user-specified curve.
    ## z_offset = (-1/( \sum len( v.edge_in ) + len( v.edge_out ) )) * \sum ( v.z_free - v.z_fixed ) * ( len( v.edge_in ) + len( v.edge_out ) )
    z_offset = 0.
    for vi in range(len( v_z_free )):
        z_offset += ( v_z_free[ vi ] - v_z_fixed[ vi ] ) * ( v_edge_in_len[ vi ] + v_edge_out_len[ vi ] )
    ## When no points on the connection curve intersect both surfaces, the lists will be empty,
    ## and the sum() will be 0.
    ## Skipping this divide-by-zero is all we need to do to keep z_offset 0., and we can continue
    ## to execute the general case.
    if len( v_edge_in_len ) + len( v_edge_out_len ) > 0:
        z_offset *= -1./( sum( v_edge_in_len ) + sum( v_edge_out_len ) )
    
    ## Offset the free region by z_offset.
    free[0].offset_in_depth( -z_offset, oblique )
    
    ## Now that we've placed the region in depth, freeze it.
    frozen_regions.add( free[0] )
    
    ## 2
    ## TODO: smooth connections

def connection_frozen_com( connection, frozen_regions, oblique = [0.,0.] ):
    '''
    Given a Connection 'connection',
    returns the center-of-mass of the connection curve intersected with the frozen region shape.
    The set 'frozen_regions' contains the set of all regions whose z depth is already frozen in z.
    One of connection.regions must be in 'frozen_regions' and one not.
    
    This assumes that the curves in 'connection' have already been smoothed and self-intersections
    removed.
    
    used
    '''
    
    ## This routine follows place_in_depth_according_to_connection() until it has collected
    ## enough data, and then simply returns the weighted average.
    
    assert len( connection.regions ) == 2
    assert ( connection.regions[0] in frozen_regions ) != ( connection.regions[1] in frozen_regions )
    
    fixed = [
        ( region, region_attaches_front )
        for region, region_attaches_front in zip( connection.regions, connection.regions_attaches_front )
        if region in frozen_regions
        ][0]
    
    ## Data for each vertex along the connection curve.
    v_fixed = []
    v_edge_in_len = []
    v_edge_out_len = []
    ## Precompute all edge lengths.
    edge_lengths = line_strip_seg_lengths( connection.connection_curve )
    ## Precompute offsets into the intersection array.
    fixed_isect_off = 0 if fixed[1] else 1
    
    ## Shoot the ray in the negative screen-space z direction.
    ## (Invert the (oblique) projection matrix to get the actual ray.)
    ray_dir = asarrayf( ( oblique[0], oblique[1], -1. ) )
    ray_z_height = 0.
    
    ## Iterate over connection curve points.
    for cpi, curve_pt in enumerate( connection.connection_curve ):
        ray_pt = append( curve_pt, ray_z_height )
        
        fixed_isects = fixed[0].ray_surface_intersection_ts( ( ray_pt, ray_dir ), ray_is_line = True )
        fixed_isects.sort()
        ## Skip the point if it doesn't intersect with one of the regions.
        if len( fixed_isects ) == 0: continue
        ## Because we have a closed surface, there should be an even number of intersections.
        expect( ( len( fixed_isects ) % 2 ) == 0 )
        ## Unless it folds back on itself, there should be exactly two intersections.
        expect( len( fixed_isects ) == 2 )
        
        v_fixed.append( ray_pt + fixed_isects[fixed_isect_off] * ray_dir )
        if cpi == 0:
            v_edge_in_len.append( 0. )
        else:
            v_edge_in_len.append( edge_lengths[cpi - 1] )
        
        if cpi == len( edge_lengths ):
            v_edge_out_len.append( 0. )
        else:
            v_edge_out_len.append( edge_lengths[cpi] )
    
    ## All 3 lists should have the same length.
    assert len( v_edge_in_len ) == len( v_edge_out_len )
    assert len( v_fixed ) == len( v_edge_in_len )
    
    if len( v_fixed ) == 0:
        #raise RuntimeError, "Connection curve doesn't intersect fixed shape."
        print("ERROR: Connection curve doesn't intersect fixed shape.")
        return append( average( connection.connection_curve, axis = 0 ), 0. )
    
    ## In the case of a connection curve consisting of a single point, simply set the only
    ## in and out edge lengths to 1.  The lengths are only for normalization, so with a single
    ## point it doesn't matter what they are so long as they are non-zero.
    if len( edge_lengths ) == 0:
        assert 1 == len( v_edge_in_len )
        assert 1 == len( v_edge_out_len )
        
        v_edge_in_len[0] = 1.
        v_edge_out_len[0] = 1.
    
    return average( v_fixed, axis = 0, weights = .5 * ( asarray( v_edge_in_len ) + asarray( v_edge_out_len ) ) )

def shortest_path_on_mesh( mesh, source, dest, **kwargs ):
    '''
    Given a TriMesh 'mesh', source vertex index 'source', and
    destination vertex index 'dest'
    returns the list of vertex indices [ boundary1, ..., boundary2 ] such that there is an edge
    in the mesh between each successive pair of vertices and the path returned is the shortest path
    according to the (by default) euclidean distance between the vertices.
    If no such path exists, returns None.
    
    Optional named parameters:
        restrict_vertices, a set of vertices to restrict the search to
        metric, a function taking two adjacent vertex indices that returns the distance between them
    
    used
    '''
    
    kwargs.setdefault( 'restrict_vertices', list(range(len( mesh.vs ))) )
    #kwargs.setdefault( 'restrict_vertices', range( len( mesh.vs ) ) )
    kwargs.setdefault( 'metric', lambda v1, v2: mag( mesh.vs[v1] - mesh.vs[v2] ) )
    metric = kwargs[ 'metric' ]
    
    ## We follow the wikipedia article on Dijkstra's algorithm:
    ##   http://en.wikipedia.org/wiki/Dijkstra's_algorithm
    dist = {}
    previous = {}
    
    dist[ source ] = 0.
    Q = set( kwargs[ 'restrict_vertices' ] )
    while len( Q ) > 0:
        distQ = [ ( dist[v], v ) for v in Q if v in dist ]
        if len( distQ ) == 0: break
        u = min( distQ )[1]
        Q.remove( u )
        if u == dest: break
        for v in mesh.vertex_vertex_neighbors( u ):
            alt = dist[u] + metric( u,v )
            if v not in dist or alt < dist[v]:
                dist[v] = alt
                previous[v] = u
    
    # return previous
    
    ## Make sure we didn't fail.
    if dest not in previous: return None
    
    S = []
    u = dest
    while u in previous:
        S.insert( 0, u )
        u = previous[u]
    S.insert( 0, u )
    return S

## A note on normal/silhouette constraints:
'''
Given two vectors < x1, y1, z1 > and < x2, y2, z2 >, the cross product is the determinant of
    | i  j  k  |
    | x1 y1 z1 | = i ( y1 z2 - y2 z1 ) + j ( z1 x2 - z2 x1 ) + k ( x1 y2 - x2 y1 )
    | x2 y2 z2 |

We are interested in the edge normal of the edge between vertex v2 and v1.
Let v4 and v3 be the vertices opposite the edge (v2,v1) and (v1,v2) as pictured:

    v1
  / | \
v4  |  v3
  \ | /
    v2

Let the (unit) normal vector along the edge between v2 and v1 be the sum of the
unit face normals (v2,v1,v4) and (v1,v2,v3):
    n * s0 = s4 (v1-v2) x (v4-v2)  +  s3 (v1-v2) x (v2-v3),
where
    s4 is the reciprocal of the magnitude of (v1-v2) x (v4-v2) and
    s3 is the reciprocal of the magnitude of (v1-v2) x (v2-v3) and
    s0 is the reciprocal of the magnitude of the right-hand-side.
Then
    nx * s0 =   s4 ( v1 - v2 ).y ( v4 - v2 ).z  -  s4 ( v1 - v2 ).z ( v4 - v2 ).y
              + s3 ( v1 - v2 ).y ( v2 - v3 ).z  -  s3 ( v1 - v2 ).z ( v2 - v3 ).y
and
    ny * s0 =   s4 ( v1 - v2 ).z ( v4 - v2 ).x  -  s4 ( v1 - v2 ).x ( v4 - v2 ).z
              + s3 ( v1 - v2 ).z ( v2 - v3 ).x  -  s3 ( v1 - v2 ).x ( v2 - v3 ).z
and
    nz * s0 =   s4 ( v1 - v2 ).x ( v4 - v2 ).y  -  s4 ( v1 - v2 ).y ( v4 - v2 ).x
              + s3 ( v1 - v2 ).x ( v2 - v3 ).y  -  s3 ( v1 - v2 ).y ( v2 - v3 ).x

If v1 and v2 have their x, y, and z coordinates constrained (absolutely or relatively)
and we assume constant s0, s3, and s4, then these equations are linear in vertices v4 and v3.

The situation is simpler if we only want to constrain the edge to be a silhouette.
In that case, our constraint is nz = 0, so s0 falls away.  Assuming s3 ~= s4, we are left with:
    0 =   ( v1 - v2 ).x ( v4 - v2 ).y  -  ( v1 - v2 ).y ( v4 - v2 ).x
        + ( v1 - v2 ).x ( v2 - v3 ).y  -  ( v1 - v2 ).y ( v2 - v3 ).x
'''
## Oblique silhouettes
'''
Oblique silhouettes are more complex than orthogonal silhouettes.  Given an oblique skew [a,b],
as in
    x.screen = x.world + a*z.world
    y.screen = y.world + b*z.world,
    z.screen = z.world
we want the skewed geometry's normal to have z = 0 along oblique silhouette edges.
So, taking the nz * s0 line from above, we replace
    vi.x with vi.x + a*vi.z
and
    vi.y with vi.y + b*vi.z
Note that a constrained vertex vi have the constraint vi.x + a*vi.z equals a constant, which
we can substitute in to keep these constraints linear.
So, we have
    0 =   ( v1.x + a*v1.z - v2.x - a*v2.z ) ( v4.y + b*v4.z - v2.y - b*v2.z )  -  ( v1.y + b*v1.z - v2.y - b*v2.z ) ( v4.x + a*v4.z - v2.x - b*v2.z )
        + ( v1.x + a*v1.z - v2.x - a*v2.z ) ( v2.y + b*v2.z - v3.y - b*v3.z )  -  ( v1.y + b*v1.z - v2.y - b*v2.z ) ( v2.x + a*v2.z - v3.x - b*v3.z )
'''
