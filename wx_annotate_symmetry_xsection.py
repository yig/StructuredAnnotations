#!/usr/bin/env python

#import OpenGL
#OpenGL.ERROR_ON_COPY = True
#OpenGL.ERROR_CHECKING = False

from wx_annotate_shell import *
## I have to import wx_annotate_shell and refer to wx_annotate_shell.kDebugging
## because the imported kDebugging symbol is not a reference to the dict entry
## with the same name in wx_annotate_shell.
import wx_annotate_shell

from annotations_rotund import Connection
from helpers import line_loop_from_tilted_circle

class SymmetryXSectionShape( object ):
    def __init__( self, silhouette = None, symmetry_xsection = None, circular_xsections = None ):
        ## NOTE: I can't make [] the default parameter, because this always uses the same instance,
        ##       so the instance is shared by all objects instantiated with the default parameter,
        ##       and when I append to that instance, it is seen by all objects sharing it.
        if circular_xsections is None: circular_xsections = []
        
        ## Closed loop of 2d points, drawn by the user.
        self.silhouette = silhouette
        ## Closed loop of 2d points, drawn by the user.
        self.symmetry_xsection = symmetry_xsection
        ## Pair of 2d points, the picked anchor point and its destination.
        ## Alternately, this is a triplet of 2d points: the first two are diametrically opposed,
        ## and the third specifies its tilt.
        self.circular_xsections = circular_xsections
        
        ## Closed loop of 3d points.
        self.symmetry_xsection_3d = None
        ## The tuple: ( 3d point for circle center, radius, 3d vector for circle normal )
        self.circular_xsections_3d = None
        
        ## The inferred mesh.
        self.mesh = None
        
        self.update_3d()
    
    def generate_circular_xsections_3d_line_loop( self ):
        return [
            line_loop_from_tilted_circle( *circle )
            for circle in self.circular_xsections_3d
            ]
    
    circular_xsections_3d_line_loop = property( generate_circular_xsections_3d_line_loop )
    
    def update_3d( self ):
        #if not hasattr( self, '_decs' ): self._decs = []
        self._decs = []
        
        #self.update_3d_simple()
        
        self.update_3d_rot_after()
        self.update_3d_mesh()
    
    def update_3d_simple( self ):
        '''
        Sets self.symmetry_xsection_3d and self.circular_xsections_3d and self.mesh based
        on self.silhouette, self.symmetry_xsection, and self.circular_xsection.
        '''
        
        ## TODO: Take oblique into account
        oblique_dir = asarray( ( 0., 0., -1. ) )
        
        self.circular_xsections_3d = []
        for circular_xsection in self.circular_xsections:
            if len( circular_xsection ) == 2:
                print("## TODO: Find the right circular cross section the same way as in wx_annotate_rotund.")
            else:
                ## Find the circular cross section as if the first two 2d points of circular_xsection
                ## establish the diameter and the third the tilt.
                assert len( circular_xsection ) == 3
                from annotations_rotund import tangent_circle_between_rays
                diam1, diam2 = tangent_circle_between_rays(
                    append( circular_xsection[0], 0. ),
                    append( circular_xsection[1], 0. ),
                    oblique_dir
                    )
                center = .5 * ( diam1 + diam2 )
                radius = .5 * mag( diam1 - diam2 )
                
                ## We want the geodesic which passes through diam1, diam2, and
                ## the ray through the anchor's destination point.
                ## 1 Find the point on the sphere under (by projection) the anchor's destination point.
                anchor_dest_sphere_pt = closest_point_on_sphere_to_ray(
                    ( center, radius ),
                    ## Pull the ray's point back in the view direction so that the closest point is the one closer to the viewer.
                    ( append( circular_xsection[2], 0. ) - ( center[2] + radius + 1 )*oblique_dir, oblique_dir ),
                    ray_is_line = True
                    )[0]
                
                ## We store in self.circular_xsections_3d:
                ## ( 3d point for circle center, radius, 3d vector for circle normal )
                self.circular_xsections_3d.append( ( center, radius, dir( cross( diam1 - center, anchor_dest_sphere_pt - center ) ) ) )
        
        if self.symmetry_xsection is not None and len( self.circular_xsections_3d ) > 0:
            ## Find the 3d symmetry cross section by projecting each point onto the plane
            ## perpendicular to the circle and passing through its diameter.
            
            ## 1 Intersect the 2d curve of symmetry_xsection with the projection of
            ##   circular_xsection_3d.  This should give us 4 intersection points.
            ## 2 Solve for the z coordinates of the 2d intersection points.
            ## 3 Assume the user drew the front-facing half of symmetry_xsection before the
            ##   backfacing half.  Then the first two intersections (following the curve
            ##   symmetry_xsection) are for the front-facing half of symmetry_xsection and the
            ##   second two are for the back-facing half.  For the front-facing half, we want the
            ##   intersection closer to the viewer, and for the back-facing half we want the
            ##   intersection farther from the viewer.
            ## 4 Find the plane containing the 3d points front-facing intersection and
            ##   back-facing intersection, and having the circle's normal.
            ## 5 Project all 2d points in symmetry_xsection onto the plane to get symmetry_xsection_3d.
            
            assert len( self.circular_xsections_3d ) == 1
            circle_center, circle_radius, circle_normal = self.circular_xsections_3d[0]
            
            ## 1
            ## TODO: Take oblique into account (rotate so that oblique_dir is (0,0,1) and then drop z).
            orthogonal1 = orthogonal( circle_normal )
            ellipse = ( circle_center, [ circle_radius * orthogonal1, circle_radius * cross( circle_normal, orthogonal1 ) ] )
            ellipse2d = ( ellipse[0][:2], [ axis[:2] for axis in ellipse[1] ] )
            intersections = []
            for i in range(len( self.symmetry_xsection )):
                ray2d = (
                    asarray( self.symmetry_xsection[i] ),
                    asarray( self.symmetry_xsection[ (i+1) % len( self.symmetry_xsection ) ] ) - asarray( self.symmetry_xsection[i] )
                    )
                ts = ray_ellipse_intersection_ts( ellipse2d, ray2d )
                
                intersections.extend( [ ray2d[0] + t * ray2d[1] for t in ts if t >= 0 and t <= 1 ] )
            
            ## TODO: For numerical reasons this might not be true; if not, remove the points epsilon
            ##       away from each other so that it is true.
            assert len( intersections ) == 4
            
            ## 2
            for i in range(len( intersections )):
                isect_pt = intersections[i]
                cos_theta, sin_theta = linalg.solve( asarray( ellipse2d[1] ).T, isect_pt - ellipse2d[0] )
                intersections[i] = ellipse[0] + cos_theta * ellipse[1][0] + sin_theta * ellipse[1][1]
            
            ## 3
            front_pt = max( intersections[:2], key = lambda el: el[2] )
            back_pt = min( intersections[2:], key = lambda el: el[2] )
            
            ## 4
            plane_pt = front_pt
            plane_normal = dir( cross( front_pt - back_pt, circle_normal ) )
            
            ## 5
            self.symmetry_xsection_3d = [
                ray_plane_intersection( ( append( pt, 0. ), oblique_dir ), ( plane_pt, plane_normal ), ray_is_line = True )
                for pt in self.symmetry_xsection
                ]
            ## Smooth the 3D curve.
            ## NOTE: We pass simplify = False because we can't simplify (remove self-intersections from)
            ##       3d curves.
            from fibermesh import smooth_curve
            self.symmetry_xsection_3d = smooth_curve( self.symmetry_xsection_3d, smooth = 11, simplify = False )
    
    def update_3d_rot_after( self ):
        '''
        Sets self.symmetry_xsection_3d and self.circular_xsections_3d and self.mesh based
        on self.silhouette, self.symmetry_xsection, and self.circular_xsections.
        '''
        
        ## TODO: Take oblique into account
        oblique_dir = asarray( ( 0., 0., -1. ) )
        
        self.circular_xsections_3d = []
        self.circular_xsections_3d_unrot = []
        #self.circular_xsections_3d_rot_axis_angle = []
        for circular_xsection in self.circular_xsections:
            if len( circular_xsection ) == 2:
                print("## TODO: Find the right circular cross section the same way as in wx_annotate_rotund.")
            else:
                ## Find the circular cross section as if the first two 2d points of circular_xsection
                ## establish the diameter and the third the tilt.
                assert len( circular_xsection ) == 3
                from annotations_rotund import tangent_circle_between_rays
                diam1, diam2 = tangent_circle_between_rays(
                    append( circular_xsection[0], 0. ),
                    append( circular_xsection[1], 0. ),
                    oblique_dir
                    )
                center = .5 * ( diam1 + diam2 )
                radius = .5 * mag( diam1 - diam2 )
                
                ## We want the geodesic which passes through diam1, diam2, and
                ## the ray through the anchor's destination point.
                ## 1 Find the point on the sphere under (by projection) the anchor's destination point.
                anchor_dest = append( circular_xsection[2], 0. )
                anchor_dest_sphere_pt = closest_point_on_sphere_to_ray(
                    ( center, radius ),
                    ## Pull the ray's point back in the view direction so that the closest point is the one closer to the viewer.
                    ( anchor_dest - ( center[2] + radius + 1 )*oblique_dir, oblique_dir ),
                    ray_is_line = True
                    )[0]
                
                ## To find the unrotated circle point, pretend like anchor_dest is the midpoint of
                ## the two fixed circle points.
                '''
                anchor_untilted_dest_sphere_pt = closest_point_on_sphere_to_ray(
                    ( center, radius ),
                    ## Pull the ray's point back in the view direction so that the closest point is the one closer to the viewer.
                    ( center - ( center[2] + radius + 1 )*oblique_dir, oblique_dir ),
                    ray_is_line = True
                    )[0]
                '''
                ## Actually, this point is just the point on the surface of the sphere passing
                ## through the sphere center along the look vector.
                anchor_untilted_dest_sphere_pt = center - radius * dir( oblique_dir )
                
                ## We store in self.circular_xsections_3d:
                ## ( 3d point for circle center, radius, 3d vector for circle normal )
                self.circular_xsections_3d.append( ( center, radius, dir( cross( diam1 - center, anchor_dest_sphere_pt - center ) ) ) )
                self.circular_xsections_3d_unrot.append( ( center, radius, dir( cross( diam1 - center, anchor_untilted_dest_sphere_pt - center ) ) ) )
                #self.circular_xsections_3d_rot_axis_angle.append( get_axis_angle_rotation( anchor_src_midpoint_sphere_pt - center, anchor_dest_sphere_pt - center ) )
        
        if self.symmetry_xsection is not None and len( self.circular_xsections_3d ) > 0:
            assert len( self.circular_xsections_3d_unrot ) == len( self.circular_xsections_3d )
            
            ## Find the 3d symmetry cross section by projecting each point onto the plane
            ## perpendicular to the circle and passing through its diameter.
            
            ## 1 Intersect the 2d curve of symmetry_xsection with the projection of
            ##   circular_xsection_3d.  This should give us 4 intersection points.
            ## 2 Solve for the z coordinates of the 2d intersection points.
            ## 3 Assume the user drew the front-facing half of symmetry_xsection before the
            ##   backfacing half.  Then the first two intersections (following the curve
            ##   symmetry_xsection) are for the front-facing half of symmetry_xsection and the
            ##   second two are for the back-facing half.  For the front-facing half, we want the
            ##   intersection closer to the viewer, and for the back-facing half we want the
            ##   intersection farther from the viewer.
            ## 4 Find the plane containing the 3d points front-facing intersection and
            ##   back-facing intersection, and having the circle's normal.
            ## 5 Project all 2d points in symmetry_xsection onto the plane to get symmetry_xsection_3d.
            
            assert len( self.circular_xsections_3d ) == 1
            circle_center, circle_radius, circle_normal = self.circular_xsections_3d[0]
            ## Actually, use the unrotated circle's normal
            circle_normal = self.circular_xsections_3d_unrot[0][2]
            
            self._decs.append( [ 'line', ( 1, 0, 0 ), self.circular_xsections_3d[0][0], self.circular_xsections_3d[0][0] + .3*self.circular_xsections_3d[0][2] ] )
            self._decs.append( [ 'line', ( 1, 0, 0 ), self.circular_xsections_3d_unrot[0][0], self.circular_xsections_3d_unrot[0][0] + .3*self.circular_xsections_3d_unrot[0][2] ] )
            
            ## 1
            ## TODO: Take oblique into account (rotate so that oblique_dir is (0,0,1) and then drop z).
            intersections = []
            front_z = circle_center[2] + circle_radius + 1
            back_z = circle_center[2] - circle_radius - 1
            for i in range(len( self.symmetry_xsection )):
                z = front_z if len( intersections ) == 0 else back_z
                ray3d = (
                    append( self.symmetry_xsection[i], z ),
                    append( self.symmetry_xsection[ (i+1) % len( self.symmetry_xsection ) ], z ) - append( self.symmetry_xsection[i], z )
                    )
                ts = [ ray_plane_intersection_t( ray3d, ( circle_center, circle_normal ) ) ]
                
                intersections.extend( [ ray3d[0] + t * ray3d[1] for t in ts if t >= 0 and t <= 1 ] )
            
            ## TODO: For numerical reasons this might not be true; if not, remove the points epsilon
            ##       away from each other so that it is true.
            assert len( intersections ) == 2
            
            self._decs.append( [ 'point', ( 0, 0, 1 ) ] + list( intersections ) )
            
            ## 2
            for i in range(len( intersections )):
                isect_pt = intersections[i]
                intersections[i] = closest_point_on_sphere_to_ray( ( circle_center, circle_radius ), ( isect_pt, (0,0,1) ), ray_is_line = True )[0]
            
            ## 3
            front_pt = max( intersections, key = lambda el: el[2] )
            ## Sometimes there are two points, both with z = 0.  In this case min() and max()
            ## return the same point, which gives us a zero vector for the plane normal, and then
            ## NaNs appear.
            intersections.remove( front_pt )
            back_pt = min( intersections, key = lambda el: el[2] )
            
            self._decs.append( [ 'point', ( 0, 0, 1 ), front_pt ] )
            self._decs.append( [ 'point', ( 0, 0, 1 ), back_pt ] )
            
            ## 4
            plane_pt = front_pt
            plane_normal = dir( cross( front_pt - back_pt, circle_normal ) )
            self.symmetry_xsection_3d_plane = ( plane_pt, plane_normal )
            
            ## 5
            self.symmetry_xsection_3d_unrot = [
                ray_plane_intersection( ( append( pt, 0. ), oblique_dir ), ( plane_pt, plane_normal ), ray_is_line = True )
                for pt in self.symmetry_xsection
                ]
            ## Smooth the 3D curve.
            ## NOTE: We pass simplify = False because we can't simplify (remove self-intersections from)
            ##       3d curves.
            from fibermesh import smooth_curve
            self.symmetry_xsection_3d_unrot = smooth_curve( self.symmetry_xsection_3d_unrot, smooth = 11, simplify = False )
            
            ## Finally, rotate all points with the rotation that takes the unrotated circle to the
            ## rotated circle.
            rot_mat = get_rotation_matrix( *get_axis_angle_rotation(
                self.circular_xsections_3d_unrot[0][2], self.circular_xsections_3d[0][2]
                ) )
            self.symmetry_xsection_3d = [
                circle_center + rot_mat @ ( pt - circle_center )
                for pt in self.symmetry_xsection_3d_unrot
                ]
            
            self._decs.append( [ 'line', ( 0, 1, 0 ) ] + list( self.symmetry_xsection ) + [self.symmetry_xsection[0]] )
    
    def update_3d_mesh( self ):
        '''
        Generate self.mesh.
        '''
        
        self.update_3d_mesh_fibermesh()
        #self.update_3d_mesh_sweep_surface()
    
    def update_3d_mesh_sweep_surface( self ):
        if self.silhouette is None or len( self.circular_xsections_3d ) == 0 or self.symmetry_xsection is None: return
        
        import sweep_shape
        ## TODO Q: Should I pass rotated or unrotated circular cross-sections
        ##         (self.circular_xsections_3d or self.circular_xsections_3d_unrot)?
        self.mesh = sweep_shape.mesh_for_shape(
            self.silhouette,
            self.circular_xsections_3d,
            self.symmetry_xsection,
            None
            )
    
    def update_3d_mesh_fibermesh( self ):
        if self.symmetry_xsection_3d is None or len( self.circular_xsections_3d ) == 0: return
        
        assert self.symmetry_xsection_3d is not None
        assert self.symmetry_xsection_3d_unrot is not None
        assert len( self.circular_xsections_3d ) > 0
        assert len( self.circular_xsections_3d_unrot ) == len( self.circular_xsections_3d )
        
        ## TODO: Real oblique.
        oblique = [0.,0.]
        
        ## Viewed from the side, perpendicular to the symmetry plane, creating a mesh is a lot like
        ## traditional surface inflation with the symmetry curve as the z = 0 curve, and the
        ## circular cross sections as circles, and enforced symmetry.
        
        from fibermesh import gen_glued_mesh, attach_curve2d_to_mesh_vertices, Constraints
        
        rot_mat = get_rotation_matrix( *get_axis_angle_rotation(
                self.symmetry_xsection_3d_plane[1], (0,0,1)
                ) )
        symmetry_xsection_xy = [ ( rot_mat @ pt )[:2] for pt in self.symmetry_xsection_3d_unrot ]
        symmetry_xsection_z = ( rot_mat @ self.symmetry_xsection_3d_unrot[0] )[2]
        mesh, boundary_vertices, front_vertices, back_vertices = gen_glued_mesh( symmetry_xsection_xy, N = 200 )
        ## These sets should be disjoint.
        assert len( set( boundary_vertices ).intersection( set( front_vertices ) ) ) == 0
        assert len( set( boundary_vertices ).intersection( set( back_vertices ) ) ) == 0
        assert len( set( front_vertices ).intersection( set( back_vertices ) ) ) == 0
        mesh.write_OBJ( 'update_3d_mesh-glued.obj' )
        mesh.vs = asarray( mesh.vs )
        
        constraints = Constraints()
        
        ## Constraint the silhouette cross section
        boundary_vertices, boundary_vi_to_curve_vi = attach_curve2d_to_mesh_vertices( symmetry_xsection_xy, mesh, boundary_vertices )
        constraints.xy_constraint_vertices = list( boundary_vertices )
        constraints.xy_constraint_positions = [ symmetry_xsection_xy[ si ] for si in boundary_vi_to_curve_vi ]
        ## The silhouette lives on a constant z plane.
        constraints.z_constraint_vertices = list( boundary_vertices )
        constraints.z_constraint_values = [0.] * len( boundary_vertices )
        
        ## Symmetry means that the sum of matching front and back vertices' z coordinates is 0.
        for front_vi, back_vi in zip( front_vertices, back_vertices ):
            constraints.z_vertex_weighted_sum_indices.append( [ front_vi, back_vi ] )
            constraints.z_vertex_weighted_sum_weights.append( [ 1., 1. ] )
            constraints.z_vertex_weighted_sum_rhss.append( 0. )
        
        ## Run 4 iterations of fibermesh inflation.
        from fibermesh_cpp import optimize_surface2 as optimize_surface
        optimize_surface( mesh, constraints, max_iterations = 4, graph_not_cotangent_weights = True, vertex_normal_motion_only = False )
        
        self._decs.append( [ 'line',  ( 0, 0, 0 ) ] + list( symmetry_xsection_xy ) + [symmetry_xsection_xy[-1]] )
        
        ## For each cross section, project intersected mesh edges onto the cylinder implied by the
        ## (unrotated) circle.
        constrained_edges = set()
        for center, radius, circle_normal in self.circular_xsections_3d_unrot:
            
            ## Rotate the circular cross-section to match the rotation of the symmetry cross-section.
            center = rot_mat @ asarray( center )
            ## The (un-tilted) circular cross section normal should remain unchanged.
            #circle_normal = rot_mat @ asarray( circle_normal )
            assert mag2( asarray( circle_normal ) - rot_mat @ asarray( circle_normal ) ) < kEps or debugger()
            ## It should also lie in the image plane.
            assert abs( circle_normal[2] ) < kEps
            
            axis1 = orthogonal( circle_normal )
            axis2 = cross( axis1, circle_normal )
            
            ## Find the intersected mesh edges.
            ## TODO Q: Should I be doing something much more direct, connecting a spine
            ##         from the circle center to the top and bottom of the symmetry cross-section
            ##         (where it intersects the silhouette) and then sweeping a circle along
            ##         the spine?
            import annotations_rotund
            try:
                boundary1, cs_front_edges, boundary2, cs_back_edges = annotations_rotund.get_cross_section_edgesA(
                    center[:2], asarray( [ -circle_normal[1], circle_normal[0] ] ),
                    mesh, boundary_vertices, front_vertices, back_vertices,
                    oblique
                    )
            except RuntimeError:
                print('''TODO:
                ## This can happen if the intersection includes an edge with one vertex a rejected
                ## boundary vertex (rejected by attach_curve2d_to_mesh_vertices()).
                ''')
                continue
            expect( len( cs_front_edges ) > 0 )
            expect( len( cs_back_edges ) > 0 )
            
            ## Make sure we're not overconstraining edges.
            cross_section_all_edges = set()
            cross_section_all_edges.add( frozenset( boundary1[0] ) )
            cross_section_all_edges.add( frozenset( boundary2[0] ) )
            cross_section_all_edges.update( [ frozenset( edge_t[0] ) for edge_t in cs_front_edges ] )
            cross_section_all_edges.update( [ frozenset( edge_t[0] ) for edge_t in cs_back_edges ] )
            if len( constrained_edges.intersection( cross_section_all_edges ) ) > 0:
                print('cross section intersects another cross section; skipping:', cross_section)
                continue
            
            constrained_edges.update( cross_section_all_edges )
            
            boundary1pos = lerp( mesh.vs[ boundary1[0][0] ], mesh.vs[ boundary1[0][1] ], boundary1[1] )
            boundary2pos = lerp( mesh.vs[ boundary2[0][0] ], mesh.vs[ boundary2[0][1] ], boundary2[1] )
            
            ## The circle may not be centered in the symmetry cross section; center it.
            center = .5 * ( boundary1pos + boundary2pos )
            
            for edge_t in ( cs_front_edges + cs_back_edges ):
                edge0 = mesh.vs[ edge_t[0][0] ]
                edge1 = mesh.vs[ edge_t[0][1] ]
                
                self._decs.append( (
                    'line',
                    ( 1, 0, 0 ),
                    #dot( rot_mat.T, edge0 ),
                    #dot( rot_mat.T, edge1 )
                    edge0, edge1
                    ) )
                
                ## Constrain the edge so that it lies on the cylinder through the circular
                ## cross section (with the circle's normal as the cylinder's axis).
                
                edge0proj = (
                    center +
                    radius * dir(
                        dot( edge0 - center, axis1 ) * axis1 + dot( edge0 - center, axis2 ) * axis2
                        ) +
                    dot( edge0 - center, circle_normal ) * circle_normal
                    )
                edge1proj = (
                    center +
                    radius * dir(
                        dot( edge1 - center, axis1 ) * axis1 + dot( edge1 - center, axis2 ) * axis2
                        ) +
                    dot( edge1 - center, circle_normal ) * circle_normal
                    )
                
                ## Constrain the edge to this vector.
                constraints.xy_vertex_weighted_sum_indices.append( edge_t[0] )
                constraints.xy_vertex_weighted_sum_weights.append( [1, -1] )
                constraints.xy_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[:2] )
                constraints.z_vertex_weighted_sum_indices.append( edge_t[0] )
                constraints.z_vertex_weighted_sum_weights.append( [1, -1] )
                constraints.z_vertex_weighted_sum_rhss.append( ( edge0proj - edge1proj )[2] )
                
                self._decs.append( (
                    'line',
                    ( .5, .5, 0 ),
                    #dot( rot_mat.T, edge0proj ),
                    #dot( rot_mat.T, edge1proj )
                    edge0proj, edge1proj
                    ) )
        
        assert len( self.circular_xsections_3d_unrot ) > 0
        ## Run 4 more fibermesh iterations.
        optimize_surface( mesh, constraints, max_iterations = 4, graph_not_cotangent_weights = True, vertex_normal_motion_only = False )
        ## Now that we have a circular cross-section, we can run a linear solve.
        from fibermesh_cpp import optimize_surface_linear
        #optimize_surface_linear( mesh, constraints, graph_not_cotangent_weights = False, vertex_normal_motion_only = False )
        
        ## Now apply the inverse rotation to the mesh.
        mesh.vs[:,2] += symmetry_xsection_z
        mesh.vs = ( rot_mat.T @ mesh.vs.T ).T
        ## Finally, rotate all points with the rotation that takes the unrotated circle to the
        ## rotated circle.
        ## TODO: Average rotation of all circles. (See sweep_shape.symmetry_plane_from_tilted_circles() for code which averages planes.)
        if len( self.circular_xsections_3d_unrot ) > 1: print("WARNING: Only looking at the first tilted circle.")
        rot_mat = get_rotation_matrix( *get_axis_angle_rotation(
            self.circular_xsections_3d_unrot[0][2], self.circular_xsections_3d[0][2]
            ) )
        mesh.vs = ( rot_mat @ ( mesh.vs - asarray( self.circular_xsections_3d[0][0] ) ).T ).T + asarray( self.circular_xsections_3d[0][0] )
        
        self.mesh = mesh
        self.mesh_original_v0 = array( self.mesh.vs[0] )
    
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
        
        self.mesh.vs[:,:] += z_offset * array( ( oblique[0], oblique[1], -1 ) )
        self.mesh.positions_changed()


def PickShapes( anns, world_ray ):
    '''
    Given a list of annotations 'anns' and a ray (3d point, 3d direction),
    returns all SymmetryXSectionShape objects picked
    by the ray, sorted in depth (first intersection is first).
    
    NOTE: To create the world ray from a point on the screen, convert the point with
          screen2world() and use 0 as the z component.  The direction should be
          ( oblique a, oblique b, -1 ).
    '''
    
    assert len( world_ray[0] ) == 3
    assert len( world_ray[1] ) == 3
    world_ray = ( asarrayf( world_ray[0] ).ravel(), asarrayf( world_ray[1] ).ravel() )
    
    from surface_intersections import ray_mesh_intersections
    t2shape = []
    for shape in anns.of_type( SymmetryXSectionShape ):
        if shape.mesh is None: continue
        
        isects = ray_mesh_intersections( world_ray, shape.mesh )
        if len( isects ) == 0: continue
        
        minisect = min( isects )
        t2shape.append( ( minisect[0], shape ) )
    t2shape.sort()
    return [ shape for t, shape in t2shape ]

def PickShape( anns, world_ray ):
    '''
    Given a list of annotations 'anns' and a ray (3d point, 3d direction),
    returns the closest SymmetryXSectionShape object
    picked by the ray, sorted in depth (first intersection is first).
    If no SymmetryXSectionShape object is picked, returns None.
    
    NOTE: To create the world ray from a point on the screen, convert the point with
          screen2world() and use 0 as the z component.  The direction should be
          ( oblique a, oblique b, -1 ).
    '''
    
    shapes = PickShapes( anns, world_ray )
    if len( shapes ) == 0: return None
    else: return shapes[0]

Color = {
    'mesh': ( 200, 200, 255 ),
    'silhouette': ( 234, 145, 248 ),
    'symmetry_xsection': ( 174//2, 237//2, 111//2 ),
    'circular_xsection': ( 213, 99, 191 ),
    'connection': ( 174//2, 237//2, 111//2 )
    }

def DrawShape( prepare2D, prepare3D, shape, what_to_draw = None ):
    if what_to_draw is None: what_to_draw = list(Color.keys())
    
    prepare3D()
    
    if 'mesh' in what_to_draw:
        ## Draw the interior
        ## TODO Q: Should I use draw it 2d, using gluTessPolygon (see my layers2d app for reference)
        if shape.mesh is not None:
            #glDisable( GL_DEPTH_TEST )
            glColor4ub( *( Color['mesh'] + (50,) ) )
            ## draw_mesh_faces_smooth() is much faster than draw_mesh_faces_flat()
            draw_mesh_faces_smooth( shape.mesh, True )
            #glEnable( GL_DEPTH_TEST )
    
    if 'silhouette' in what_to_draw:
        if shape.mesh is not None:
            glDisable( GL_LIGHTING )
            glLineWidth( 8 )
            glColor4ub( *( Color['silhouette'] + (200,) ) )
            
            ## TODO: Use weakref.WeakKeyDictionary()
            if not hasattr( shape.mesh, '_silhouette_edges' ) or shape.mesh._silhouette_edges[0] != shape.mesh.lifetime_counter:
                from surface_intersections import mesh_silhouette_edges
                mesh_skew = TriMesh()
                mesh_skew.faces = array( shape.mesh.faces )
                skew = identity(3)
                oblique = [ 0., 0. ]
                skew[0,2] = oblique[0]
                skew[1,2] = oblique[1]
                mesh_skew.vs = dot( skew, asarray( shape.mesh.vs ).T ).T
                mesh_skew.topology_changed()
                shape.mesh._silhouette_edges = shape.mesh.lifetime_counter, mesh_silhouette_edges( mesh_skew )
            
            pts = []
            for edge in shape.mesh._silhouette_edges[1]:
                pts.append( shape.mesh.vs[ edge[0] ] )
                pts.append( shape.mesh.vs[ edge[1] ] )
            draw_lines( pts )
            glEnable( GL_LIGHTING )
    
    ## Draw the boundary stroke.
    if shape.silhouette is not None:
        assert len( shape.silhouette ) > 2
        
        prepare2D()
        
        if 'silhouette' in what_to_draw:
            ## Draw silhouette edges
            glLineWidth( 8 )
            glColor4ub( *( Color['silhouette'] + (50,) ) )
            draw_lineloop( shape.silhouette )
    
    
    cs3ds = []
    if 'circular_xsection' in what_to_draw and shape.circular_xsections_3d:
        for circular_xsection_3d_line_loop in shape.circular_xsections_3d_line_loop:
            cs3ds.append( ( 'circular_xsection', circular_xsection_3d_line_loop ) )
    
    if 'symmetry_xsection' in what_to_draw and shape.symmetry_xsection_3d:
        cs3ds.append( ( 'symmetry_xsection', shape.symmetry_xsection_3d ) )
    
    ## Draw cross sections.
    if cs3ds:
        prepare3D()
        glDisable( GL_LIGHTING )
        
        ## Offset the way the mesh was offset (by placement in depth)
        glMatrixMode( GL_MODELVIEW )
        glPushMatrix()
        if shape.mesh is not None and hasattr( shape, "mesh_original_v0" ):
            glTranslate( *( shape.mesh.vs[0] - shape.mesh_original_v0 ) )
        
        glDisable( GL_DEPTH_TEST )
        ## This doesn't work so well at making only the front of the cross section visible,
        ## because the cross section is independent of the geometry it's intersecting, so
        ## we might need quite a bit of an offset.
        ## Instead, we'll disable the depth test and set the alpha to the z of each point.
        #glMatrixMode( GL_MODELVIEW )
        #glPushMatrix()
        #glTranslated( *( -.05 * append( self.annotations.oblique, -1 ) ) )
        
        glLineWidth( 3 )
        #glColor4ub( *( Color['circular_xsection'] + (200,) ) )
        #for cs3d in region.cross_sections3d: draw_lineloop( cs3d )
        
        for kind, cs3d in cs3ds:
            colors = empty( ( len(cs3d), 4 ), dtype = ubyte )
            colors[:,0:3] = Color[kind]
            
            cs3d = asarray( cs3d )
            minz = cs3d[:,2].min()
            maxz = cs3d[:,2].max()
            if maxz == minz: colors[:,3] = 200
            else: colors[:,3] = 200 * ( (cs3d[:,2] - minz) / (maxz-minz) )
            
            glEnableClientState( GL_COLOR_ARRAY )
            glColorPointer( len(colors[0]), GL_UNSIGNED_BYTE, 0, asarray( colors ) )
            
            ## OK, I know this is implemented by with glVertexPointer, but it was
            ## until a moment ago implemented using GL_LINE_STRIP and tacking the first
            ## point onto the end.  Better to just put the four lines of code here,
            ## to ensure consistency.
            #draw_lineloop( cs3d )
            glEnableClientState( GL_VERTEX_ARRAY )
            glVertexPointer( len(cs3d[0]), GL_DOUBLE, 0, asarray( cs3d ) )
            glDrawArrays( GL_LINE_LOOP, 0, len(cs3d) )
            glDisableClientState( GL_VERTEX_ARRAY )
            
            glDisableClientState( GL_COLOR_ARRAY )
        
        #glPopMatrix()
        
        glMatrixMode( GL_MODELVIEW )
        glPopMatrix()
    

def DrawConnection( prepare2D, prepare3D, connection, alpha = 200 ):
    glLineWidth( 8 )
    glColor4ub( *( Color['connection'] + (alpha,) ) )
    if connection.connection_curve_precise is not None:
        prepare3D()
        glDisable( GL_LIGHTING )
        glDisable( GL_DEPTH_TEST )
        draw_linestrip( connection.connection_curve_precise )
    else:
        prepare2D()
        draw_linestrip( connection.connection_curve )

def DrawAnnotations( prepare2D, prepare3D, anns, opts ):
    ## First draw the rotund regions, which sets the zbuffer
    for shape in anns.of_type( SymmetryXSectionShape ):
        DrawShape( prepare2D, prepare3D, shape, what_to_draw = 'mesh' )
    
    ## Now draw the decorating lines
    decs = set( Color.keys() )
    decs.remove( 'mesh' )
    for shape in anns.of_type( SymmetryXSectionShape ):
        DrawShape( prepare2D, prepare3D, shape, what_to_draw = decs )
    
    ## Draw the connections either lightly or darkly, depending on whether we're in
    ## connecting mode
    connection_alpha = 200 if type( opts.mode ) is ConnectionCurveDrawMode else 50
    for connection in anns.of_type( Connection ):
        DrawConnection( prepare2D, prepare3D, connection, alpha = connection_alpha )
    
    if wx_annotate_shell.kDebugging:
        DrawDebugging( prepare2D, prepare3D, anns, opts )
    
def DrawDebugging( prepare2D, prepare3D, anns, opts = None ):
    prepare3D()
    
    glEnable( GL_DEPTH_TEST )
    glDisable( GL_LIGHTING )
    for shape in anns:
        if hasattr( shape, '_decs' ) and shape._decs is not None:
            for dec in shape._decs:
                #print 'dec:', dec
                if dec[0] == 'line':
                    glLineWidth( 4 )
                    [ None, None, None, glColor3f, glColor4f ][ len( dec[1] ) ]( *dec[1] )
                    draw_linestrip( dec[2:] )
                elif dec[0] == 'point':
                    glPointSize( 8 )
                    [ None, None, None, glColor3f, glColor4f ][ len( dec[1] ) ]( *dec[1] )
                    glBegin( GL_POINTS )
                    for pt in dec[2:]:
                        [ None, None, glVertex2d, glVertex3d, glVertex4d ][ len( pt ) ]( *pt )
                    glEnd()


def PlaceInDepth( anns, oblique = [0.,0.] ):
    ## TODO: proper oblique
    
    from annotations_rotund import place_in_depth_according_to_connection
    
    frozen_regions = set()
    
    all_regions = []
    for connection in anns.of_type( Connection ):
        all_regions.extend( connection.regions )
    
    ## A map from a region to connections touching it.
    region2connections = {}
    for region in all_regions:
        region2connections.setdefault( region, [] )
    for connection in anns.of_type( Connection ):
        for region in connection.regions:
            region2connections[ region ].append( connection )
    
    def dfs( region ):
        assert region in frozen_regions
        
        for connection in region2connections[ region ]:
            assert len( connection.regions ) == 2
            if connection.regions[0] in frozen_regions and connection.regions[1] in frozen_regions: continue
            place_in_depth_according_to_connection( connection, frozen_regions, oblique )
            dfs( [ child for child in connection.regions if child != region ][0] )
    
    for region in all_regions:
        if region not in frozen_regions:
            frozen_regions.add( region )
            dfs( region )


def AddConnection( anns, connection ):
    assert connection not in anns
    assert False not in [ shape in anns for shape in connection.regions ]
    
    anns.append( connection )

def DeleteConnection( anns, connection ):
    assert connection in anns
    anns.remove( connection )

def AddSymmetryXSectionShape( anns, silhouette, *args, **kwargs ):
    ## Reset if not enough points in the stroke.
    if len( silhouette ) < 3:
        print("AddSymmetryXSectionShape(): not enough points in the stroke.")
        return
    
    from helpers import smooth_curve
    silhouette = smooth_curve( silhouette, smooth = 11 )
    _AddPerfectSymmetryXSectionShape( anns, silhouette, *args, **kwargs )

def DeleteSymmetryXSectionShape( anns, shape ):
    anns.remove( shape )
    
    ## Also delete any connections which reference this region.
    connections_delete_me = []
    for connection in anns.of_type( Connection ):
        if shape in connection.regions:
            connections_delete_me.append( connection )
    for con in connections_delete_me: anns.remove( con )

def _AddPerfectSymmetryXSectionShape( anns, *args, **kwargs ):
    shape = SymmetryXSectionShape( *args, **kwargs )
    
    anns.append( shape )
    #self.GetTopLevelParent().Refresh()

def AddPerfectCircleRegion( anns ):
    
    silhouette = [ ( .5 * cos( theta ), .5 * sin( theta ) ) for theta in linspace( 0, 2*pi, 100 ) ]
    
    _AddPerfectSymmetryXSectionShape( anns, silhouette )

def AddOvalRegion( anns ):
    
    silhouette = [ ( cos( theta ), .5 * sin( theta ) ) for theta in linspace( 0, 2*pi, 100 ) ]
    ## Scale by .9
    silhouette = ( .9 * asarrayf( silhouette ) ).tolist()
    
    _AddPerfectSymmetryXSectionShape( anns, silhouette )

def AddRoundedCylinderRegion( anns ):
    
    silhouette = (
        [ ( .5 * cos( theta ) - .5, .5 * sin( theta ) ) for theta in linspace( pi/2, 3*pi/2, 100 ) ] +
        list(zip( linspace( -.5, .5, 10 ), -.5 * ones(10) ))[1:-1] +
        [ ( .5 * cos( theta ) + .5, .5 * sin( theta ) ) for theta in linspace( -pi/2, pi/2, 100 ) ] +
        list(zip( linspace( .5, -.5, 10 ), .5 * ones(10) ))[1:-1]
        )
    ## Scale by .9
    silhouette = ( .9 * asarrayf( silhouette ) ).tolist()
    
    _AddPerfectSymmetryXSectionShape( anns, silhouette )


class AnnotationCreateMode( ModeProxy ):
    '''
    Mode to begin a new annotation; user draws a silhouette.
    '''
    
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.stroke = None
        self.stroke_handler = None
    
    def doMouseDown( self, evt ):
        ## Forward events to the stroke drawing class.
        self.stroke = []
        if self.stroke_handler: unforward_events( self.stroke_handler )
        self.stroke_handler = StrokeHandler(
            self.stroke,
            lambda pt: self.Refresh(),
            self.doStrokeHandlerFinished,
            self.event2world,
            debugging = self.debugging
            )
        forward_events( self, self.stroke_handler )
        ## Now that we've forwarded events, call doMouseDown() again.  This time it will forward.
        return self.doMouseDown( evt )
    
    def doStrokeHandlerFinished( self, stroke ):
        self.stroke = None
        unforward_events( self.stroke_handler )
        self.stroke_handler = None
        
        AddSymmetryXSectionShape( self.annotations.anns, stroke )
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
        ## in wxPython.
        if keycode in (wx.WXK_DELETE,8):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            try:
                ## TODO: Use proper oblique direction.
                shape = PickShape( self.annotations.anns, ( append( self.screen2world( screen_pos ), 0. ), (0,0,-1) ) )
                print('deleting a shape')
                DeleteSymmetryXSectionShape( self.annotations.anns, shape )
                
                PlaceInDepth( self.annotations.anns )
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        else:
            evt.Skip()
        
        return self
    
    def doDraw( self ):
        def draw_stroke():
            if self.stroke is None: return
            
            glColor4ub( *( Color[ 'silhouette' ] + (200,) ) )
            glLineWidth( 8 )
            draw_linestrip( self.stroke )
        
        self.DrawScene( callback = draw_stroke )


class AnnotationCreateMode2( ModeProxy ):
    '''
    Mode to begin a new annotation; user draws an ellipse which is the circle, and then draws the
    symmetry cross section.
    '''
    
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.stroke = None
        self.stroke_handler = None
        self.shape = None
    
    def doMouseDown( self, evt ):
        if self.shape is None:
            ## Append the new, empty shape.
            self.shape = SymmetryXSectionShape()
            self.annotations.anns.append( self.shape )
            
            ## Forward events to the stroke drawing class.
            self.stroke = []
            if self.stroke_handler: unforward_events( self.stroke_handler )
            self.stroke_handler = StrokeHandler(
                self.stroke,
                self.doCircularPointAdded,
                self.doCircularStrokeHandlerFinished,
                self.event2world,
                debugging = self.debugging
                )
        else:
            ## Draw the symmetry cross section.
            self.shape.symmetry_xsection = []
            ## Forward events to the stroke drawing class.
            self.stroke = self.shape.symmetry_xsection
            if self.stroke_handler: unforward_events( self.stroke_handler )
            self.stroke_handler = StrokeHandler(
                self.stroke,
                self.doSymmetryPointAdded,
                self.doSymmetryStrokeHandlerFinished,
                self.event2world,
                debugging = self.debugging
                )
        
        forward_events( self, self.stroke_handler )
        ## Now that we've forwarded events, call doMouseDown() again.  This time it will forward.
        return self.doMouseDown( evt )
    
    def doCircularPointAdded( self, pt ):
        ellipse = lsq_ellipse( self.stroke )
        if ellipse is None: return
        
        ellipse = list(map( asarray, ellipse ))
        
        ## We always want the ellipse tilted down.
        s = 1 if ellipse[1][1][1] < 0 else -1
        self.shape.circular_xsections = [(
            ellipse[0] - ellipse[1][0], ellipse[0] + ellipse[1][0], ellipse[0] + s*ellipse[1][1]
            )]
        self.shape.update_3d()
        
        self.GetTopLevelParent().Refresh()
    
    def doCircularStrokeHandlerFinished( self, stroke ):
        self.stroke = None
        unforward_events( self.stroke_handler )
        self.stroke_handler = None
        
        self.GetTopLevelParent().Refresh()
    
    def doSymmetryPointAdded( self, pt ):
        #self.shape.update_3d()
        self.GetTopLevelParent().Refresh()
    
    def doSymmetryStrokeHandlerFinished( self, stroke ):
        self.stroke = None
        unforward_events( self.stroke_handler )
        self.stroke_handler = None
        
        from helpers import smooth_curve
        self.shape.symmetry_xsection = smooth_curve( self.shape.symmetry_xsection, smooth = 1 )
        self.shape.update_3d()
        ## Now forget about the shape so next time we start drawing a new one.
        #self.shape = None
        
        self.GetTopLevelParent().Refresh()
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
        ## in wxPython.
        if keycode in (wx.WXK_DELETE,8):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            try:
                ## TODO: Use proper oblique direction.
                shape = PickShape( self.annotations.anns, ( append( self.screen2world( screen_pos ), 0. ), (0,0,-1) ) )
                print('deleting a shape')
                DeleteSymmetryXSectionShape( self.annotations.anns, shape )
                
                PlaceInDepth( self.annotations.anns )
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        else:
            evt.Skip()
        
        return self
    
    def doDraw( self ):
        def draw_stroke():
            if self.stroke is None: return
            
            if self.shape.symmetry_xsection is None:
                glColor4ub( *( Color[ 'circular_xsection' ] + (50,) ) )
                glLineWidth( 1 )
            else:
                glColor4ub( *( Color[ 'symmetry_xsection' ] + (100,) ) )
                glLineWidth( 4 )
            draw_linestrip( self.stroke )
        
        self.DrawScene( callback = draw_stroke )

class ConnectionCurveDrawMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.stroke = None
        self.stroke_handler = None
    
    def doMouseDown( self, evt ):
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return self
        
        screen_pos = self.event2screen( self.GetMousePosition() )
        ## If we don't pick exactly two regions under the mouse, don't draw a connection.
        shapes = PickShapes( self.annotations.anns, ( append( self.screen2world( screen_pos ), 0. ), (0,0,-1) ) )
        if len( shapes ) != 2: return False
        ## If there's already a connection between the two regions, overdraw it.
        ## TODO Q: Should this drag the connection curve instead?
        previous_connection = [ con for con in self.annotations.anns.of_type( Connection ) if set( con.regions ) == set( shapes ) ]
        assert len( previous_connection ) <= 1
        previous_connection = None if len( previous_connection ) == 0 else previous_connection[0]
        
        ## TODO: Check if we would create a cycle in the connection graph; if so, bail.
        
        ## Forward events to the stroke drawing class.
        self.stroke = []
        if self.stroke_handler: unforward_events( self.stroke_handler )
        self.stroke_handler = StrokeHandler(
            self.stroke,
            lambda pt: self.Refresh(),
            lambda stroke: self.doStrokeHandlerFinished( shapes, previous_connection, stroke ),
            self.event2world,
            debugging = self.debugging
            )
        forward_events( self, self.stroke_handler )
        ## Now that we've forwarded events, call doMouseDown() again.  This time it will forward.
        return self.doMouseDown( evt )
    
    def doStrokeHandlerFinished( self, shapes, previous_connection, stroke ):
        assert len( shapes ) == 2
        
        self.stroke = None
        unforward_events( self.stroke_handler )
        self.stroke_handler = None
        
        if len( stroke ) > 1 and tuple( stroke[0] ) != tuple( stroke[-1] ):
            
            stroke = smooth_open_stroke( stroke, resample_N = 10 )
            if previous_connection:
                previous_connection.ReplaceConnectionCurve( stroke )
            else:
                from annotations_rotund import Connection
                connection = Connection( stroke, shapes[0], True, shapes[1], True )
                AddConnection( self.annotations.anns, connection )
            
            PlaceInDepth( self.annotations.anns )
            self.GetTopLevelParent().Refresh()
    
    def doContextMenu( self, evt ):
        evt_pos = self.GetMousePosition()
        #print 'context menu:', self.event2screen( evt_pos )
        
        screen_pos = self.event2screen( evt_pos )
        shapes = PickShapes( self.annotations.anns, ( append( self.screen2world( screen_pos ), 0. ), (0,0,-1) ) )
        #print 'context menu picked regions:', regions
        ## If we don't pick exactly two regions, bail.
        if len( shapes ) != 2:
            evt.Skip()
            return self
        
        connection = None
        for con in self.annotations.anns.of_type( Connection ):
            if set( con.regions ) == set( shapes ):
                connection = con
                break
        #print 'context menu picked connection:', connection
        ## If the regions don't have a connection between them, bail.
        if connection is None:
            evt.Skip()
            return self
        
        
        ## Show a contextual menu for the closest two overlapping annotations.
        menu = wx.Menu()
        menu_ids = []
        
        def disconnect( evt ):
            DeleteConnection( self.annotations.anns, connection )
            
            PlaceInDepth( self.annotations.anns )
            self.GetTopLevelParent().Refresh()
        
        item = menu.Append( wx.ID_ANY, "Disconnect" )
        self.Bind( wx.EVT_MENU, disconnect, id = item.GetId() )
        menu_ids.append( item.GetId() )
        menu.AppendSeparator()
        
        def connect_with_signature( signature ):
            connection.regions_attaches_front = signature
            
            PlaceInDepth( self.annotations.anns )
            self.GetTopLevelParent().Refresh()
        
        def connect_front_front( evt ):
            print('connect_front_front')
            connect_with_signature( (True,True) )
        def connect_back_back( evt ):
            print('connect_back_back')
            connect_with_signature( (False,False) )
        
        item = menu.Append( wx.ID_ANY, "Connect, Front" )
        self.Bind( wx.EVT_MENU, connect_front_front, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        item = menu.Append( wx.ID_ANY, "Connect, Front/Back" )
        self.Bind( wx.EVT_MENU, lambda evt: connect_with_signature( (True,False) ), id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        item = menu.Append( wx.ID_ANY, "Connect, Back/Front" )
        self.Bind( wx.EVT_MENU, lambda evt: connect_with_signature( (False,True) ), id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        item = menu.Append( wx.ID_ANY, "Connect, Back" )
        self.Bind( wx.EVT_MENU, connect_back_back, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        
        self.PopupMenu( menu )
        menu.Destroy()
        for id in menu_ids: self.Unbind( wx.EVT_MENU, id = id )
        
        return self
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
        ## in wxPython.
        if keycode in (wx.WXK_DELETE,8):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            try:
                ## TODO: Use proper oblique direction.
                shape = PickShape( self.annotations.anns, ( append( self.screen2world( screen_pos ), 0. ), (0,0,-1) ) )
                print('deleting a shape')
                DeleteSymmetryXSectionShape( self.annotations.anns, shape )
                
                PlaceInDepth( self.annotations.anns )
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        else:
            evt.Skip()
        
        return self
    
    def doDraw( self ):
        def draw_stroke():
            if self.stroke is None or len( self.stroke ) < 2: return
            
            glColor4ub( *( Color[ 'connection' ] + (200,) ) )
            glLineWidth( 8 )
            draw_linestrip( self.stroke )
        
        self.DrawScene( callback = draw_stroke )


class AnnotateCanvasSymmetryXSection( AnnotateCanvas ):
    def __init__( self, *args, **kwargs ):
        AnnotateCanvas.__init__( self, *args, **kwargs )
        
        self.RegisterMode( ObliqueMode, "Specify the oblique camera", "Oblique" )
        self.RegisterMode( AnnotateViewRotateMode, "Rotate the camera, temporarily", "Rotate" )
        #self.RegisterMode( AnnotationCreateMode, "Create shapes, silhouette first", "Create" )
        mode_id = self.RegisterMode( AnnotationCreateMode2, "Create shapes, no silhouette", "Create" )
        self.RegisterMode( ConnectionCurveDrawMode, "Draw connection curves", "Connect" )
        
        self.SetModeByID( mode_id )
    
    def DrawScene( self, **kwargs ):
        glClearColor( .9,.9,.9,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        ## Draw the gltex if there is one.
        if self.annotations.gltex is not None:
            self.prepare2D()
            
            glColor4f( 1,1,1, .5 )
            draw_gltexture( self.annotations.gltex, self.annotations.tex_pos, self.annotations.tex_size )
        
        if self.show_annotations:
            DrawAnnotations( self.prepare2D, self.prepare3D, self.annotations.anns, self )
        
        if 'callback' in kwargs: kwargs['callback']()
        
        self.SwapBuffers()

class ResultCanvasSymmetryXSection( ResultCanvas ):
    def doDraw( self ):
        glClearColor( .5,.5,.5,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        self.prepare3D()
        
        #self.drawTexture()
        ## Clear the depth buffer.
        ## Actually, this looks weird.
        # glClear( GL_DEPTH_BUFFER_BIT )
        
        
        ## Draw only the resulting meshes.
        for shape in self.annotations.anns.of_type( SymmetryXSectionShape ):
            DrawShape( self.prepare2D, self.prepare3D, shape, what_to_draw = 'mesh' )
        
        if wx_annotate_shell.kDebugging:
            DrawDebugging( self.prepare2D, self.prepare3D, self.annotations.anns )
        
        
        ## Draw the gltex if there is one.
        if self.annotations.gltex is not None:
            self.drawTexture(
                self.annotations.gltex,
                self.annotations.tex_pos,
                self.annotations.tex_size
                )
        
        self.SwapBuffers()

def main():
    import sys
    argv = list( sys.argv[1:] )
    parse_options( argv )
    
    import wx
    #app = wx.App() ## <-- on OS X tracebacks and stdout/stderr appear in a dialog
    app = wx.PySimpleApp()
    
    ## load the PyCrust shell
    ## NOTE: This generates the following errors and is very slow switching between the app window
    ##       and the PyCrust window.
    ## Sun Jun 15 20:19:36 yotronagon.local Python[32825] <Error>: CGContextRestoreGState: invalid context
    if False: #wx_annotate_shell.kDebugging:
        from wxBasicGLCanvas import load_PyCrust
        crust = load_PyCrust( app )
    
    annotations = Annotations()
    if wx_annotate_shell.kDebugging:
        shape = SymmetryXSectionShape()
        shape.circular_xsections = [(
            ( -.5, 0 ), ( .5, 0 ), ( 0, -.25 )
            )]
        shape.update_3d()
        shape.symmetry_xsection = (
            [ ( .2 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 20 )[:-1] + pi/2 ]
            )
        shape.update_3d()
        #from helpers import smooth_curve
        #shape.symmetry_xsection_3d = smooth_curve( shape.symmetry_xsection_3d, smooth = 11, simplify = False )
        annotations.anns.append( shape )
    
    frame = AnnotateFrame(
        None, annotations,
        annotate_canvas_factory = AnnotateCanvasSymmetryXSection,
        result_canvas_factory = ResultCanvasSymmetryXSection
        )
    frame.Show(True)
    
    for fname in argv:
        if fname.lower().endswith( '.annscene' ):
            annotations.Load( fname )
            frame.Refresh()
        elif fname.lower().endswith( '.png' ):
            annotations.SetImage( fname )
            frame.Refresh()
    
    app.annotateframe = frame
    app.MainLoop()
    app.Destroy()

if __name__ == '__main__': main()
