from illustrator_mode import *

class SweepShapeSilhouetteManip( StrokeDrawer ):
        def __init__( self, canvas, shape, t = None, right_turn = None ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.SaveStateForUndo()
            
            self.shape = shape
            self.spine_undeformed = array( self.shape.spine2d )
            self.dense_radii_undeformed = array( self.shape.dense_radii )
            
            self.t = t
            self.right_turn = right_turn
            
            if self.t is not None:
                
                assert self.t in self.shape.radii
                assert self.t >= 0.
                assert self.t <= 1.
                
                ## Right turn is a boolean.
                assert self.right_turn is not None
                assert bool( self.right_turn ) == self.right_turn
                
                self.vertex_index = self.shape.vertex_index_at_t( self.t )
                self.opposite_pos = self.shape.boundary_perpendicular_to_t( self.t )[ 1 if self.right_turn else 0 ]
        
        def doPointAdded( self ):
            ## Is this the first point?  If so, pick the closest spine point.
            if self.t is None:
                assert self.right_turn is None
                
                assert len( self.stroke ) == 1
                ## Cast a ray perpendicular to the closest point on the shape's boundary?
                offsets = annotations_sweep.boundary_offsets2d( self.shape.spine2d, self.shape.dense_radii )
                spine2d = asarray( self.shape.spine2d )
                ( distsqr, pt, edge_index ), right_turn = min(
                    ( closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
                        spine2d + offsets, self.stroke[0]
                        ),
                        True
                        ),
                    ( closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
                        spine2d - offsets, self.stroke[0]
                        ),
                        False
                        )
                    )
                ## Make sure the closest edge doesn't already have a cross section.
                ## TODO Q: Should I assign one regardless, and leave it up to the shape object
                ##         to tessellate itself sufficiently?
                if distsqr < EdgeDistanceThreshold( self )**2 and edge_index not in [ self.shape.edge_index_at_t( t ) for t in self.shape.radii.keys() ]:
                    self.t = annotations_sweep.edge_index2t( edge_index, self.shape.spine2d )
                    self.vertex_index = self.shape.vertex_index_at_t( self.t )
                    self.right_turn = right_turn
                    self.opposite_pos = self.shape.boundary_perpendicular_to_t( self.t )[ 1 if self.right_turn else 0 ]
                    ## Initialize the cross section to None.
                    ## (It will be set properly below.)
                    self.shape.radii[ self.t ] = None
                else:
                    unforward_events( self )
                    return
            
            assert self.opposite_pos is not None
            assert self.vertex_index is not None
            
            ## Update the circle's radius.
            mouse_pos = asarray( self.stroke[-1] )
            #print "mouse_pos:", mouse_pos
            #print "self.opposite_pos:", self.opposite_pos
            center = .5 * ( mouse_pos + self.opposite_pos )
            #print "center:", center
            radius = .5 * mag( mouse_pos - self.opposite_pos )
            self.shape.radii[ self.t ] = radius
            
            ## Update the spine.
            constraints = {}
            
            constrain_start_index = True
            tangent_constrain_start_index = False
            constrain_end_index = True
            tangent_constrain_end_index = False
            
            if self.vertex_index == 0:
                constrain_start_index = False
                tangent_constrain_start_index = False
            
            if self.vertex_index == len( self.shape.spine2d ) - 1:
                constrain_end_index = False
                tangent_constrain_end_index = False
            
            if constrain_start_index:
                constraints[ 0 ] = self.spine_undeformed[ 0 ]
            if tangent_constrain_start_index:
                constraints[ 1 ] = self.spine_undeformed[ 1 ]
            if constrain_end_index:
                constraints[ len( self.spine_undeformed ) - 1 ] = self.spine_undeformed[ -1 ]
            if tangent_constrain_end_index:
                constraints[ len( self.spine_undeformed ) - 2 ] = self.spine_undeformed[ -2 ]
            
            constraints[ self.vertex_index ] = center
            
            kConstrainEdge = True
            if kConstrainEdge:
                ## Constrain the edge at vertex_index so that self.opposite_pos stays
                ## in exactly the same place.
                edge_index = self.shape.edge_index_at_t( self.t )
                vertex_index2 = set([ edge_index, edge_index+1 ])
                vertex_index2.remove( self.vertex_index )
                assert 1 == len( vertex_index2 )
                vertex_index2 = list( vertex_index2 )[0]
                perp = dir( mouse_pos - self.opposite_pos )
                perp *= mag( self.spine_undeformed[ edge_index+1 ] - self.spine_undeformed[ edge_index ] )
                edge = zeros(2)
                ## a 90-degree right turn (CW rotation)
                edge[0] = perp[1]
                edge[1] = -perp[0]
                if self.right_turn: edge = -edge
                vertex_pos2 = center + (vertex_index2 - self.vertex_index) * edge
                constraints[ vertex_index2 ] = vertex_pos2
            
            from laplacian_curve_editing import new_curve_from_curve_and_constraints
            new_curve = new_curve_from_curve_and_constraints(
                self.spine_undeformed,
                curve_is_closed = False,
                constraints = list(constraints.items())
                )
            #print "|old spine - new spine|^2:", ( ( self.shape.spine2d - new_curve )**2 ).sum()
            self.shape.set_spine2d( new_curve )
            
            kSmoothRadii = False
            if kSmoothRadii:
                self.SolveForSmoothRadii1( radius )
                #self.SolveForSmoothRadii2( radius )
            
            self.Refresh()
        
        def SolveForSmoothRadii1( self, radius_at_vertex_index ):
            ## Minimize ( ( laplacian of radii ) - ( laplacian of undeformed radii ) )^2
            
            constraints = {}
            
            constrain_start_index = True
            tangent_constrain_start_index = False
            constrain_end_index = True
            tangent_constrain_end_index = False
            
            if self.vertex_index == 0:
                constrain_start_index = False
                tangent_constrain_start_index = False
            
            if self.vertex_index == len( self.shape.spine2d ) - 1:
                constrain_end_index = False
                tangent_constrain_end_index = False
            
            if constrain_start_index:
                constraints[ 0 ] = self.dense_radii_undeformed[ 0 ]
            if tangent_constrain_start_index:
                constraints[ 1 ] = self.dense_radii_undeformed[ 1 ]
            if constrain_end_index:
                constraints[ len( self.spine_undeformed ) - 1 ] = self.dense_radii_undeformed[ -1 ]
            if tangent_constrain_end_index:
                constraints[ len( self.spine_undeformed ) - 2 ] = self.dense_radii_undeformed[ -2 ]
            
            constraints[ self.vertex_index ] = radius_at_vertex_index
            
            radii_undeformed = asarray( [
                ( annotations_sweep.vertex_index2t( vi, self.spine_undeformed ), self.dense_radii_undeformed[ vi ] )
                for vi in range(len( self.dense_radii_undeformed ))
                ] )
            print("radii_undeformed:")
            print(radii_undeformed)
            from laplacian_curve_editing import L_for_curve
            L = L_for_curve( radii_undeformed, curve_is_closed = False )
            L2 = dot( L, L )
            rhs = dot( L2, radii_undeformed )
            print("rhs:")
            print(rhs)
            print("constraints:")
            print(constraints)
            for vi, val in constraints.items():
                
                L2[ vi, : ] = 0.
                L2[ vi, vi ] = 1.
                rhs[ vi ] = val
            
            new_radii = linalg.solve( L2, rhs )[:,1]
            print("new_radii:")
            print(new_radii)
            self.shape.radii.clear()
            for t, radius in zip( radii_undeformed[:,0], new_radii ):
                self.shape.radii[ t ] = radius
            
            '''
            new_radii = new_curve_from_curve_and_constraints(
                radii_undeformed,
                curve_is_closed = False,
                constraints = constraints.items()
                )
            print "new_radii:"
            print new_radii
            self.shape.radii.clear()
            for t,radius in new_radii:
                self.shape.radii[ t ] = radius
            '''
        
        def SolveForSmoothRadii2( self, radius_at_vertex_index ):
            ## Our energy is
            ##   ( laplacian magnitudes of radii )^2
            ##   + ( boundary2d_right_deformed - boundary_2d_right_undeformed )^2
            ##   + ( boundary2d_left_deformed - boundary_2d_left_undeformed )^2
            ## Where boundary2d_right_deformed is calculated as: spine2d_deformed + radii * unit_offsets_deformed
            ## TODO Q: Don't add energy term for self.right_turn boundary at self.vertex_index?
            ## A1: I think that's a good idea.
            ## TODO Q: Why not just minimize ( ( laplacian of radii ) - ( laplacian of undeformed radii ) )^2?
            
            
            
            ## The next six lines are from annotations_sweep.boundary_offsets2d().
            normals_deformed = tangents_and_normals_along_line_strip_2d( self.shape.spine2d )[1]
            normals_undeformed = tangents_and_normals_along_line_strip_2d( self.spine_undeformed )[1]
            ## Repeat the last normal, because there is one less edge than point.
            unit_offsets_deformed = append( normals_deformed, normals_deformed[-1:], 0 )
            unit_offsets_undeformed = append( normals_undeformed, normals_undeformed[-1:], 0 )
            
            ## annotations_sweep.gen_offset_boundary2d() uses positive offsets for right
            ## and negative offsets for left.
            boundary2d_right_deformed = self.shape.spine2d + self.shape.dense_radii[:,newaxis] * unit_offsets_deformed
            boundary2d_left_deformed = self.shape.spine2d - self.shape.dense_radii[:,newaxis] * unit_offsets_deformed
            
            boundary2d_right_undeformed = self.spine_undeformed + self.dense_radii_undeformed[:,newaxis] * unit_offsets_undeformed
            boundary2d_left_undeformed = self.spine_undeformed - self.dense_radii_undeformed[:,newaxis] * unit_offsets_undeformed
        
        def doStrokeFinished( self ):
            ## Call smooth_spine2d(), because it resamples the spine evenly.
            ## UPDATE: It makes things jump around too much.
            #self.shape.smooth_spine2d( 0 )
            self.shape.update_3d()
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            self.shape = None
        
        def doDraw( self ):
            ## We override doDraw() so that we don't draw the stroke.  Instead, draw the circle.
            ## Actually, adding the tilted circle to self.shape.tilted_circles means that the
            ## regular shape drawing code draws the cross section.
            self.DrawSweepShapeHandles( [ self.shape ] )
