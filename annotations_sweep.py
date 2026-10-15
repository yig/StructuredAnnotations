from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *
from helpers import *
from vecutil import *
kEps = 1e-5

#kMinRadius = .04
kMinRadius = .01
kMaxAbsTiltRadians = .5*pi * .9
kNoPinch = True

## Number of samples for arbitrary cross sections.
kCrossSectionSamples = 21

unit_circle = asarray( [ ( cos( theta ), sin( theta ) ) for theta in linspace( 0, 2*pi, kCrossSectionSamples+1 )[:-1] ] )

def boundary_offsets2d( spine2d, radii ):
    '''
    Given a sequence 'spine2d' of 2d points for a spine and
    a sequence 'radii' of scalar radii for each point,
    returns the offset vectors 'offsets' for each point in spine2d
    such that 'spine2d' +/- 'offsets' are the offset curves by radii
    in the right-turn and left-turn curve normal directions.
    '''
    
    assert len( spine2d ) == len( radii )
    assert len( spine2d[0] ) == 2
    
    spine2d = asarray( spine2d )
    radii = asarray( radii )
    
    unit_edges, normals = tangents_and_normals_along_line_strip_2d( spine2d )
    ## Repeat the last normal, because there is one less edge than point.
    offsets = append( normals, normals[-1:], 0 )
    offsets *= radii[:,newaxis]
    
    return offsets

def gen_offset_boundary2d( spine2d, radii, caps ):
    '''
    Given a sequence 'spine2d' of 2d points for a spine and
    a sequence 'radii' of scalar radii for each point,
    returns a line loop for the boundary of the spine offset
    by the radii in the +/- normal directions.
    '''
    
    assert len( spine2d ) == len( radii )
    assert len( spine2d[0] ) == 2
    
    import sweep_surface
    spine2d, radii, start_pt, end_pt = sweep_surface.path_and_radii_with_caps(
        spine2d, radii, caps = caps, cap_segments = 10
        )
    
    spine2d = asarray( spine2d )
    radii = asarray( radii )
    
    ## Walk the boundary by walking along the spine offset to the right and left by radius.
    offsets = boundary_offsets2d( spine2d, radii )
    
    ## TODO Q: If kNoPinch, should I create the boundary using the smoothed distance field?
    ## A1: No, this will cause problems with boundary picking code.
    ## A2: Yes.  Boundary picking code just looks for the closest spine point, so it won't
    ##     be much affected.
    right = spine2d + offsets
    left = (spine2d - offsets)[::-1]
    
    boundary_loop = list(right) + [end_pt] + list(left) + [start_pt]
    return boundary_loop

def t2edge_index( t, line_strip ):
    '''
    Returns the index of the edge in line_strip that 't' lies in.
    '''
    
    ## Clip to [0,1]
    t = min( 1., max( 0., t ) )
    
    N = len( line_strip )
    ei = min( int(floor( t * ( N-1 ) )), N-2 )
    assert ei >= 0 and ei < N-1
    return ei

def edge_index2t( edge_index, line_strip ):
    '''
    Returns t along the line_strip corresponding to the given edge index.
    '''
    
    ## This maps the first edge to 0 and the last edge to 1.
    return float( edge_index ) / ( len( line_strip )-2 )
    ## This maps each edge to the t halfway along it.
    #return ( edge_index + .5 ) / ( len( line_strip )-1 )

def vertex_index2s( vertex_index, line_strip ):
    '''
    Returns s along the line_strip corresponding to the given vertex index.
    
    NOTE: s is an arc length parameterization of the curve, whereas t is the existing (vertex)
          parameterization of the curve.
    '''
    
    ## Get the length of each edge.  ('seg_lengths' has size #vertices - 1.)
    seg_lengths = line_strip_seg_lengths( line_strip )
    
    ## Construct the edge-to-s map, normalizing values from 0 to 1.
    edge2s = seg_lengths.cumsum()
    edge2s /= edge2s[-1]
    
    if vertex_index == 0:
        return 0.
    else:
        return edge2s[ vertex_index - 1 ]

def t2vertex_index( t, line_strip ):
    '''
    Returns the index of the vertex in line_strip that 't' lies in.
    '''
    
    ## Clip to [0,1]
    t = min( 1., max( 0., t ) )
    
    N = len( line_strip )
    vi = int(round( t * ( N-1 ) ))
    assert vi >= 0 and vi < N
    return vi

def vertex_index2t( vertex_index, line_strip ):
    '''
    Returns t along the line_strip corresponding to the given vertex index.
    '''
    
    ## This maps the first vertex to 0 and the last vertex to 1.
    return float( vertex_index ) / ( len( line_strip )-1 )

def tilt_radians2tangent3d( tilt_radians, tangent2d ):
    ## Clamp tilt_radians.
    tilt_radians = max( -kMaxAbsTiltRadians, min( kMaxAbsTiltRadians, tilt_radians ) )
    
    tangent3d = zerosf( 3 )
    tangent3d[:2] = tangent2d
    tangent3d[2] = tan( -tilt_radians ) * mag( tangent2d )
    
    return tangent3d

def normalize_arbitrary_cross_section( xsection ):
    '''
    Given a sequence of kCrossSectionSamples 2d points, returns the same sequence
    normalized to be roughly the same size as a unit circle (currently max distance is 1,
    may be changed to average distance).
    Also shifts the points so the first point is the one closest to the x-axis.
    Reverses the points if they are oriented clock-wise (CW) instead of counter-clockwise (CCW).
    '''
    
    assert len( xsection ) == kCrossSectionSamples
    
    xsection = asarray( xsection )
    
    ## Translate so that it's centered.
    xsection -= average( xsection, axis = 0 )
    ## Scale so that its max distance from 0 is 1.
    ## TODO Q: Scale so that its maximum or average distance is 1?
    ## TODO: Scale so that its maximum in the plane is 1.  This would make boundaries computed
    ##       with boundary_offsets2d() correct.  Unfortunately it requires knowledge of the
    ##       symmetry plane direction (which changes, too).
    ##       (Given the symmetry plane direction, rotate the cross section so that its x direction
    ##       is aligned with the symmetry plane, then normalize so that the maximum |x| coordinate
    ##       is 1.  This is still inaccurate, because it doesn't center the average x at 0 and
    ##       because the symmetry plane can change along the spine.)
    mags = sqrt( ( xsection**2 ).sum(1) )
    xsection *= 1./max( mags )
    
    ## TODO: Symmetrize the cross section.
    
    ## Shift points so that the first point is the one closest to the x-axis.
    ## 1 Find adjacent, successive (CCW) points which cross the x axis (y < 0, y >= 0).
    ## 2 Cycle the minimum amount so that the -1-th and 0-th point cross the x-axis.
    
    ## We'll store edge indices.
    crossings = []
    N = len( xsection )
    for ei in range( N ):
        if xsection[ei][1] < 0 and xsection[ (ei+1) % N ][1] >= 0:
            crossings.append( ei )
    
    ## There must be at least one crossings, because the loop has its center-of-mass at the origin.
    best_crossing = min( [ ( min( ei + 1, N-1 - ei ), ei ) for ei in crossings ] )[1]
    ## 'best_crossing' is the index of the edge which we want to straddle the end-to-beginning
    ## of the loop.
    ## This means we want the loop to start with the vertex at best_crossing+1 and end
    ## with the vertex at best_crossing.
    xsection = append( xsection[ best_crossing+1 : ], xsection[ : best_crossing+1 ], 0 )
    
    if not line_loop_is_CCW( xsection ):
        print("Flipping clock-wise loop.")
        xsection = xsection[::-1]
    
    return xsection

class SweepShape( object ):
    def __init__( self, spine2d = None, tilt_radians = None, radii = None, cross_sections = None, symmetry_radians = None, cap_radii = None ):
        '''
        Given a simple, open path (sequence) of 2d points for the spine 'spine2d',
        and dictionaries mapping percents along the path to attributes:
            'tilt_radians', tilt in radians about the axis perpendicular (right turn) to spine2d
                            (0 corresponds to a circle with one of its axes the z axis; in other
                            words, its normal is the 2d edge vector)
            'radii', the radius of the cross section,
            'cross_sections', a 2d line loop that will be scaled by the shape's radius at t,
                              oriented such that its x axis points along the symmetry direction
                              and its z axis is tangent to the 3d spine,
            'symmetry_radians', the angle of the symmetry plane (sheet, really) in radians
                                (0 corresponds to the z axis)
        and an optional tuple of cap distances 'cap_radii' at the start and end of the spine,
        creates a SweepShape object.
        '''
        
        ## NOTE: I can't make {} a default parameter, because this always uses the same instance,
        ##       so the instance is shared by all objects instantiated with the default parameter,
        ##       and when I append to that instance, it is seen by all objects sharing it.
        ## These are sensible defaults:
        #if tilt_radians is None: tilt_radians = { 0: 0., 1: 0. }
        #if radii is None: radii = { 0: [ .1 ], 1: [ .1 ] }
        #if cross_sections is None: cross_sections = { .5: unit_circle }
        #if symmetry_radians is None: symmetry_radians = { .5: -.5*pi }
        if tilt_radians is None: tilt_radians = {}
        if radii is None: radii = {}
        if cross_sections is None: cross_sections = {}
        if symmetry_radians is None: symmetry_radians = {}
        if cap_radii is None: cap_radii = [None,None]
        
        self.spine2d = spine2d
        ## Copy elements so we can modify them.
        ## The tilt radians can be turned into a slope, which we normalize to get the cross
        ## section's normal.
        ## Viewed in a uv plane, where the plane is defined as perpendicular to
        ## (and viewed along this vector, negative) the right turn normal to the edge
        ## (this normal is actually in the xy plane) and having as its horizontal +u axis
        ## the 2d edge and its vertical +v axis the z direction,
        ## we see that the v-axis value (the world-space z value) is the tangent of
        ## -tilt_radians, scaled by the 2d edge length.
        ## For a diagram, see "discussion/SweepShape-tilt_radians.pdf".
        self.tilt_radians = dict( tilt_radians )
        self.radii = dict( radii )
        self.cross_sections = dict( [
            ( t, asarrayf( xsection ) ) for t, xsection in cross_sections.items()
            ] )
        ## The symmetry plane is the sheet (2-manifold) defined as the union of all lines
        ## passing through self.spine3d in the direction determined by self.symmetry_plane_radians
        ## as follows:
        ##   Let x be the -look vector, let y be the 2d right-turn spine2d normal vector with z = 0,
        ##   then the symmetry plane vector at a spine3d point spine3d[i] is
        ##   cos( symmetry_plane_radians ) * -look + sin( symmetry_plane_radians ) * normal[i]
        self.symmetry_radians = dict( symmetry_radians )
        self.cap_radii = list( cap_radii )
        
        import uuid
        self.uuid = uuid.uuid4()
        
        ## Generated attributes.
        ## TODO: Change these to properties and generate them on-the-fly as needed. (Reset them when update_3d() is called.)
        self.everything_changed()
        
        self.assert_self()
    
    def everything_changed( self ):
        ## Reset generated attributes to None.
        self._mesh = None
        self._spine3d = None
        self._spine3d_diffuser = None
        self._dense_radii = None
        self._symmetry_directions = None
        self._boundary2d = None
        
        self._mesh_hash = None
        self._spine2d_hash = None
        self._spine3d_hash = None
        self._spine3d_diffuser_hash = None
        self._dense_radii_hash = None
        self._dense_cross_sections_hash = None
        self._symmetry_directions_hash = None
        self._boundary2d_hash = None
    
    def tilt_radians_hash( self ):
        return hash( tuple( self.tilt_radians.items() ) )
    def radii_hash( self ):
        return hash( tuple( self.radii.items() ) )
    def cross_sections_hash( self ):
        #return hash( tuple( [ ( t, tuple( map( tuple, xsection ) ) ) for t, xsection in self.cross_sections.iteritems() ] ) )
        return hash( tuple( [ ( t, asarray( xsection ).tobytes() ) for t, xsection in self.cross_sections.items() ] ) )
    def symmetry_radians_hash( self ):
        return hash( tuple( self.symmetry_radians.items() ) )
    def cap_radii_hash( self ):
        return hash( tuple( self.cap_radii ) )
    
    def spine2d_hash( self ):
        if self._spine2d_hash is None and self.spine2d is not None:
            #self._spine2d_hash = hash( tuple( map( tuple, self.spine2d ) ) )
            self._spine2d_hash = hash( asarray( self.spine2d ).tobytes() )
        
        return ( self._spine2d_hash, )
    
    def __setstate__( self, unpickled ):
    
        ## Convert older-version objects.
        if 'tilted_circles' in unpickled:
            
            ## Check for *really* old classes, which didn't have symmetry_plane_radians.
            if not 'symmetry_plane_radians' in unpickled:
                if 'symmetry_plane_angle' in unpickled:
                    unpickled['symmetry_plane_radians'] = unpickled['symmetry_plane_angle']
                    del unpickled['symmetry_plane_angle']
                else:
                    unpickled['symmetry_plane_radians'] = -.5*pi
            
            self.everything_changed()
            
            self.spine2d = unpickled['spine2d']
            self.symmetry_radians = {}
            self.radii = {}
            self.cross_sections = {}
            self.tilt_radians = {}
            self.cap_radii = [None, None]
            
            import uuid
            self.uuid = uuid.uuid4()
            
            ## This is ugly but it won't be saved and it's necessary when unpickling annotations.
            self._tilted_circle_index2t = {}
            
            ## The handle used to be located 1/3 of the way along the spine.
            self.symmetry_radians[ 1./3. ] = unpickled['symmetry_plane_radians']
            for i, ( edge_index, radius, tilt_radians ) in enumerate( unpickled['tilted_circles'] ):
                t = edge_index2t( edge_index, unpickled['spine2d'] )
                self.radii[ t ] = radius
                self.tilt_radians[ t ] = tilt_radians
                
                self._tilted_circle_index2t[ i ] = t
            
        else:
            
            ## I recently added self.cap_radii.
            if 'cap_radii' not in unpickled:
                unpickled['cap_radii'] = [None,None]
            
            if '_spine2d_hash' not in unpickled:
                unpickled['_spine2d_hash'] = None
            
            if '_dense_cross_sections_hash' not in unpickled:
                unpickled['_dense_cross_sections_hash'] = None
            
            if 'uuid' not in unpickled:
                import uuid
                unpickled['uuid'] = uuid.uuid4()
            
            self.__dict__.update( unpickled )
    
    def __getstate__( self ):
        '''
        This method is called by the pickler, and we use it to not pickle generated properties.
        '''
        
        state = dict( self.__dict__ )
        for key in list(state.keys()):
            if key[0] == '_':
                state[key] = None
        return state
    
    def fast__deepcopy__( self, memodict ):
        result = SweepShape()
        
        result.spine2d = array( self.spine2d )
        result.tilt_radians = dict( self.tilt_radians )
        result.radii = dict( self.radii )
        result.cross_sections = dict( [ ( t, array( xsection ) ) for t, xsection in self.cross_sections.items() ] )
        result.symmetry_radians = dict( self.symmetry_radians )
        result.cap_radii = list( self.cap_radii )
        result.uuid = self.uuid
        
        result.everything_changed()
        
        return result
    
    def assert_self( self ):
        ## All attributes must be attached to the curve, which is defined from 0 to 1.
        assert False not in [ 0 <= t <= 1 for t in self.tilt_radians.keys() ]
        assert False not in [ 0 <= t <= 1 for t in self.radii.keys() ]
        assert False not in [ 0 <= t <= 1 for t in self.cross_sections.keys() ]
        assert False not in [ 0 <= t <= 1 for t in self.symmetry_radians.keys() ]
        
        ## All specified cross sections must have the same number of samples.
        assert len( set( [ len( xsection ) for xsection in self.cross_sections.values() ] ) ) <= 1
        ## All radii must be positive.
        assert False not in [ radius > 0 for radius in self.radii.values() ]
        
        ## Tilt radians should never be larger than pi/2.
        ## NOTE: Clamping to kMaxAbsTiltRadians happens when solving, not here.
        assert False not in [ abs( radians ) < .5*pi + kEps for radians in self.tilt_radians.values() ]
        ## The symmetry plane angle should also never be larger that pi/2.
        ## TODO Q: Or never larger than pi?
        #assert False not in [ abs( radians ) < .5*pi + kEps for radians in self.symmetry_radians.itervalues() ]
        assert False not in [ abs( radians ) < pi + kEps for radians in self.symmetry_radians.values() ]
        
        ## There should be two cap radii, both floating point numbers.
        assert len( self.cap_radii ) == 2
        assert self.cap_radii[0] is None or float( self.cap_radii[0] ) == self.cap_radii[0]
        assert self.cap_radii[1] is None or float( self.cap_radii[1] ) == self.cap_radii[1]
    
    def spine2d_at_t( self, t ):
        '''
        Returns the 2d point on the spine parameterized by 0 <= t <= 1.
        '''
        
        ## For now, piecewise constant.
        return self.spine2d[ self.vertex_index_at_t( t ) ]
    
    def spine3d_at_t( self, t ):
        '''
        Returns the 3d point on the spine parameterized by 0 <= t <= 1.
        '''
        
        ## For now, piecewise constant.
        return self.spine3d[ self.vertex_index_at_t( t ) ]
    
    def tangent3d_at_t( self, t ):
        '''
        Returns the 3d tangent vector on the spine parameterized by 0 <= t <= 1.
        '''
        
        ## Always return the specifically set value if present.
        ## TODO Q: Should I have a different interface for this?
        if t in self.tilt_radians:
            return tilt_radians2tangent3d( self.tilt_radians[ t ], self.tangent2d_at_t( t ) )
        
        ei = self.edge_index_at_t( t )
        return self.spine3d[ ei+1 ] - self.spine3d[ ei ]
    
    def tangent2d_at_t( self, t ):
        '''
        Returns the 2d tangent vector on the spine parameterized by 0 <= t <= 1.
        '''
        
        ei = self.edge_index_at_t( t )
        return asarray( self.spine2d[ ei+1 ] ) - asarray( self.spine2d[ ei ] )
    
    def normal2d_at_t( self, t ):
        '''
        Returns the 2d normal vector (90-degree right turn) to the spine
        parameterized by 0 <= t <= 1.
        '''
        
        tangent2d = self.tangent2d_at_t( t )
        return asarray( ( tangent2d[1], -tangent2d[0] ) )
    
    def radius_at_t( self, t ):
        '''
        Returns the radius of the cross section at t.
        '''
        
        ## Always return the specifically set value if present.
        ## TODO Q: Should I have a different interface for this?
        if t in self.radii: return self.radii[ t ]
        
        ## For now, piecewise constant.
        return self.dense_radii[ self.vertex_index_at_t( t ) ]
    
    def cross_section_at_t( self, t ):
        '''
        Returns the normalized cross section at t as a 2d line loop.
        '''
        
        ## Always return the specifically set value if present.
        ## TODO Q: Should I have a different interface for this?
        if t in self.cross_sections: return self.cross_sections[ t ]
        
        if len( self.dense_cross_sections ) > 0:
            ## For now, piecewise constant.
            return self.dense_cross_sections[ self.vertex_index_at_t( t ) ]
        else:
            return unit_circle
    
    def symmetry_radians_at_t( self, t ):
        '''
        Returns the symmetry radians at t as a scalar.
        '''
        
        ## Always return the specifically set value if present.
        ## TODO Q: Should I have a different interface for this?
        if t in self.symmetry_radians: return self.symmetry_radians[ t ]
        
        vertex_index = self.vertex_index_at_t( t )
        
        ups, normals, tangents = self.get_frame_components()
        assert len( ups ) == len( self.spine2d )
        assert len( normals ) == len( self.spine2d )
        assert len( tangents ) == len( self.spine2d )
        
        sym_dir = self.symmetry_directions_at_vertices()[ vertex_index ]
        sym_dir = dir( sym_dir )
        
        ## In gen_symmetry_directions(), ups are used as the x directions,
        ## and normals are used as the y directions.
        return atan2( dot( sym_dir, normals[ vertex_index ] ), dot( sym_dir, ups[ vertex_index ] ) )
    
    def cap_radius_at_t( self, t ):
        '''
        Returns the cap radius at t.
        NOTE: t must be either 0 or 1, since the cap radii only exist at the spine endpoints.
        '''
        
        assert t == 0 or t == 1
        result = self.cap_radii[ int(t) ]
        if result is None:
            result = self.radius_at_t( t )
        
        return result
    
    def edge_index_at_t( self, t ):
        '''
        Returns the index of the edge in self.spine2d/3d that 't' lies in.
        '''
        
        return t2edge_index( t, self.spine2d )
    
    def vertex_index_at_t( self, t ):
        '''
        Returns the index of the vertex in self.spine2d/3d that 't' is closest to.
        '''
        
        return t2vertex_index( t, self.spine2d )
    
    def boundary_perpendicular_to_t( self, t ):
        '''
        Returns two 2d points, the boundary of the surface in the perpendicular directions of the
        spine.  The first point is a 90-degree right turn from the tangent,
        and the second a 90-degree left turn.
        '''
        
        point = asarray( self.spine2d_at_t( t ) )
        radius = self.radius_at_t( t )
        
        tangent = self.tangent2d_at_t( t )
        perp = zerosf(2)
        perp[0] = tangent[1]
        perp[1] = -tangent[0]
        perp = dir( perp )
        
        return point + radius*perp, point - radius*perp
    
    def tilt_ray_at_t( self, t ):
        '''
        Returns a 3d ray originating from the point t along the spine in the tilted direction
        (the ray points along z if the spine is not tilted at t).
        '''
        
        ray_origin = self.spine3d_at_t( t )
        
        ## The 3d spine tangent is 90-degrees off from the one we want.
        tangent3d = self.tangent3d_at_t( t )
        axis = cross( tangent3d, (0,0,1) )
        ## These two lines compute the same thing:
        #ray_direction = rotate( tangent3d, axis, pi/2 )
        ray_direction = cross( axis, tangent3d )
        
        return ray_origin, ray_direction
    
    def compute_tilt_radians_from_point_at_t( self, t, point ):
        '''
        Given a parameter along the spine t and a 2d point, returns tilt radians such that
        tilting ( spine2d( t ) + (0,0,-radius) ) by the result is as close as possible to 'point'.
        '''
        
        center = asarray( self.spine2d_at_t( t ) )
        unit_tangent2d = dir( self.tangent2d_at_t( t ) )
        radius = self.radius_at_t( t )
        x = dot( asarrayf( point ) - center, unit_tangent2d ) / radius
        ## Following the comments in __init__() and "discussion/SweepShape-tilt_radians.pdf":
        tilt_radians = asin( max( -1., min( 1., x ) ) )
        return tilt_radians
    
    def cross_section_line_loop3d_at_t( self, t ):
        '''
        Returns a 3d line loop representing the cross section at t along the spine.
        
        NOTE: t must appear in self.cross_sections.
        '''
        
        if len( self.dense_cross_sections ) == 0:
            center = self.spine3d_at_t( t )
            normal = self.tangent3d_at_t( t )
            radius = self.radius_at_t( t )
            
            line_loop = line_loop_from_tilted_circle( center, radius, normal )
            return line_loop
        
        else:
            line_loop = self.cross_section_at_t( t )
            radius = self.radius_at_t( t )
            
            assert len( line_loop ) > 2
            assert line_loop.shape[1] == 2
            
            ## I could use the origin of the ray self.symmetry_plane_ray_at_t(), but I'm afraid it
            ## won't be sampled the same way as self.spine3d_at_t().
            center = self.spine3d_at_t( t )
            xvec = dir( self.symmetry_plane_ray_at_t( t )[1] )
            zvec = dir( self.tangent3d_at_t( t ) )
            yvec = cross( zvec, xvec )
            
            ## Now transform the line loop from a 2d line loop to one whose axes are xvec, yvec, zvec,
            ## and scale by radius.
            line_loop = concatenate( ( line_loop, zeros( ( len(line_loop), 1 ) ) ), axis=1 )
            ## TODO: assert that this rotation matrix is orthonormal.
            rot = vstack( ( xvec, yvec, zvec ) ).T
            line_loop3d = center + radius * dot( rot, line_loop.T ).T
            return line_loop3d
    
    def symmetry_plane_ray_at_t( self, t ):
        '''
        Returns the symmetry plane ray at t along the spine.
        '''
        vi = self.vertex_index_at_t( t )
        return self.symmetry_plane_rays_at_vertices()[ vi ]
    
    def tilt_radians_at_t( self, t ):
        '''
        Given t along the spine, returns the tilted radians.
        
        untested
        '''
        
        ## Following the comments in SweepShape.__init__() and
        ## "discussion/SweepShape-tilt_radians.pdf":
        normal = self.tangent3d_at_t( t )
        #tilt_radians = asin( max( -1., min( 1., dir(normal)[2] ) ) )
        tilt_radians = -atan( normal[2] / mag(normal[:2]) )
        print('atan vs asin:', tilt_radians - asin( max( -1., min( 1., dir(normal)[2] ) ) ))
        
        return tilt_radians
    
    def symmetry_directions_at_vertices( self ):
        return self.symmetry_directions
    
    def get_caps( self ):
        #kEndCap = ('spherical',)
        #return kEndCap,kEndCap
        
        caps = [
            [ 'fixed', self.cap_radii[0] ],
            [ 'fixed', self.cap_radii[1] ]
            ]
        
        ## If unspecified, use end-caps the same radius as at the end of the spine.
        #if caps[0][1] is None: caps[0][1] = self.dense_radii[0]
        #if caps[-1][1] is None: caps[-1][1] = self.dense_radii[-1]
        
        ## If unspecified, use spherical end-caps.
        if caps[0][1] is None: caps[0] = ['spherical']
        if caps[-1][1] is None: caps[-1] = ['spherical']
        
        #print "caps:", caps
        return caps
    
    def set_radii_from_silhouette( self, silhouette ):
        '''
        Given a silhouette curve (sequence of 2d points forming a line loop),
        sets the radii of tilted circles along the spine.
        
        used
        '''
        
        ## TODO: take a max 2d arc length parameter that inserts a new circle between two others
        ##       (using tilt_radians_at_t()) if the distance is greater than the
        ##       parameter.
        
        silhouette = list(silhouette)
        
        if line_strips_intersect( silhouette + [silhouette[0]], self.spine2d ):
            print("Refusing to set radii from silhouette: silhouette intersects spine.")
            return
        
        if not contour_contains( silhouette, self.spine2d[0] ):
            print("Refusing to set radii from silhouette: silhouette does not enclose spine.")
            return
        
        for t in self.radii.keys():
            radius = self.radii[ t ]
            
            point = self.spine2d_at_t( t )
            unit_tangent2d = dir( self.tangent2d_at_t( t ) )
            perp_ray = ( point, ( unit_tangent2d[1], -unit_tangent2d[0] ) )
            
            isects = ray_lineloop_intersection_ts( perp_ray, silhouette, ray_is_line = True )
            abs_isect_pos_t = [ abs( t_etc[0] ) for t_etc in isects if t_etc[0] >= 0. ]
            abs_isect_neg_t = [ abs( t_etc[0] ) for t_etc in isects if t_etc[0] < 0. ]
            assert len( abs_isect_pos_t ) > 0 and len( abs_isect_neg_t ) > 0
            
            ## Set the radius to the average of the two closest +- intersection points.
            self.radii[ t ] = .5 * ( min( abs_isect_pos_t ) + min( abs_isect_neg_t ) )
    
    def set_spine_and_radii_from_silhouette( self, silhouette ):
        '''
        Given a silhouette curve (sequence of 2d points forming a line loop),
        sets the 2d spine positions and radii of tilted circles from the silhouette.
        
        used
        '''
        
        ## TODO: take a max 2d arc length parameter that inserts a new circle between two others
        ##       (using tilt_radians_at_t()) if the distance is greater than the
        ##       parameter.
        
        silhouette = list(silhouette)
        self.debug_draw = [ ( 'lineloop', silhouette ) ]
        
        if line_strips_intersect( silhouette + [silhouette[0]], self.spine2d ):
            print("Refusing to set radii from silhouette: silhouette intersects spine.")
            return
        
        if not contour_contains( silhouette, self.spine2d[0] ):
            print("Refusing to set radii from silhouette: silhouette does not enclose spine.")
            return
        
        
        ## Get a list of normal (perpendicular) directions for each point.
        unit_edges, normals = tangents_and_normals_along_line_strip_2d( self.spine2d )
        ## Repeat the last normal, because there is one less edge than point.
        normals = append( normals, normals[-1:], 0 )
        
        new_spine = self.spine2d
        
        ## Update positions, moving them in the normal direction.
        radii = emptyf( len( new_spine ) )
        for i, pt in enumerate( new_spine ):
            normal = normals[i]
            
            isects = ray_lineloop_intersection_ts( ( pt, normal ), silhouette, ray_is_line = True )
            abs_isect_pos_t = [ abs( t_etc[0] ) for t_etc in isects if t_etc[0] >= 0. ]
            abs_isect_neg_t = [ abs( t_etc[0] ) for t_etc in isects if t_etc[0] < 0. ]
            assert len( abs_isect_pos_t ) > 0 and len( abs_isect_neg_t ) > 0
            
            pt += .5 * ( min( abs_isect_pos_t ) - min( abs_isect_neg_t ) ) * normal
            
            radii[i] = .5 * ( min( abs_isect_pos_t ) + min( abs_isect_neg_t ) )
        
        
        ## Update circle radii.
        ## NOTE: We must do this before we smooth/resample the spine, because
        ##       the 'radii' array is tied to the current spine.
        for t in self.radii.keys():
            ## Set the radius to the average of the two closest +- intersection points.
            self.radii[ t ] = radii[ self.edge_index_at_t( t ) ]
        
        
        ## Q: Smooth spine2d?
        ## A: Yes, and resample.  Not resampling creates a lot of jitter in the outline.
        self.set_spine2d( smooth_open_stroke( new_spine, smooth = 4, resample_N = len( new_spine ) ) )
    
    def set_spine_and_radii_from_silhouette( self, silhouette, extra = 0 ):
        '''
        Given a silhouette curve (sequence of 2d points forming a line loop),
        sets the 2d spine positions and radii of tilted circles from the silhouette.
        
        used
        '''
        
        ## TODO: take a max 2d arc length parameter that inserts a new circle between two others
        ##       (using tilt_radians_at_t()) if the distance is greater than the
        ##       parameter.
        
        silhouette = list(silhouette)
        self.debug_draw = [ ( 'lineloop', silhouette ) ]
        
        
        ## These checks have been moved deeper.
        '''
        if line_strips_intersect( silhouette + [silhouette[0]], self.spine2d ):
            print "Refusing to set radii from silhouette: silhouette intersects spine."
            return
        
        if not contour_contains( silhouette, self.spine2d[0] ):
            print "Refusing to set radii from silhouette: silhouette does not enclose spine."
            return
        '''
        
        
        ## Get a list of normal (perpendicular) directions for each point.
        unit_edges, normals = tangents_and_normals_along_line_strip_2d( self.spine2d )
        ## Repeat the last normal, because there is one less edge than point.
        normals = append( normals, normals[-1:], 0 )
        
        
        ## Transform the spine so that its start and end tangent directions
        ## are close to the silhouette.
        new_start_isects = ray_lineloop_intersection_ts( ( self.spine2d[0], -unit_edges[0] ), silhouette, ray_is_line = False )
        if len( new_start_isects ) == 0:
            print("Extension of spine from start doesn't hit silhouette.")
            return
        assert all([ t_etc[0] >= 0. for t_etc in new_start_isects ])
        #new_start = self.spine2d[0] + ( min( new_start_isects )[0] - self.radius_at_t(0) ) * -unit_edges[0]
        new_start = self.spine2d[0] + ( min( new_start_isects )[0] - .1 ) * -unit_edges[0]
        new_start_cap_radius = .1
        
        new_end_isects = ray_lineloop_intersection_ts( ( self.spine2d[-1], unit_edges[-1] ), silhouette, ray_is_line = False )
        if len( new_end_isects ) == 0:
            print("Extension of spine from end doesn't hit silhouette.")
            return
        assert all([ t_etc[0] >= 0. for t_etc in new_start_isects ])
        #new_end = self.spine2d[-1] + ( min( new_end_isects )[0] - self.radius_at_t(1) ) * unit_edges[-1]
        new_end = self.spine2d[-1] + ( min( new_end_isects )[0] - .1 ) * unit_edges[-1]
        new_end_cap_radius = .1
        
        new_spine = solve_for_fixed_points(
            self.spine2d,
            ( self.spine2d[0], self.spine2d[-1] ),
            ( new_start, new_end )
            )
        self.debug_draw.append( ( 'linestrip', new_spine, ( 0, 255, 0 ) ) )
        if ( extra % 10 ) != 0:
        #if True:
        #if False:
            new_spine = array( self.spine2d )
            new_start_cap_radius = self.cap_radii[ 0 ]
            new_end_cap_radius = self.cap_radii[ 1 ]
        
        ## Updated tangents, normals
        unit_edges, normals = tangents_and_normals_along_line_strip_2d( new_spine )
        ## Repeat the last normal, because there is one less edge than point.
        normals = append( normals, normals[-1:], 0 )
        
        
        ## Update positions, moving them in the normal direction.
        radii = emptyf( len( new_spine ) )
        for i, pt in enumerate( new_spine ):
            normal = normals[i]
            
            isects = ray_lineloop_intersection_ts( ( pt, normal ), silhouette, ray_is_line = True )
            
            pt += .5 * ( max( isects )[0] + min( isects )[0] ) * normal
            
            radii[i] = .5 * ( max( isects )[0] - min( isects )[0] )
        
        
        if polyline_self_intersects( new_spine ):
            print('Resulting spine self-intersects!  Aborting.')
            return
        
        
        ## Update circle radii.
        ## NOTE: We must do this before we smooth/resample the spine, because
        ##       the 'radii' array is tied to the current spine.
        for t in self.radii.keys():
            ## Set the radius to the average of the two closest +- intersection points.
            self.radii[ t ] = radii[ self.edge_index_at_t( t ) ]
        
        
        self.cap_radii[ 0 ] = new_start_cap_radius
        self.cap_radii[ 1 ] = new_end_cap_radius
        
        
        ## Q: Smooth spine2d?
        ## A: Yes, and resample.  Not resampling creates a lot of jitter in the outline.
        #resample_N = len( new_spine )
        resample_N = max( 3, int(ceil( 20 * line_strip_length( new_spine ) )) )
        self.set_spine2d( smooth_open_stroke( new_spine, smooth = 4, resample_N = resample_N ) )
    
    def get_all_particular_ts( self ):
        '''
        Returns a list of all particular t values at which values are affixed, such as tilts,
        radii, and cross sections.
        '''
        
        result = set()
        result.update( iter(self.tilt_radians.keys()) )
        result.update( iter(self.radii.keys()) )
        result.update( iter(self.cross_sections.keys()) )
        result.update( iter(self.symmetry_radians.keys()) )
        return list( result )
    
    def set_spine2d( self, spine2d, t_map = None ):
        '''
        Replaces the object's current 'spine2d' with the parameter.
        If 't_map' is not None, it must have as keys all elements from self.get_all_particular_ts();
        it will be used to map all t's to which values are affixed along the curve, such as
        tilts, radii, and cross sections.
        
        used
        '''
        
        N = len( spine2d )
        ## Any line must have more than 1 point.
        assert N > 1
        
        ## Re-map attributes' t's.
        if t_map is None:
            tilt_radians = self.tilt_radians
            radii = self.radii
            cross_sections = self.cross_sections
            symmetry_radians = self.symmetry_radians
        else:
            tilt_radians = {}
            for t, val in self.tilt_radians.items(): tilt_radians[ t_map[t] ] = val
            radii = {}
            for t, val in self.radii.items(): radii[ t_map[t] ] = val
            cross_sections = {}
            for t, val in self.cross_sections.items(): cross_sections[ t_map[t] ] = val
            symmetry_radians = {}
            for t, val in self.symmetry_radians.items(): symmetry_radians[ t_map[t] ] = val
        
        ## If 'N' is so small that user constraints collide, return immediately (silently fail).
        ## TODO Q: Should this be handled automatically by the optimizer (increasing resolution automatically)?
        for ts in [ list(tilt_radians.keys()), list(radii.keys()), list(cross_sections.keys()), list(symmetry_radians.keys()) ]:
            if len( ts ) != len( set( [
                t2edge_index( t, spine2d )
                for t in ts
                ] ) ):
                ## It is an error to try to shrink the number of vertices in the spine if it would lead
                ## to collisions; if there are already collisions, it's not an error but a warning.
                if N < len( self.spine2d ):
                    print("Can't SweepShape.set_spine2d() with only %d vertices; tilted circles collide." % (N,))
                    return
                else:
                    print("WARNING: Setting spine2d with only %d vertices; tilted circles collide." % (N,))
        
        ## Reset generated attributes.
        self.everything_changed()
        
        ## Replace the spine.
        self.spine2d = spine2d
        
        ## Update attributes.
        if t_map is not None:
            self.tilt_radians = tilt_radians
            self.radii = radii
            self.cross_sections = cross_sections
            self.symmetry_radians = symmetry_radians
        
        self.assert_self()
    
    def resample( self, N ):
        '''
        Resamples the spine2d so that it has 'N' vertices.
        Silently does nothing if 'N' is so small that tilted circles can't be placed
        independently of each other.
        
        NOTE: Callers can detect success by checking whether len( self.spine2d ) == N afterwards.
        
        used
        '''
        
        #print "SweepShape.resample() from %d to %d vertices." % ( len(self.spine2d), N )
        
        ## Any line must have more than 1 point.
        assert N > 1
        
        ## If we already have N vertices, we're done.
        if N == len( self.spine2d ): return
        
        ## TODO Q: Smooth spine2d?
        new_spine = smooth_open_stroke( self.spine2d, resample_N = N )
        #new_spine = resample_line_strip( self.spine2d, N )
        ## In case smoothing decided to change the number of samples.
        N = len( new_spine )
        
        ## Convert t's (parameterization of curve by the vertices)
        ## to
        ## s's (arc length parameterization of the curve).
        ## TODO Q: Why does this sometimes causes t's at the start or end (I'm guessing the end) to shift one place towards the middle?
        old_ts = self.get_all_particular_ts()
        new_ts = [
            vertex_index2s(
                max( 0, min( N-1, int(round(
                    t2vertex_index( t, self.spine2d ) * float(N-1) / float( len(self.spine2d)-1 )
                    )) ) ), new_spine ) for t in old_ts
            ]
        
        self.set_spine2d( new_spine, dict( list(zip( old_ts, new_ts )) ) )
    
    def smooth_spine2d( self, smooth_iterations = 10 ):
        '''
        Smoothes the spine2d 'smooth_iterations' times.
        
        used
        '''
        
        ## Smooth and resample.  Resampling also removes jitter in the outline.
        new_spine = smooth_open_stroke( self.spine2d, smooth = smooth_iterations, resample_N = len( self.spine2d ) )
        #new_spine = smooth_open_stroke( self.spine2d, smooth = smooth_iterations, resample_N = None )
        
        ## Convert t's (parameterization of curve by the vertices)
        ## to
        ## s's (arc length parameterization of the curve).
        old_ts = self.get_all_particular_ts()
        new_ts = [
            vertex_index2s( t2vertex_index( t, self.spine2d ), new_spine ) for t in old_ts
            ]
        
        self.set_spine2d( new_spine, dict( list(zip( old_ts, new_ts )) ) )
    
    def translate_2d( self, delta2d ): return self.translate( delta2d )
    def translate( self, delta ):
        '''
        Given a 2d or 3d offset vector 'delta', translate this shape.
        '''
        
        assert len( delta ) in (2,3)
        
        delta3d = zerosf(3)
        delta3d[:len(delta)] = delta
        
        self.spine2d = asarrayf( self.spine2d )
        self.spine2d += delta3d[:2]
        #self._spine2d_hash = None
        
        if self.mesh is not None:
            self.mesh.vs = asarrayf( self.mesh.vs )
            self.mesh.vs += delta3d
        
        if self.spine3d is not None:
            self.spine3d[:] += delta3d
        
        if self.boundary2d is not None:
            self.boundary2d[:] += delta3d[:2]
    
    def gen_spine3d( self ):
        #print "====== gen_spine3d() ======="
        
        assert len( self.spine2d ) > 0
        assert len( self.tilt_radians ) > 0
        self.assert_self()
        
        ## Create a 3d spine from the 2d spine.
        ## 1 Create a 3d array out of the 2d spine.
        ## 2 Set edge vector constraints from the tilt radians.
        ## 3 Fix all xy and set a single z = 0 constraint to fix the 3d spine in space.
        ## 4 Solve for the z values along the spine.
        
        from sweep_spine import optimize_shell_with_energy_and_constraint_functions
        import sweep_spine
        
        ## 1
        dim = 3
        spine3d = zerosf( ( len( self.spine2d ), 3 ) )
        spine3d[:,:2] = self.spine2d
        ## Use the last optimized z values as a starting guess for our new optimization.
        if self._spine3d is not None and len( self.spine2d ) == len( self._spine3d ):
            spine3d[:,2] = self._spine3d[:,2]
        
        ## 2
        ## We will use soft (least square) z difference constraints.
        eqs = []
        for t, tilt_radians in self.tilt_radians.items():
            center_edge_index = self.edge_index_at_t( t )
            ## Following the comments in __init__().
            tangent2d = spine3d[center_edge_index+1][:2] - spine3d[center_edge_index][:2]
            zdiff = tilt_radians2tangent3d( tilt_radians, tangent2d )[2]
            #print 'spine3d[%s][2]: %s' % ( center_edge_index+1, spine3d[center_edge_index+1][2] )
            eqs.append( (
                [
                    ( (center_edge_index + 1)*dim + 2, 1. ),
                    ( center_edge_index*dim + 2, -1. )
                ],
                zdiff
                ) )
        
        ## 3
        def constrain_grad_hess( grad, hess ):
            ## Apply constraints, by setting constrained degrees-of-freedom (all xy's and
            ## the first point's z)
            ## to have rows in the hessian from the identity matrix and gradients of 0.
            ## NOTE: Settings gradients to zero assumes that the current values already have the
            ##       constraints baked in, as hard constraints.
            for constraint_vertex_index in range(len(spine3d)):
                for c in range(2):
                    if hess is not None:
                        hess[ constraint_vertex_index*dim + c, : ] = 0.
                        hess[ constraint_vertex_index*dim + c, constraint_vertex_index*dim + c ] = 1.
                    grad[ constraint_vertex_index*dim + c ] = 0.
            
            ## Constrain the first point's z.
            if hess is not None:
                hess[ 2, : ] = 0.
                hess[ 2, 2 ] = 1.
            grad[ 2 ] = 0.
        
        ## 4
        self._spine3d = optimize_shell_with_energy_and_constraint_functions(
            spine3d,
            [
                (1., sweep_spine.Elaplacian_editing_undeformed( spine3d )),
                (1e7, sweep_spine.Elinear_lsq_constraints( eqs ))
            ],
            [ constrain_grad_hess ],
            linear = True
            )
    
    def gen_spine3d_hash( self ):
        return ( self.spine2d_hash(), self.tilt_radians_hash(), )
    
    def get_spine3d( self ):
        cur_hash = self.gen_spine3d_hash()
        #print "self._spine3d_hash:", self._spine3d_hash
        #print "cur_hash:", cur_hash
        #print "self._spine3d_hash != cur_hash:", self._spine3d_hash != cur_hash
        if self._spine3d_hash != cur_hash: self._spine3d = None
        if self._spine3d is None:
            self.gen_spine3d()
            self._spine3d_hash = cur_hash
        return self._spine3d
    
    spine3d = property( get_spine3d )
    
    
    def gen_spine3d_diffuser( self ):
        from sweep_spine import diffuse_values_along_sweep_spine_factored
        self._spine3d_diffuser = diffuse_values_along_sweep_spine_factored( self.spine3d )
    
    def gen_spine3d_diffuser_hash( self ):
        return ( self.spine2d_hash(), self.tilt_radians_hash(), )
    
    def get_spine3d_diffuser( self ):
        cur_hash = self.gen_spine3d_diffuser_hash()
        #print "self._spine3d_hash:", self._spine3d_hash
        #print "cur_hash:", cur_hash
        #print "self._spine3d_hash != cur_hash:", self._spine3d_hash != cur_hash
        if self._spine3d_diffuser_hash != cur_hash: self._spine3d_diffuser = None
        if self._spine3d_diffuser is None:
            self.gen_spine3d_diffuser()
            self._spine3d_diffuser_hash = cur_hash
        return self._spine3d_diffuser
    
    spine3d_diffuser = property( get_spine3d_diffuser )
    
    
    def gen_dense_radii( self ):
        #print "====== gen_dense_radii() ======="
        
        assert len( self.radii ) > 0
        self.assert_self()
        
        ## Smooth radii along the spine.
        radii_constraints = [
            ( self.edge_index_at_t( t ), (radius,) )
            for t, radius in self.radii.items()
            ]
        ## TODO Q: Follow mesh_for_shape() to add hemisphere caps?
        self._dense_radii = self.spine3d_diffuser( radii_constraints ).ravel()
        self._dense_radii = maximum( kMinRadius, self._dense_radii )
    
    def gen_dense_radii_hash( self ):
        return ( self.spine2d_hash(), self.tilt_radians_hash(), self.radii_hash() )
    
    def get_dense_radii( self ):
        cur_hash = self.gen_dense_radii_hash()
        #print "self._dense_radii_hash:", self._dense_radii_hash
        #print "cur_hash:", cur_hash
        #print "self._dense_radii_hash != cur_hash:", self._dense_radii_hash != cur_hash
        if self._dense_radii_hash != cur_hash: self._dense_radii = None
        if self._dense_radii is None:
            self.gen_dense_radii()
            self._dense_radii_hash = cur_hash
        return self._dense_radii
    
    dense_radii = property( get_dense_radii )
    
    def gen_dense_cross_sections( self ):
        #print "====== gen_dense_cross_sections() ======="
        
        self.assert_self()
        
        if len( self.cross_sections ) == 0:
            self._dense_cross_sections = zeros( ( 0,kCrossSectionSamples,2 ) )
            return
        
        sparse_xsections = [
                ( self.edge_index_at_t( t ), asarray( xsection ).ravel() )
                for t, xsection in self.cross_sections.items()
            ]
        
        ## Duplicate cross sections so we have gradient = 0 constraints.
        sparse_xsections = dict( sparse_xsections )
        extra_xsections = []
        for edge_index, xsection in sparse_xsections.items():
            if edge_index+1 not in sparse_xsections:
                extra_xsections.append( ( edge_index+1, xsection ) )
            elif edge_index-1 not in sparse_xsections:
                extra_xsections.append( ( edge_index-1, xsection ) )
            else:
                print("WARNING: Unable to set cross section gradient constraint.")
        
        sparse_xsections = list(sparse_xsections.items()) + extra_xsections
        
        num_cross_section_samples = len( list(self.cross_sections.values())[0] )
        assert num_cross_section_samples == kCrossSectionSamples
        dense_xsections = self.spine3d_diffuser( sparse_xsections )
        dense_xsections = dense_xsections.reshape( ( dense_xsections.shape[0], num_cross_section_samples, 2 ) )
        self._dense_cross_sections = dense_xsections
    
    def gen_dense_cross_sections_hash( self ):
        return ( self.spine2d_hash(), self.cross_sections_hash() )
    
    def get_dense_cross_sections( self ):
        cur_hash = self.gen_dense_cross_sections_hash()
        if self._dense_cross_sections_hash != cur_hash: self._dense_cross_sections = None
        if self._dense_cross_sections is None:
            self.gen_dense_cross_sections()
            self._dense_cross_sections_hash = cur_hash
        return self._dense_cross_sections
    
    dense_cross_sections = property( get_dense_cross_sections ) 
    
    def gen_mesh( self ):
        #print "====== gen_mesh() ======="
        
        self.assert_self()
        
        if 0 == len( self.cross_sections ):
            if kNoPinch:
                from sweep_surface import mesh_from_sweep_circles_no_pinch
                self._mesh = mesh_from_sweep_circles_no_pinch( self.spine3d, self.dense_radii, caps = self.get_caps() )
            else:
                from sweep_surface import mesh_from_sweep_circles
                self._mesh = mesh_from_sweep_circles( self.spine3d, self.dense_radii, caps = self.get_caps() )
        else:
            from sweep_surface import mesh_from_sweep_cross_sections
            
            dense_xsections = self.dense_cross_sections * self.dense_radii[:,newaxis,newaxis]
            
            assert len( self.symmetry_directions ) == len( dense_xsections )
            self._mesh = mesh_from_sweep_cross_sections(
                self.spine3d,
                list(zip( self.symmetry_directions, dense_xsections )),
                caps = self.get_caps(),
                no_pinch = kNoPinch
                )
    
    def gen_mesh_hash( self ):
        ## The mesh only depends on self.symmetry_radians if self.cross_sections is non-empty.
        return ( self.spine2d_hash(), self.tilt_radians_hash(), self.radii_hash(), self.cross_sections_hash(), self.symmetry_radians_hash(), self.cap_radii_hash() )
    
    def get_mesh( self ):
        cur_hash = self.gen_mesh_hash()
        if self._mesh_hash != cur_hash: self._mesh = None
        if self._mesh is None:
            self.gen_mesh()
            self._mesh_hash = cur_hash
        return self._mesh
    
    mesh = property( get_mesh )
    
    
    def gen_boundary2d( self ):
        #print "====== gen_boundary2d() ======="
        
        self._boundary2d = asarray( gen_offset_boundary2d( self.spine2d, self.dense_radii, self.get_caps() ) )
    
    def gen_boundary2d_hash( self ):
        return ( self.spine2d_hash(), self.tilt_radians_hash(), self.radii_hash(), self.cap_radii_hash() )
    
    def get_boundary2d( self ):
        cur_hash = self.gen_boundary2d_hash()
        if self._boundary2d_hash != cur_hash: self._boundary2d = None
        if self._boundary2d is None:
            self.gen_boundary2d()
            self._boundary2d_hash = cur_hash
        return self._boundary2d
    
    boundary2d = property( get_boundary2d )
    
    
    def gen_symmetry_directions( self ):
        #print "====== gen_symmetry_directions() ======="
        
        #sym_dir = 'old'
        #sym_dir = 'minimal_rotation_frame'
        sym_dir = 'hybrid'
        if sym_dir == 'old':
            ## 1 Find the distance and closest point to 'pt' on all lines passing through points
            ##   along self.spine3d in the 3d direction of its symmetry plane.
            ##   This 3d direction lies in the plane spanned by the 2d normal (with z = 0)
            ##   and the cross product of that vector with the 3d tangent.
            
            ## 1
            ## Pass spine3d, so we get proper tangents and 2d normals with bogus z values.
            tangents, normals = tangents_and_normals_along_line_strip_2d( self.spine3d )
            ## Set the normals' z values to 0, so we can use them as 2d normals.
            normals[:,2] = 0.
            ## Repeat the last normal and tangent, because there is one less edge than point.
            tangents = append( tangents, tangents[-1:], 0 )
            normals = append( normals, normals[-1:], 0 )
            
            ## Get the vector orthogonal to normals and tangents; the 3d symmetry direction lies in
            ## this plane.
            ups = cross( normals, tangents )
            
            cos_theta = cos( list(self.symmetry_radians.values())[0] )
            sin_theta = sin( list(self.symmetry_radians.values())[0] )
            assert len( ups ) == len( normals )
            self._symmetry_directions = [
                cos_theta * up + sin_theta * normal
                for up, normal in zip( ups, normals )
                ]
        
        elif sym_dir == 'minimal_rotation_frame':
            from sweep_surface import smooth_orthogonal_vectors_from_constraints
            
            sparse_orthogonal_vectors = []
            for t, theta in self.symmetry_radians.items():
                unit_tangent3d = dir( self.tangent3d_at_t( t ) )
                normal = cross( unit_tangent3d, (0,0,1) )
                assert abs( normal[2] ) < kEps
                up = cross( normal, unit_tangent3d )
                
                symmetry_dir = cos( theta ) * up + sin( theta ) * normal
                sparse_orthogonal_vectors.append( ( self.edge_index_at_t( t ), symmetry_dir ) )
            
            self._symmetry_directions = smooth_orthogonal_vectors_from_constraints(
                self.spine3d,
                sparse_orthogonal_vectors,
                path_diffuser = self.spine3d_diffuser
                )
        
        elif sym_dir == 'hybrid':
            from sweep_surface import smooth_orthogonal_vectors_from_constraints
            
            ups, normals, tangents = self.get_frame_components()
            frames = list(zip( ups, normals, tangents ))
            
            sparse_orthogonal_vectors = []
            for t, theta in self.symmetry_radians.items():
                edge_index = self.edge_index_at_t( t )
                
                ## edge_index is OK because ups and normals simply repeat the last value.
                symmetry_dir = cos( theta ) * ups[ edge_index ] + sin( theta ) * normals[ edge_index ]
                sparse_orthogonal_vectors.append( ( edge_index, symmetry_dir ) )
            
            self._symmetry_directions = smooth_orthogonal_vectors_from_constraints(
                self.spine3d,
                sparse_orthogonal_vectors,
                sweep_path_frames = frames,
                path_diffuser = self.spine3d_diffuser
                )
        
        else:
            raise RuntimeError("Unknown method of calculating symmetry directions: " + sym_dir)
    
    def gen_symmetry_directions_hash( self ):
        return ( self.spine2d_hash(), self.tilt_radians_hash(), self.symmetry_radians_hash() )
    
    def get_symmetry_directions( self ):
        cur_hash = self.gen_symmetry_directions_hash()
        if self._symmetry_directions_hash != cur_hash: self._symmetry_directions = None
        if self._symmetry_directions is None:
            self.gen_symmetry_directions()
            self._symmetry_directions_hash = cur_hash
        return self._symmetry_directions
    
    symmetry_directions = property( get_symmetry_directions )
    
    def get_frame_components( self ):
        ## Pass spine3d, so we get proper 3d tangents and 2d normals with bogus z values.
        tangents, normals = tangents_and_normals_along_line_strip_2d( self.spine3d )
        ## Set the normals' z values to 0, so we can use them as 2d normals.
        normals[:,2] = 0.
        ## Normalize the planar normals
        normals = normals / sqrt( ( normals**2 ).sum(1) )[:,newaxis]
        ## Repeat the last normal and tangent, because there is one less edge than point.
        tangents = append( tangents, tangents[-1:], 0 )
        normals = append( normals, normals[-1:], 0 )
        
        ## Get the vector orthogonal to normals and tangents; the 3d symmetry direction lies in
        ## this plane.
        ups = cross( normals, tangents )
        
        return ( ups, normals, tangents )
    
    def update_3d( self ):
        '''
        used
        '''
        
        #print "SweepShape.update_3d()"
        
        assert len( self.spine2d ) > 0
        assert len( self.tilt_radians ) > 0
        assert len( self.radii ) > 0
        assert len( self.symmetry_radians ) > 0
        self.assert_self()
        
        ## This call is now a no-op, and properties are regenerated as-needed.
        
        #self.gen_spine3d()
        #self.gen_dense_radii()
        #self.gen_mesh()
        #self.gen_boundary2d()
        #self.gen_symmetry_directions()
        
        #self.everything_changed()
        
        '''
        #if not hasattr( self, '_decs' ): self._decs = []
        self._decs = []
        self._decs.extend( [
            ( 'line', (232./255, 146./255, 26./255), ray[0], ray[0] + ray[1] )
            for ray in self.symmetry_plane_rays_at_vertices()
            ] )
        '''
    
    def mirror_plane( self, pt ):
        '''
        Given a 3d point 'pt', returns the tangent plane at the closest point
        on the shape's symmetry sheet.
        
        used
        '''
        
        pt = asarrayf( pt ).ravel()
        assert len( pt ) == 3
        
        ## 1 Find the closest cross section line loop to pt.
        ## 2 Return the local symmetry plane at the spine point for that cross section.
        
        ## 1
        ## Get 3D tangents
        tangents = zeros( self.spine3d.shape )
        tangents[:-1] = self.spine3d[1:] - self.spine3d[:-1]
        tangents[-1] = tangents[-2]
        
        distsqrs = []
        if len( self.dense_cross_sections ) > 0:
            ## untested
            assert len( self.dense_cross_sections ) == len( self.spine3d )
            for spine_pt, tangent, cross_section in zip( self.spine3d, tangents, self.dense_cross_sections ):
                ## Find the rotation matrix that takes (0,0,1) to tangent.
                R = get_rotation_matrix( *get_axis_angle_rotation( tangent, (0,0,1) ) )
                ## Apply the rotation matrix to the cross section.
                cross_section3d = zeros( ( len( cross_section ), 3 ) )
                cross_section3d[:,:2] = cross_section
                cross_section3d = ( R @ cross_section3d.T ).T
                cross_section3d += spine_pt
                
                distsqr = closest_distsqr_and_point_and_edge_index_on_line_loop_to_point( cross_section3d, pt )[0]
                distsqrs.append( distsqr )
        else:
            ## used
            for spine_pt, tangent, radius in zip( self.spine3d, tangents, self.dense_radii ):
                ## Find the rotation matrix that takes the tangent to (0,0,1).
                R = asarray( get_rotation_matrix( *get_axis_angle_rotation( (0,0,1), tangent ) ) )
                ## Apply the rotation matrix to the spine point and pt
                Rspine = dot( R, spine_pt )
                Rpt = dot( R, pt )
                ## The distance^2 is the 2D point-to-circle distance^2 (ignore the z component),
                ## plus the difference^2 of z components.
                ## The distance^2 from a point to a circle is the absolute value of the
                ## point-to-point distance^2 minus the radius^2.
                d = Rspine[:2] - Rpt[:2]
                d = abs( dot( d, d ) - radius**2 )
                d += ( Rspine[2] - Rpt[2] )**2
                distsqrs.append( d )
        
        closest_pt_i = argmin( distsqrs )
        
        
        ## 2
        closest_plane = self.symmetry_plane_tangent_plane_at_vertex_index( closest_pt_i )
        return closest_plane
        
        
        ########################
        
        '''
        ## 1 Find the distance and closest point to 'pt' on all symmetry plane rays.
        ## 2 Reflect 'pt' across the plane created by its closest point and that point's
        ##   neighbor.
        
        ## 1
        symmetry_plane_rays = self.symmetry_plane_rays_at_vertices()
        lines = [ ( ray[0], ray[0] + ray[1] ) for ray in symmetry_plane_rays ]
        assert len( lines ) == len( self.spine3d )
        
        closest_pts = [
            closest_point_on_line_to_point( line, pt, line_is_segment = False )
            for line in lines
            ]
        
        ## 2
        closest_pt_i = argmin( [ mag2( pt - close_pt ) for close_pt in closest_pts ] )
        closest_plane = self.symmetry_plane_tangent_plane_at_vertex_index( closest_pt_i )
        return closest_plane
        '''
    
    def mirror_point( self, pt ):
        '''
        Given a 3d point 'pt', returns the reflection of 'pt' across this shape's symmetry sheet.
        
        used
        '''
        closest_plane = self.mirror_plane( pt )
        return reflect_point_across_plane( pt, closest_plane )
    
    def symmetry_plane_rays_at_vertices( self ):
        '''
        Returns 3d rays on the symmetry plane sheet, emanating from the 3d spine,
        one per point on self.spine3d.
        NOTE: The rays are not normalized.
        
        used
        '''
        
        return list(zip( self.spine3d, self.symmetry_directions_at_vertices() ))
    
    def symmetry_plane_ray_at_t( self, t ):
        '''
        Given t along the spine curve,
        returns a ray originating from the spine3d point in the symmetry plane.
        NOTE: The ray is not normalized.
        
        untested
        '''
        
        vi = self.vertex_index_at_t( t )
        return self.symmetry_plane_rays_at_vertices()[ vi ]
    
    def symmetry_plane_tangent_plane_at_vertex_index( self, index ):
        '''
        Given the index of a point along the spine curve,
        returns the tangent plane of the symmetry manifold passing through the spine3d point.
        
        used
        '''
        
        symmetry_plane_rays = self.symmetry_plane_rays_at_vertices()
        ## Clip 'index' so that it can be used as an edge index.
        edge_index = min( len( symmetry_plane_rays )-2, index )
        ## Find the 3d tangent vector.
        tangent = symmetry_plane_rays[ edge_index + 1 ][0] - symmetry_plane_rays[ edge_index ][0]
        ## The plane is the plane spanned by the symmetry plane ray and the tangent vector.
        sym_plane = (
            symmetry_plane_rays[ index ][0],
            cross( tangent, symmetry_plane_rays[ index ][1] )
            )
        return sym_plane
    
    def candidate_origins( self ):
        '''
        Returns a sequence of possible origin points, namely the 3d spine endpoints.
        '''
        
        return [ self.spine3d[0], self.spine3d[-1], self.spine3d[ len( self.spine3d ) // 2 ] ]
    
    def center( self ):
        '''
        Returns the center of the shape.
        '''
        
        #return average( self.spine3d, axis = 0 )
        
        ## This is always inside the shape, which is usually concave.
        return self.spine3d[ len(self.spine3d) // 2 ]
    
    def length3d( self ):
        '''
        Returns the 3d length of the spine.
        '''
        return line_strip_length( self.spine3d )
    
    def set_length3d( self, length, origin = (0,0,0) ):
        '''
        Sets the 3d length of the spine.
        '''
        assert length > 0.
        
        assert len( origin ) == 3
        origin = asarrayf( origin ).ravel()
        
        scale = length / self.length3d()
        assert scale > 0.
        
        ## We can't directly scale the 3d spine, but scaling just the 2d size by the desired factor
        ## should also scale the 3d spine by the same factor, since the z coordinates should be
        ## scale independent (or linearly scale dependent).
        #old_spine2d0 = array( self.spine2d[0] )
        #old_spine3d0 = array( self.spine3d[0] )
        spine3d0 = scale * ( self.spine3d[0] - origin ) + origin
        spine2d = scale * ( asarrayf( self.spine2d ) - origin[:2] ) + origin[:2]
        self.set_spine2d( spine2d )
        #self.shape.resample( max( 3, int(round( scale * len( self.shape.spine2d ) )) ) )
        
        ## For some reason accessing self.spine3d fixed the depth placement bug.
        ## TODO Q: Why is this?  Is there some private value being accessed that's
        ##         not being invalidated when spine2d changes?
        spine3d0cur = self.spine3d[0]
        ## This can fail if spine3d was not quite spine2d. (Which seems to happen when turning on turbo mode.)
        #expect( mag2( spine3d0cur[:2] - spine3d0[:2] ) < 1e-5 )
        expect( mag2( spine3d0cur[:2] - spine2d[0] ) < 1e-5 )
        #print "spine3d0cur:", spine3d0cur
        #print "spine3d0:", spine3d0
        ## This allows us to connect first and then apply same size constraints.
        self.offset_in_depth( spine3d0cur[2] - spine3d0[2] )
    
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
    
    def offset_in_depth( self, z_offset, oblique = (0.,0.) ):
        '''
        Translate the object along the look vector (oblique[0], oblique[1], -1) by 'z_offset'.
        Inverts the oblique projection, so the object maintains the same (projected) appearance.
        '''
        
        self.mesh.vs = asarray( self.mesh.vs )
        self.mesh.vs[:,:] += z_offset * array( ( oblique[0], oblique[1], -1 ) )
        self.mesh.positions_changed()
        
        ## We need to keep spine3d up-to-date because it is used for things such as the
        ## symmetry plane.
        self.spine3d[:] += z_offset * array( ( oblique[0], oblique[1], -1 ) )
