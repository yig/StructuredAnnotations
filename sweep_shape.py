from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *
from vecutil import *

from pprint import pprint

try: from pydb import debugger
except ImportError:
    import pdb
    def debugger(): pdb.set_trace()

def mesh_for_shape( silhouette, tilted_circles, symmetry_xsection, symmetry_xsection_front_to_back ):
    '''
    Given a simple, closed sequence of 2d points for the silhouette 'silhouette',
    a sequence of tilted circles 'tilted_circles' in the form of 3-tuples ( 2/3d circle center, radius, 3d circle normal ),
    a sequence of 2d points for the symmetry cross-section curve 'symmetry_xsection' and optional index 'symmetry_xsection_front_to_back',
        where the symmetry curve is assumed to start at one of its poles, travel across the
        front of the shape, then, at index 'symmetry_xsection_front_to_back', switch to
        travelling across the back of the shape,
    returns a TriMesh object for the surface of the shape.
    
    If 'symmetry_xsection_front_to_back' is None, it will be automatically determined from
    'symmetry_xsection' and 'silhouette'.
    
    lightly tested
    '''
    ## This code doesn't actually look at the silhouette.
    #assert len( silhouette ) > 2
    assert len( tilted_circles ) > 0
    assert len( symmetry_xsection ) > 2
    assert symmetry_xsection_front_to_back is None or symmetry_xsection_front_to_back < len( symmetry_xsection )
    
    ## Automatically determine symmetry_xsection_front_to_back.
    symmetry_xsection_front_to_back = None
    if symmetry_xsection_front_to_back is None:
        assert len( silhouette ) > 2
        
        ## To automatically determine 'symmetry_xsection_front_to_back', use the
        ## local minima of distance between symmetry_xsection and the silhoutte that is farthest
        ## arc-length-wise from symmetry_xsection[0].
        
        ## NOTE: I am going to assume that the symmetry curve has been uniformly resampled,
        ##       so that I can just indices as distance.
        ## 1 Calculate the distance from the symmetry curve to the silhouette curve everywhere.
        ## 2 Find all local minima of distances everywhere.
        ## 3 Choose the local minima whose index's modular distance from 0 is greatest.
        
        ## 1
        from edge_distances import min_distance_to_edges
        edges = []
        for i in range( len( silhouette )-1 ):
            edges.append( ( silhouette[i], silhouette[i+1] ) )
        edges.append( ( silhouette[i], silhouette[0] ) )
        distances = min_distance_to_edges( asarray( symmetry_xsection ).T, edges )
        assert len( distances ) == len( symmetry_xsection )
        
        ## 2
        local_mins = []
        
        from helpers import topo_less
        N = len( distances )
        for i in range( N ):
            if topo_less( distances, i, (i+1) % N ) and topo_less( distances, i, (N+i-1) % N ):
                local_mins.append( i )
        assert len( local_mins ) > 0
        expect( local_mins > 1 )
        
        ## 3
        mod_dist_local_mins = []
        for i in local_mins:
            mod_dist_local_mins.append( ( min( i, N-i ), i ) )
        del N
        
        symmetry_xsection_front_to_back = max( mod_dist_local_mins )[1]
    
    ## 1 Determine the symmetry plane based on the tilted circles.
    ## 2 Connect a chain of vertices from the two poles of the symmetry curve through the tilted
    ##   circle constraints.
    ## 3 Solve for a smooth curve passing through the constraints with sweep_spine.optimize_sweep_spine2d().
    ## 4 Solve for circle radii with 0's at the poles and tilted circle radii at the tilted circle
    ##   constraints.  Or truncate the spine by a distance equal to the first and last tilted circle
    ##   radii at the start/end of the spine.
    ## 5 Build a mesh with sweep_surface.mesh_from_sweep_circles()
    
    ## 1
    plane = symmetry_plane_from_tilted_circles( tilted_circles, symmetry_xsection, symmetry_xsection_front_to_back )
    if plane is None:
        raise RuntimeError("Symmetry curve does not intersect tilted circles; can't generate mesh.")
    plane_pt, plane_normal = plane
    rot_mat = asarray( get_rotation_matrix( *get_axis_angle_rotation( plane_normal, (0,0,1) ) ) )
    
    ## 2
    ## Project the symmetry curve's poles and the circle centers onto the symmetry plane.
    ## TODO: Shouldn't this take into account the symmetry curve itself?
    oblique_dir = (0,0,-1)
    start = ray_plane_intersection( ( append( symmetry_xsection[0], 0. ), oblique_dir ), ( plane_pt, plane_normal ), ray_is_line = True )
    end = ray_plane_intersection( ( append( symmetry_xsection[symmetry_xsection_front_to_back], 0. ), oblique_dir ), ( plane_pt, plane_normal ), ray_is_line = True )
    
    ## TODO Q: This assumes the circle centers are ordered from start to end.
    ## A: Sort circle centers according to their projection onto the straight line from start to end.
    circle_centers = [
        ray_plane_intersection( ( append( circle_center[:2], 0. ), oblique_dir ), ( plane_pt, plane_normal ), ray_is_line = True )
        for circle_center, circle_radius, circle_normal in tilted_circles
        ]
    
    start2d = dot( rot_mat, start )[:2]
    end2d = dot( rot_mat, end )[:2]
    circle_centers2d = [ dot( rot_mat, center )[:2] for center in circle_centers ]
    
    ## Five points for each segment.
    N = 10
    constraints = {}
    line = []
    
    constraints[0] = tuple( start2d )
    total_length = 0.
    segment_length = mag( asarray( start2d ) - asarray( circle_centers2d[0] ) ) if len( circle_centers2d ) > 0 else 1.
    line.extend( [ (x, 0.) for x in linspace( 0, segment_length, N ) ] )
    total_length += segment_length
    
    prev_pos2d = start2d
    circle_constraint_indices = []
    for i, pos2d in enumerate( circle_centers2d ):
        next_pos2d = asarray( circle_centers2d[i+1] if i+1 < len( circle_centers2d ) else end2d )
        
        line_i = len(line)-1
        circle_constraint_indices.append( line_i )
        
        constraints[ line_i ] = tuple( pos2d )
        segment_length = mag( asarray( pos2d ) - next_pos2d )
        line.extend( [ (x, 0.) for x in linspace( total_length, total_length + segment_length, N )[1:] ] )
        total_length += segment_length
        
        ## Also constrain the circle normal.
        ## TODO Q: What is the right amount of gradient magnitude to use here?
        line_i_minus_1_mag = mag( asarray( line[ line_i-1 ] ) - asarray( line[ line_i ] ) )
        line_i_plus_1_mag = mag( asarray( line[ line_i+1 ] ) - asarray( line[ line_i ] ) )
        
        circle_normal2d = dir( dot( rot_mat, tilted_circles[i][2] )[:2] )
        
        if dot( circle_normal2d, prev_pos2d - pos2d ) < 0.: circle_normal2d = -circle_normal2d
        constraints[ line_i - 1 ] = tuple( pos2d + line_i_minus_1_mag * circle_normal2d )
        
        if dot( circle_normal2d, next_pos2d - pos2d ) < 0.: circle_normal2d = -circle_normal2d
        #constraints[ line_i + 1 ] = tuple( pos2d + line_i_plus_1_mag * circle_normal2d )
        
        prev_pos2d = pos2d
    
    constraints[ len(line)-1 ] = tuple( end2d )
    
    #pprint( "constraints:" )
    #pprint( constraints )
    
    ## 3
    import sweep_spine
    spine2d = sweep_spine.optimize_sweep_spine2d( line, list(constraints.items()) )
    spine3d = [ dot( rot_mat.T, append( pos2d, 0. ) ) for pos2d in spine2d ]
    
    ## 4
    constraints_radii = {}
    ## Constrain y coordinates circle centers.
    for i, index in enumerate( circle_constraint_indices ):
        constraints_radii[index] = [tilted_circles[i][1]]
    
    caps = 'hemisphere'
    if caps == 'zero':
        ## Constrain the start and end to 0.
        constraints_radii[0] = [0.]
        constraints_radii[ len(line)-1 ] = [0.]
    
    elif caps == 'hemisphere':
        ## Constrain the start and end of the line to be hemispheres
        ## 1 For the start of the line, constrain the first few radii to increase as a square root
        ##   function, like perpendicular slices of a sphere.  Do this up to the radius of the
        ##   first circular cross section.  If there isn't enough arc length before the first
        ##   circular cross section, then scale the square root function by the constant function
        ##   so that the radii reach the first cross section radius.
        ## 2 Do similar operations for the end of the line.
        
        assert len( tilted_circles ) > 0
        assert len( tilted_circles ) == len( circle_constraint_indices )
        
        ## 1
        assert circle_constraint_indices[0] > 0
        
        start_radius = tilted_circles[0][1]
        assert start_radius > 0.
        
        start_section = asarrayf( spine3d[ : circle_constraint_indices[0]+1 ] )
        ## Parameterize by arc length.
        ts = cumsum( [0.] + list( sqrt(((start_section[1:] - start_section[:-1] )**2).sum(1)) ) ) / start_radius
        ## If we didn't reach 1, force it by scaling.
        if ts[-1] < 1.:
            ts /= ts[-1]
        ## Otherwise, clip it at 1.
        else:
            ts = ts[ : list( ts > 1 ).index(True) ]
        
        ## We want to go from 1. to 0., not 0. to 1.
        ts = 1. - ts
        vals = start_radius * sqrt( 1. - ts**2 )
        
        for i in range(len( ts )):
            constraints_radii[i] = [vals[i]]
        
        assert abs( constraints_radii[0][0] ) < kEps
        
        
        ## 2
        assert circle_constraint_indices[-1] > 0
        
        start_radius = tilted_circles[-1][1]
        assert start_radius > 0.
        
        start_section = asarrayf( spine3d[ circle_constraint_indices[-1]: ] )
        ## Parameterize by arc length.
        ## This line differs from the start-cap case, because we accumulate the length along the reverse direction.
        ts = cumsum( ( list( sqrt(((start_section[1:] - start_section[:-1] )**2).sum(1)) ) + [0.] )[::-1] ) / start_radius
        ## If we didn't reach 1, force it by scaling.
        if ts[-1] < 1.:
            ts /= ts[-1]
        ## Otherwise, clip it at 1.
        else:
            ts = ts[ : list( ts > 1 ).index(True) ]
        
        ## We want to go from 1. to 0., not 0. to 1.
        ts = 1. - ts
        vals = start_radius * sqrt( 1. - ts**2 )
        
        for i in range(len( ts )):
            ## This line differs from the start-cap case, because we undo the reversal of length accumulation.
            constraints_radii[ len(line)-1 - i ] = [vals[i]]
        
        assert abs( constraints_radii[ len(line)-1 ][0] ) < kEps
    
    else:
        raise RuntimeError("Unknown cap type: %s" % (caps,))
    
    radii = sweep_spine.diffuse_values_along_sweep_spine( line, list(constraints_radii.items()) ).ravel()
    
    assert abs( radii[0] ) < kEps and abs( radii[-1] ) < kEps
    ## Q: The first and last radii are 0.  Should I drop them?
    ## A1: Yes, sweep_surface.mesh_from_sweep_circles() can't handle start or end caps with 0 radius.
    ##     UPDATE: No need to drop them; sweep_surface.mesh_from_sweep_circles() can now handle
    ##             start or end caps with 0 radius.
    #spine3d = spine3d[1:-1]
    #radii = radii[1:-1]
    #expect( abs( radii[0] ) > kEps and abs( radii[-1] ) > kEps )
    
    ## TODO Q: What if there are radii less than 0?  I must drop those!
    assert False not in [ radius > -kEps for radius in radii ] or debugger()
    
    import sweep_surface
    return sweep_surface.mesh_from_sweep_circles( spine3d, radii, caps = 'flat' )


def symmetry_plane_from_tilted_circles( tilted_circles, symmetry_xsection, symmetry_xsection_front_to_back ):
    '''
    Given a sequence of tilted circles 'tilted_circles' in the form of 3-tuples ( 2/3d circle center, radius, 3d circle normal ),
    a sequence of 2d points for the symmetry cross-section curve 'symmetry_xsection' and index 'symmetry_xsection_front_to_back',
        where the symmetry curve is assumed to start at one of its poles, travel across the
        front of the shape, then, at index 'symmetry_xsection_front_to_back', switch to
        travelling across the back of the shape,
    returns a plane (3d point, 3d normal).
    Returns None if the symmetry curve and tilted circles do not intersect.
    
    lightly tested
    '''
    
    assert len( tilted_circles ) > 0
    
    ## Average all plane points and normals.
    plane_pts = []
    plane_normals = []
    for tilted_circle in tilted_circles:
        plane = symmetry_plane_from_tilted_circle( tilted_circle, symmetry_xsection, symmetry_xsection_front_to_back )
        if plane is None:
            print("Symmetry curve did not intersect tilted circle, skipping.")
            continue
        pt, normal = plane
        plane_pts.append( pt )
        plane_normals.append( normal )
    
    assert len( plane_pts ) == len( plane_normals )
    if len( plane_pts ) == 0: return None
    
    plane_pt = asarrayf( plane_pts ).sum(0) / len( plane_pts )
    plane_normal = dir( asarray( plane_normals ).sum(0) )
    
    return plane_pt, plane_normal

def symmetry_plane_from_tilted_circle( tilted_circle, symmetry_xsection, symmetry_xsection_front_to_back ):
    '''
    Given a tilted circle 'tilted_circle' in the form of a 3-tuple ( 2/3d circle center, radius, 3d circle normal ),
    a sequence of 2d points for the symmetry cross-section curve 'symmetry_xsection' and index 'symmetry_xsection_front_to_back',
        where the symmetry curve is assumed to start at one of its poles, travel across the
        front of the shape, then, at index 'symmetry_xsection_front_to_back', switch to
        travelling across the back of the shape,
    returns a plane (3d point, 3d normal).
    Returns None if the symmetry curve and tilted circle do not intersect.
    
    lightly tested
    '''
    
    circle_center, circle_radius, circle_normal = tilted_circle
    ## Actually, use the unrotated circle's normal
    circle_normal = arrayf( circle_normal )
    circle_normal[2] = 0.
    circle_normal = dir( circle_normal )
    
    ## 1
    ## TODO: Take oblique into account (rotate so that oblique_dir is (0,0,1) and then drop z).
    intersections = []
    front_z = circle_center[2] + circle_radius + 1
    back_z = circle_center[2] - circle_radius - 1
    for i in range(len( symmetry_xsection )):
        on_front = ( i < symmetry_xsection_front_to_back )
        z = front_z if on_front else back_z
        ray3d = (
            append( symmetry_xsection[i], z ),
            append( symmetry_xsection[ (i+1) % len( symmetry_xsection ) ], z ) - append( symmetry_xsection[i], z )
            )
        ts = [ ray_plane_intersection_t( ray3d, ( circle_center, circle_normal ) ) ]
        
        intersections.extend( [ ray3d[0] + t * ray3d[1] for t in ts if t >= 0 and t <= 1 ] )
    
    if len( intersections ) == 0: return None
    
    ## TODO: For numerical reasons this might not be true; if not, remove the points epsilon
    ##       away from each other so that it is true.
    ## TODO: This might also not be true if we don't have any back intersections.
    expect( len( intersections ) == 2 )
    
    ## 2
    for i in range(len( intersections )):
        isect_pt = intersections[i]
        intersections[i] = closest_point_on_sphere_to_ray( ( circle_center, circle_radius ), ( isect_pt, (0,0,1) ), ray_is_line = True )[0]
    
    ## 3
    front_pt = max( intersections, key = lambda el: el[2] )
    #back_pt = min( intersections, key = lambda el: el[2] )
    back_pt = circle_center
    
    ## 4
    plane_pt = front_pt
    plane_normal = dir( cross( front_pt - back_pt, circle_normal ) )
    
    return plane_pt, plane_normal

def test1():
    
    silhouette = [ ( .5 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] ]
    tilted_circles = [ ( (0.,0.,0.), .5, (0.,1.,0.) ) ]
    symmetry_xsection = [ ( .2 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] + pi/2 ]
    symmetry_xsection_front_to_back = len( symmetry_xsection ) // 2
    mesh = mesh_for_shape( silhouette, tilted_circles, symmetry_xsection, symmetry_xsection_front_to_back )
    mesh.write_OBJ( "sweep_shape-test1.obj" )

def test1rot():
    
    silhouette = [ ( .5 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] ]
    tilted_circles = [ ( (0.,0.,0.), .5, (0.,1.,0.) ) ]
    symmetry_xsection = [ ( .2 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] + pi/2 ]
    symmetry_xsection_front_to_back = len( symmetry_xsection ) // 2
    
    ## A CCW rotation in the xy plane by pi/2 radians.
    rot_mat = asarray( get_rotation_matrix( (0,0,1), pi/4 ) )
    silhouette = [ dot( rot_mat[:2,:2], pt ) for pt in silhouette ]
    symmetry_xsection = [ dot( rot_mat[:2,:2], pt ) for pt in symmetry_xsection ]
    tilted_circles = [
        ( dot( rot_mat, circle_center ), circle_radius, dot( rot_mat, circle_normal ) )
        for circle_center, circle_radius, circle_normal in tilted_circles
        ]
    
    mesh = mesh_for_shape( silhouette, tilted_circles, symmetry_xsection, symmetry_xsection_front_to_back )
    mesh.write_OBJ( "sweep_shape-test1rot.obj" )

def test2():
    
    silhouette = [ ( .5 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] ]
    tilted_circles = [
        ( (0.,.5,0.), .5, (0.,1.,0.) ),
        ( (0.,0.,0.), .5, (0.,1.,0.) ),
        ( (0.,-.5,0.), .5, (0.,1.,0.) ),
        ]
    symmetry_xsection = [ ( .2 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] + pi/2 ]
    symmetry_xsection_front_to_back = len( symmetry_xsection ) // 2
    mesh = mesh_for_shape( silhouette, tilted_circles, symmetry_xsection, symmetry_xsection_front_to_back )
    mesh.write_OBJ( "sweep_shape-test2.obj" )

def test2rot():
    
    silhouette = [ ( .5 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] ]
    tilted_circles = [
        ( (0.,.5,0.), .5, (0.,1.,0.) ),
        ( (0.,0.,0.), .5, (0.,1.,0.) ),
        ( (0.,-.5,0.), .5, (0.,1.,0.) ),
        ]
    symmetry_xsection = [ ( .2 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] + pi/2 ]
    symmetry_xsection_front_to_back = len( symmetry_xsection ) // 2
    
    ## A CCW rotation in the xy plane by pi/2 radians.
    rot_mat = asarray( get_rotation_matrix( (0,0,1), pi/4 ) )
    silhouette = [ dot( rot_mat[:2,:2], pt ) for pt in silhouette ]
    symmetry_xsection = [ dot( rot_mat[:2,:2], pt ) for pt in symmetry_xsection ]
    tilted_circles = [
        ( dot( rot_mat, circle_center ), circle_radius, dot( rot_mat, circle_normal ) )
        for circle_center, circle_radius, circle_normal in tilted_circles
        ]
    
    mesh = mesh_for_shape( silhouette, tilted_circles, symmetry_xsection, symmetry_xsection_front_to_back )
    mesh.write_OBJ( "sweep_shape-test2rot.obj" )

def test3():
    
    silhouette = [ ( .5 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] ]
    tilted_circles = [ ( (0.,0.,0.), .5, (1.,1.,0.) ) ]
    symmetry_xsection = [ ( .2 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] + pi/2 ]
    symmetry_xsection_front_to_back = len( symmetry_xsection ) // 2
    mesh = mesh_for_shape( silhouette, tilted_circles, symmetry_xsection, symmetry_xsection_front_to_back )
    mesh.write_OBJ( "sweep_shape-test3.obj" )

def test4():
    
    silhouette = [ ( .5 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] ]
    tilted_circles = [
        ( (0.,.5,0.), .2, (0.,1.,0.) ),
        ( (0.,0.,0.), .5, (0.,1.,0.) ),
        ( (0.,-.5,0.), .1, (0.,1.,0.) ),
        ]
    symmetry_xsection = [ ( .2 * cos( theta ) + kEps, .7 * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] + pi/2 ]
    symmetry_xsection_front_to_back = len( symmetry_xsection ) // 2
    mesh = mesh_for_shape( silhouette, tilted_circles, symmetry_xsection, symmetry_xsection_front_to_back )
    mesh.write_OBJ( "sweep_shape-test4.obj" )

def test5():
    
    ## A CCW rotation in the xy plane by pi/2 radians.
    rot_mat = asarray( get_rotation_matrix( (0,0,1), pi/4 ) )
    
    silhouette = [ ( .75 * cos( theta ) + kEps, 1. * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] ]
    silhouette = [ dot( rot_mat[:2,:2], pt ) for pt in silhouette ]
    tilted_circles = [
        ( (-.5,.5,0.), .2, (-1.,1.,0.) ),
        ( (0.,0.,0.), .5, (-1.,1.,0.) ),
        ( (.5,-.5,0.), .2, (-1.,1.,0.) ),
        ]
    symmetry_xsection = [ ( .2 * cos( theta ) + kEps, 1. * sin( theta ) + kEps ) for theta in linspace( 0, 2*pi, 21 )[:-1] + pi/2 ]
    symmetry_xsection = [ dot( rot_mat[:2,:2], pt ) for pt in symmetry_xsection ]
    symmetry_xsection_front_to_back = len( symmetry_xsection ) // 2
    mesh = mesh_for_shape( silhouette, tilted_circles, symmetry_xsection, symmetry_xsection_front_to_back )
    mesh.write_OBJ( "sweep_shape-test5.obj" )

def main():
    test1()
    test1rot()
    test2()
    test2rot()
    test3()
    test4()
    test5()

if __name__ == "__main__": main()
