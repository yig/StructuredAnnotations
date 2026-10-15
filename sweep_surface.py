from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *
from vecutil import *
#from helpers import *

from pprint import pprint

kFixedStyleTangentMagSimple = True

## This controls the number of samples in each circular cross section.
NUM_CIRCLE_SAMPLES = None
NUM_SPHERE_CAP_SLICES = None
## This controls the size of steps when tracing points in the smoothed distance field.
kNoPinchTraceDt = None

kLowRes = None

kGenUVs = False

def set_lowres( lowres ):
    global NUM_CIRCLE_SAMPLES, NUM_SPHERE_CAP_SLICES, kNoPinchTraceDt, kLowRes
    
    kLowRes = bool( lowres )
    if kLowRes:
        NUM_CIRCLE_SAMPLES = 8
        NUM_SPHERE_CAP_SLICES = 2
        kNoPinchTraceDt = 6e-2
    else:
        NUM_CIRCLE_SAMPLES = 20
        NUM_SPHERE_CAP_SLICES = 5
        kNoPinchTraceDt = 3e-2

set_lowres( False )


def path_frames( path ):
    '''
    Given an open 3d path (sequence of 3d points) 'path',
    returns a tuple containing a frame (x,y,z unit vectors) for every vertex of the path.
    '''
    
    sweep_path = asarray( path )
    
    ## Use the forward difference as the tangent (and repeat the last one)
    sweep_tangents = sweep_path[1:] - sweep_path[:-1]
    sweep_tangents = sweep_tangents / sqrt( ( sweep_tangents**2 ).sum(1) )[:,newaxis]
    sweep_tangents = append( sweep_tangents, sweep_tangents[-1:], 0 )
    
    ## A frame of rotation
    ortho0 = orthogonal( sweep_tangents[0] )
    sweep_frames = [ ( ortho0, cross( sweep_tangents[0], ortho0 ), sweep_tangents[0] ) ]
    del ortho0
    for i in range( 1, len( sweep_path ) ):
        #rot = asarray( get_rotation_matrix( *get_axis_angle_rotation( sweep_tangents[i-1], sweep_tangents[i] ) ) )
        #sweep_frames.append( ( dot( rot, sweep_frames[-1][0] ), dot( rot, sweep_frames[-1][1] ), sweep_tangents[i] ) )
        axis, radians = get_axis_angle_rotation( sweep_tangents[i-1], sweep_tangents[i] )
        sweep_frames.append( ( rotate( sweep_frames[-1][0], axis, radians ), rotate( sweep_frames[-1][1], axis, radians ), sweep_tangents[i] ) )
    
    assert len( sweep_path ) == len( sweep_tangents )
    assert len( sweep_path ) == len( sweep_frames )
    
    #print 'path_frames non-orthogonality max abs:', abs( array( [ ( dot(x,y), dot(y,z), dot(z,x) ) for x,y,z in sweep_frames ] ) ).max()
    
    return sweep_frames

def generate_cap( radius_0, radius_0_prime, cap_style, cap_args, cap_segments ):
    '''
    Given a radius 'radius_0' and derivative of radius 'radius_0_prime' at x = 0
    and a cap style (described below) along with a list or tuple of style parameters 'cap_args',
    returns a sequence of ( x, radius ) tuples to append to the end of the path.
    The last element of the returned sequence is always guaranteed to have radius 0.
    
    Cap style must be one of:
        'flat': the mesh is closed at both ends with a flat circle.
        'hemisphere': the mesh is closed at both ends with a hemisphere.
        'spherical': the mesh is closed at both ends with a portion of a sphere
                     matching the gradients.
        'spherical-clipped': the mesh is closed at both ends with a portion of a sphere
                             matching the gradients, up to but not exceeding a hemisphere.
        'fixed': the mesh is closed smoothly.  The 'cap_args' must be a single floating point
                 value, the x at which the cap should end.
    '''
    
    if cap_style == 'flat':
        return [( 0., 0. )]
    
    elif cap_style in ('hemisphere', 'spherical', 'spherical-clipped'):
        ## Now adjust these values to add a piece of the sphere which matches
        ## the gradients, too.
        ## radius_0 is the y value at x = 0.
        ## radius_0_prime is the y' at x = 0.
        y = radius_0
        y_prime = radius_0_prime
        
        c,r = circle_from_point_and_gradient( 0., y, y_prime )
        assert r >= 0.
        assert r > abs(c)
        
        if cap_style == 'hemisphere' or ( cap_style == 'spherical-clipped' and c > 0 ):
            c = 0.
            r = y
        
        ## sqrt() our linear spacing so our output rings don't bunch up.
        result = []
        for t in sqrt( linspace( 0, 1, cap_segments+2 )[1:-1] ):
            x = (r+c) * t
            ri = sqrt( r**2 - (x-c)**2 )
            result.append( ( x, ri ) )
        result.append( ( r+c, 0. ) )
    
    elif cap_style == 'fixed':
        cap_radius = cap_args[0]
        assert float( cap_radius ) == cap_radius
        ## If there is no cap, just return the end-point.
        if cap_radius <= kEps:
            return [ ( cap_radius, 0. ) ]
        
        ## Now adjust these values to add a piece of the sphere which matches
        ## the gradients, too.
        ## radius_0 is the y value at x = 0.
        ## radius_0_prime is the y' at x = 0.
        ## The endpoint:
        y0 = asarray( ( 0., radius_0 ) )
        ## The tangent:
        y0_prime = asarray( ( 1., radius_0_prime ) )
        
        ## The tip of the cap.
        y1 = asarray( ( cap_radius, 0. ) )
        ## The cap tip's tangent is flat.
        y1_prime = asarray( ( 0., -1. ) )
        
        ## Now evaluate the bezier with the control points: ( y0, y0 + y0_prime, y1 - y1_prime, y1 ).
        ## The magnitude of y0_prime and y1_prime should be adjusted to give a nice curve.
        kFixedStyleTangentMag = 'medium'
        if kFixedStyleTangentMag == 'simple':
            tangent_mag = .5 * cap_radius
            y0_prime = tangent_mag * dir( y0_prime )
            y1_prime = tangent_mag * dir( y1_prime )
        
        elif kFixedStyleTangentMag == 'medium':
            ## Intersect the line passing through y0 with tangent y0_prime with the x = cap_radius
            ## line.
            ## Set y0_prime to the vector from y0 to the intersection point,
            ## and set y1_prime to the vector from y1 to the intersection point.
            #y0_prime *= cap_radius / y0_prime[0]
            y0_prime *= cap_radius
            #y1_prime *= -( y0[1] + y0_prime[1] ) / y1_prime[1]
            y1_prime[1] = -( y0[1] + y0_prime[1] )
            ## UPDATE: This produces boxy results.
            ## Q1: Halve the tangents?
            ## A1: The end-cap gets very sharp.
            y0_prime *= .5
            y1_prime *= .5
            ## Q2: Halve only y1_prime?
            ## A2: It's too boxy at default.
            #y1_prime *= .5
            
            ## Finally, clamp so that y0 + y0_prime doesn't intersect the x axis.
            if y1_prime[1] > 0.:
                y1_prime[1] = 0.
                y0_prime *= -y0[1] / y0_prime[1]
                cap_radius = y0_prime[0]
                y1[0] = cap_radius
        
        elif kFixedStyleTangentMag == 'fonts':
            ## Code from the fonts work which picks nice magnitudes.
            from helpers import tEquidistant
            y0_prime = dir( y0_prime )
            y1_prime = dir( y1_prime )
            #debugger()
            #print "AARG this doesn't work because we get a negative discriminant."
            t = tEquidistant( y0, y0_prime, y1, y1_prime )
            print(t)
            y0_prime = t * y0_prime
            y1_prime = t * y1_prime
        
        else:
            raise RuntimeError("Unknown font cap tangent magnitudes")
        
        
        from bezier import Bezier_Curve
        bez = Bezier_Curve( (
            y0,
            y0 + y0_prime,
            y1 - y1_prime,
            y1
            ) )
        result = []
        for t in linspace( 0, 1, cap_segments+2 )[1:-1]:
            p = bez( t )
            result.append( ( p[0], p[1] ) )
        result.append( ( cap_radius, 0. ) )
    
    else:
        raise RuntimeError("Unknown caps style: %s" % (caps,))
    
    assert 0. == result[-1][1]
    return result

def stitch_cross_sections( cross_sections ):
    '''
    Given a sequence of cross sections (closed loops of 3d vertices),
    returns an open mesh (homeomorphic to a cylinder open at both ends).
    
    The vertices in the mesh are exactly the flattened vertices in
    'cross_sections': the vertices of the first cross section come first,
    then the vertices of the second cross section, and so on.
    
    Assumes each cross section has the same number of vertices, and that
    the i-th vertex in each cross section should be attacked to the i-th
    vertex in the next and previous cross sections.
    
    Assumes each cross section is oriented counter-clockwise (CCW)
    when viewed from above, where above means that cross_sections[i]
    is closer to the viewer than cross_sections[i-1].
    '''
    
    from mesh_gui import TriMesh
    mesh = TriMesh()
    if kGenUVs: mesh.uvs = []
    
    ## Handle the degenerate case.
    if len( cross_sections ) == 0: return mesh
    
    ## All cross sections must have the same number of vertices.
    num_cross_sections = len( cross_sections )
    num_cross_section_verts = len( cross_sections[0] )
    assert len( set( map( len, cross_sections ) ) ) == 1
    '''
    for csi in xrange( num_cross_sections - 1 ):
        for vi in xrange( num_cross_section_verts ):
            v00 = csi * num_cross_section_verts + vi
            v10 = csi * num_cross_section_verts + ( (vi + 1) % num_cross_section_verts )
            v01 = (csi + 1) * num_cross_section_verts + vi
            v11 = (csi + 1) * num_cross_section_verts + ( (vi + 1) % num_cross_section_verts )
            mesh.faces.append( ( v00, v10, v11 ) )
            mesh.faces.append( ( v00, v11, v01 ) )
    '''
    from face_stitch_pattern_c import face_stitch_pattern
    mesh.faces = list(map( tuple, face_stitch_pattern( num_cross_sections, num_cross_section_verts ) ))
    
    for v, cross_section in enumerate( cross_sections ):
        mesh.vs.extend( list(map( tuple, cross_section )) )
        if kGenUVs:
            mesh.uvs.extend( [
                ( u, float(v+1) / ( len( cross_sections ) + 1 ) )
                for u in linspace( 0, 1, len( cross_section )+1 )[::-1]
                ] )
    
    return mesh

def path_and_radii_with_caps( sweep_path, sweep_radii, caps = (('spherical-clipped',),('spherical-clipped',)), cap_segments = 5 ):
    '''
    Given a 3d path (sequence of 3d points) 'sweep_path',
    a corresponding sequence of radii along the path (1 radius for each vertex),
    a parameter 'caps' specifying the type of start- and end-cap as a 2-tuple,
    and a parameter 'cap_segments' specifying the number of line segments each cap should add,
    returns a tuple of four items ( path, radii, start cap point, end cap point ),
    where path and radii are modified 'sweep_path' and 'sweep_radii' with end-caps added
    in such a way that attaching circles to the path with corresponding radii gives the
    desired end-caps, and start/end cap point are the terminating points (radius = 0)
    at the beginning and end of the path.
    
    'caps' is a 2-tuple where the first element specifies the type of start cap and the
    second element specifies the type of end cap.  Each specification must be a tuple or list
    containing the name of the cap style followed by zero-or-more values, depending on the cap style.
    Valid cap types are those that can be passed to generate_cap() as the 'cap_style' parameter.
    
    Optional parameter cap_segments determines the number of start/end-cap samples.
    '''
    
    ## Attaching hemispherical start and end caps is best accomplished by
    ## adding points to the beginning and end of the sweep_path in the tangent direction
    ## with decreasing radii.
    
    ## We need two points for a line.
    assert len( sweep_path ) >= 2
    assert len( sweep_path ) == len( sweep_radii )
    ## We want lists of arrays/floats, and copies, because we will insert elements for caps.
    sweep_path = list(map( arrayf, sweep_path ))
    sweep_radii = list(map( float, sweep_radii ))
    
    ## We can't have negative radii.
    assert False not in [ radius > -kEps for radius in sweep_radii ]
    
    ## cap_segments should be a non-negative integer.
    assert int( cap_segments ) == cap_segments
    assert cap_segments >= 0
    
    assert len( caps ) == 2
    start_cap = caps[0]
    end_cap = caps[1]
    
    ## Let's also handle the case where a radius of 0 is passed for the first point.
    if sweep_radii[0] < kEps:
        start_cap_pt = sweep_path[0]
        sweep_path = sweep_path[1:]
        sweep_radii = sweep_radii[1:]
    
    else:
        ## sweep_radii[0] is the y value at x = 0.
        ## y' is ( sweep_radii[0] - sweep_radii[1] ) / mag( sweep_path[0] - sweep_path[1] )
        y = sweep_radii[0]
        y_prime = ( sweep_radii[0] - sweep_radii[1] ) / mag( sweep_path[0] - sweep_path[1] )
        
        x_rs = generate_cap( y, y_prime, start_cap[0], start_cap[1:], cap_segments )
        assert x_rs[-1][1] == 0.
        
        axis = dir( sweep_path[0] - sweep_path[1] )
        origin = sweep_path[0]
        for x, r in x_rs[:-1]:
            sweep_path.insert( 0, origin + x * axis )
            sweep_radii.insert( 0, r )
        
        start_cap_pt = origin + x_rs[-1][0] * axis
    
    
    ## Let's also handle the case where a radius of 0 is passed for the last point.
    if sweep_radii[-1] < kEps:
        end_cap_pt = sweep_path[-1]
        sweep_path = sweep_path[:-1]
        sweep_radii = sweep_radii[:-1]
    
    else:
        ## sweep_radii[-1] is the y value at x = 0.
        ## y' is ( sweep_radii[-1] - sweep_radii[-2] ) / mag( sweep_path[-1] - sweep_path[-2] )
        y = sweep_radii[-1]
        y_prime = ( sweep_radii[-1] - sweep_radii[-2] ) / mag( sweep_path[-1] - sweep_path[-2] )
        
        x_rs = generate_cap( y, y_prime, end_cap[0], end_cap[1:], cap_segments )
        assert x_rs[-1][1] == 0.
        
        axis = dir( sweep_path[-1] - sweep_path[-2] )
        origin = sweep_path[-1]
        for x, r in x_rs[:-1]:
            sweep_path.append( origin + x * axis )
            sweep_radii.append( r )
        
        end_cap_pt = origin + x_rs[-1][0] * axis
    
    return sweep_path, sweep_radii, start_cap_pt, end_cap_pt

def path_and_cross_sections_with_caps( sweep_path, sweep_cross_sections, caps = (('spherical-clipped',),('spherical-clipped',)), cap_segments = 5 ):
    '''
    Given a 3d path (sequence of 3d points) 'sweep_path',
    a corresponding sequence of cross sections 'sweep_cross_sections' along the path, and
    a parameter 'caps' specifying the type of start- and end-cap as a 2-tuple,
    returns a tuple of four items ( path, cross_sections, start cap point, end cap point ),
    where path and cross_sections are modified 'sweep_path' and 'sweep_cross_sections' with
    end-caps added in such a way that attaching cross sections to the path as given yields
    the desired end-caps, and start/end cap point are the terminating points (radius = 0)
    at the beginning and end of the path.
    
    Each cross section is defined by the tuple ( frame x-vector, 2d line loop ),
    where the 2d line loop is defined to lie on the
        'frame x-vector', cross( path tangent, frame x-vector ) plane
    and the 2d line loops are correspondingly sampled.
    
    'caps' is a 2-tuple where the first element specifies the type of start cap and the
    second element specifies the type of end cap.  Each specification must be a tuple or list
    containing the name of the cap style followed by zero-or-more values, depending on the cap style.
    Valid cap types are those that can be passed to generate_cap() as the 'cap_style' parameter.
    
    Optional parameter cap_segments determines the number of start/end-cap samples.
    
    lightly tested
    '''
    
    ## Attaching hemispherical start and end caps is best accomplished by
    ## adding points to the beginning and end of the sweep_path in the tangent direction
    ## with decreasing circular radii.
    
    ## We need two points for a line.
    assert len( sweep_path ) >= 2
    ## We want one cross section per sweep path sample.  (This is a low level function.)
    assert len( sweep_path ) == len( sweep_cross_sections )
    ## Frame vectors are 3d.
    assert False not in [ 3 == len( framevec ) for framevec, xsection in sweep_cross_sections ]
    ## Cross sections must have the number of samples.
    assert 1 == len( set([ len( xsection ) for framevec, xsection in sweep_cross_sections ]) )
    ## Cross sections are line loops, so they must have at least 3 points.
    assert len( sweep_cross_sections[0][1] ) >= 3
    
    ## We want lists of arrays/floats, and copies, because we will insert elements for caps.
    sweep_path = list(map( arrayf, sweep_path ))
    sweep_cross_sections = [ ( asarray( framevec ), asarray( xsection ) ) for framevec, xsection in sweep_cross_sections ]
    
    ## cap_segments should be a non-negative integer.
    assert int( cap_segments ) == cap_segments
    assert cap_segments >= 0
    
    assert len( caps ) == 2
    start_cap = caps[0]
    end_cap = caps[1]
    
    ## Let's also handle the case where a single point is passed for the first cross section.
    if len( sweep_cross_sections[0][1] ) == 1:
        start_cap_pt = sweep_cross_sections[0][1]
        sweep_path = sweep_path[1:]
        sweep_cross_sections = sweep_cross_sections[1:]
    
    else:
        sweep_radii_0 = average( list(map( mag, sweep_cross_sections[0][1] )) )
        sweep_radii_1 = average( list(map( mag, sweep_cross_sections[1][1] )) )
        
        ## sweep_radii[0] is the y value at x = 0.
        ## y' is ( sweep_radii[0] - sweep_radii[1] ) / mag( sweep_path[0] - sweep_path[1] )
        y = sweep_radii_0
        y_prime = ( sweep_radii_0 - sweep_radii_1 ) / mag( sweep_path[0] - sweep_path[1] )
        
        x_rs = generate_cap( y, y_prime, start_cap[0], start_cap[1:], cap_segments )
        assert x_rs[-1][1] == 0.
        
        axis = dir( sweep_path[0] - sweep_path[1] )
        origin = sweep_path[0]
        ## lerp() each cross section towards a perfect circle with given radius.
        sweep_cross_sections_0 = sweep_cross_sections[0]
        for x, ri in x_rs[:-1]:
            t = x / x_rs[-1][0]
            sweep_path.insert( 0, origin + x * axis )
            sweep_cross_sections.insert( 0, ( sweep_cross_sections_0[0], [
                ## TODO: lerp to an actual circle---don't just project each point onto a circle---
                ##       so that a cross section with fold-overs doesn't lead to self-intersections.
                ##       Another way I can do this is by normalizing the cross section,
                ##       then placing it and a unit circle on either ends of a line with
                ##       cap_segments values, and diffusing values between.  Hmm, since there's only
                ##       two samples for each value, it will be the same as a lerp.
                # ri*dir( pt )
                # lerp( pt, ri*dir( pt ), t )
                lerp( (ri/sweep_radii_0)*pt, ri*dir( pt ), t )
                ## This doesn't work:
                # ri * pt * lerp( 1./sweep_radii_0, 1./len( pt ), t )
                for pt in sweep_cross_sections_0[1]
                ] ) )
        start_cap_pt = origin + x_rs[-1][0] * axis
    
    
    ## Let's also handle the case where a single point is passed for the last cross section.
    if len( sweep_cross_sections[-1][1] ) == 1:
        end_cap_pt = sweep_cross_sections[-1][1]
        sweep_path = sweep_path[:-1]
        sweep_cross_sections = sweep_cross_sections[:-1]
    
    else:
        sweep_radii_n1 = average( list(map( mag, sweep_cross_sections[-1][1] )) )
        sweep_radii_n2 = average( list(map( mag, sweep_cross_sections[-2][1] )) )
        
        ## sweep_radii[-1] is the y value at x = 0.
        ## y' is ( sweep_radii[-1] - sweep_radii[-2] ) / mag( sweep_path[-1] - sweep_path[-2] )
        y = sweep_radii_n1
        y_prime = ( sweep_radii_n1 - sweep_radii_n2 ) / mag( sweep_path[-1] - sweep_path[-2] )
        
        x_rs = generate_cap( y, y_prime, end_cap[0], end_cap[1:], cap_segments )
        assert x_rs[-1][1] == 0.
        
        axis = dir( sweep_path[-1] - sweep_path[-2] )
        origin = sweep_path[-1]
        ## lerp() each cross section towards a perfect circle with given radius.
        sweep_cross_sections_n1 = sweep_cross_sections[-1]
        for x, ri in x_rs[:-1]:
            t = x / x_rs[-1][0]
            sweep_path.append( origin + x * axis )
            sweep_cross_sections.append( ( sweep_cross_sections_n1[0], [
                ## TODO: lerp to an actual circle---don't just project each point onto a circle---
                ##       so that a cross section with fold-overs doesn't lead to self-intersections.
                ##       Another way I can do this is by normalizing the cross section,
                ##       then placing it and a unit circle on either ends of a line with
                ##       cap_segments values, and diffusing values between.  Hmm, since there's only
                ##       two samples for each value, it will be the same as a lerp.
                # ri*dir( pt )
                # lerp( pt, ri*dir( pt ), t )
                lerp( (ri/sweep_radii_n1)*pt, ri*dir( pt ), t )
                ## This doesn't work:
                # ri * pt * lerp( 1./sweep_radii_n1, 1./len( pt ), t )
                for pt in sweep_cross_sections_n1[1]
                ] ) )
        end_cap_pt = origin + x_rs[-1][0] * axis
    
    
    return sweep_path, sweep_cross_sections, start_cap_pt, end_cap_pt

def mesh_from_sweep_circles( sweep_path, sweep_radii, caps = (('spherical-clipped',),('spherical-clipped',)) ):
    '''
    Given a 3d path (sequence of 3d points) 'sweep_path' and a corresponding
    sequence of radii along the path (1 radius for each vertex),
    returns a TriMesh object created by approximating the surface defined
    by the union of spheres placed at each vertex.
    
    Optional parameter 'caps' determines the end caps, and may be any value accepted by
    path_and_radii_with_caps().
    '''
    
    ## We need three points (and no two tangents to be equal) so that we have curve normals.
    assert len( sweep_path ) >= 3
    assert len( sweep_path ) == len( sweep_radii )
    
    ## 1 Build a sequence of cross sections.
    ## 2 Call stitch_cross_sections().
    ## 3 Attach hemisphere caps to the beginning and end of the stitched mesh.
    
    ## 3
    ## Actually, attaching hemispherical start and end caps is best accomplished by
    ## adding points to the beginning and end of the sweep_path in the tangent direction
    ## with decreasing radii.
    sweep_path, sweep_radii, start_cap_pt, end_cap_pt = path_and_radii_with_caps(
        sweep_path, sweep_radii, caps, cap_segments = NUM_SPHERE_CAP_SLICES
        )
    
    ## We can't have negative radii.
    assert False not in [ radius > -kEps for radius in sweep_radii ]
    
    
    ## 1
    sweep_frames = path_frames( sweep_path )
    assert len( sweep_path ) == len( sweep_frames )
    
    '''
    print "sweep_path:"
    pprint( sweep_path )
    
    print "sweep_frames:"
    pprint( sweep_frames )
    '''
    
    cross_sections = []
    for ( frame0, frame1, frame2 ), ( center, radius ) in zip( sweep_frames, zip( sweep_path, sweep_radii ) ):
        cross_section = []
        
        ## Sample the circular cross section.
        for theta in linspace( 0, 2*pi, NUM_CIRCLE_SAMPLES+1 )[:-1]:
            cross_section.append( center + radius * cos( theta ) * frame0 + radius * sin( theta ) * frame1 )
        
        cross_sections.append( cross_section )
    
    assert len( cross_sections ) == len( sweep_path )
    
    
    ## 2
    mesh = stitch_cross_sections( cross_sections )
    
    
    ## 3
    ## The start cap.  Attach all points in the first cross section to start_cap_pt.
    mesh.vs.append( tuple( start_cap_pt ) )
    if kGenUVs: mesh.uvs.append( ( 0.,0. ) )
    for i in range( NUM_CIRCLE_SAMPLES ):
        mesh.faces.append( ( i, len( mesh.vs )-1, (i+1) % NUM_CIRCLE_SAMPLES ) )
    
    ## The end cap.  Attach all points in the last cross section to end_cap_pt.
    mesh.vs.append( tuple( end_cap_pt ) )
    if kGenUVs: mesh.uvs.append( ( 0.,1. ) )
    offset = ( len( cross_sections )-1 ) * NUM_CIRCLE_SAMPLES
    for i in range( NUM_CIRCLE_SAMPLES ):
        mesh.faces.append( ( len( mesh.vs )-1, offset + i, offset + ( (i+1) % NUM_CIRCLE_SAMPLES ) ) )
    
    mesh.topology_changed()
    
    return mesh

def mesh_from_sweep_circles_no_pinch( sweep_path, sweep_radii, caps = (('spherical-clipped',),('spherical-clipped',)) ):
    '''
    Given a 3d path (sequence of 3d points) 'sweep_path' and a corresponding
    sequence of radii along the path (1 radius for each vertex),
    returns a TriMesh object created by approximating the surface defined
    by the union of spheres placed at each vertex.
    
    Optional parameter 'caps' determines the end caps, and may be any value accepted by
    path_and_radii_with_caps().
    '''
    
    ## We need three points (and no two tangents to be equal) so that we have curve normals.
    assert len( sweep_path ) >= 3
    assert len( sweep_path ) == len( sweep_radii )
    
    ## 1 Build a sequence of cross sections.
    ## 2 Call stitch_cross_sections().
    ## 3 Attach hemisphere caps to the beginning and end of the stitched mesh.
    
    ## 3
    ## Actually, attaching hemispherical start and end caps is best accomplished by
    ## adding points to the beginning and end of the sweep_path in the tangent direction
    ## with decreasing radii.
    sweep_path, sweep_radii, start_cap_pt, end_cap_pt = path_and_radii_with_caps(
        sweep_path, sweep_radii, caps, cap_segments = NUM_SPHERE_CAP_SLICES
        )
    
    ## We can't have negative radii.
    assert False not in [ radius > -kEps for radius in sweep_radii ]
    
    
    ## 1
    sweep_frames = path_frames( sweep_path )
    assert len( sweep_path ) == len( sweep_frames )
    
    '''
    print "sweep_path:"
    pprint( sweep_path )
    
    print "sweep_frames:"
    pprint( sweep_frames )
    '''
    
    ## Initialize the cross sections to be offset only by trace_eps from the spine.
    trace_eps = 1e-3
    circle_samples = linspace( 0, 2*pi, NUM_CIRCLE_SAMPLES+1 )[:-1][...,newaxis]
    cross_sections = emptyf( ( len( sweep_frames ), len( circle_samples ), 3 ) )
    
    off = 0
    for ( frame0, frame1, frame2 ), center in zip( sweep_frames, sweep_path ):
    #for off in xrange( len( sweep_frames ) ):
    #    frame0, frame1, frame2 = sweep_frames[off]
    #    center = sweep_path[off]
        
        ## Sample the circular cross section.
        '''
        theta_i = 0
        for theta in circle_samples:
            cross_sections[ off, theta_i ] = center + trace_eps * cos( theta ) * frame0 + trace_eps * sin( theta ) * frame1
            theta_i += 1
        '''
        
        cross_sections[ off ] = center + trace_eps * ( cos( circle_samples ) * frame0 + sin( circle_samples ) * frame1 )
        
        #cross_sections[ off ][:] = cross_section
        off += 1
    
    assert len( cross_sections ) == len( sweep_path )
    
    
    ## 2
    mesh = stitch_cross_sections( cross_sections )
    
    ## 3
    ## The start cap.  Attach all points in the first cross section to start_cap_pt.
    mesh.vs.append( tuple( start_cap_pt ) )
    if kGenUVs: mesh.uvs.append( ( 0.,0. ) )
    for i in range( NUM_CIRCLE_SAMPLES ):
        mesh.faces.append( ( i, len( mesh.vs )-1, (i+1) % NUM_CIRCLE_SAMPLES ) )
    
    ## The end cap.  Attach all points in the last cross section to end_cap_pt.
    mesh.vs.append( tuple( end_cap_pt ) )
    if kGenUVs: mesh.uvs.append( ( 0.,1. ) )
    offset = ( len( cross_sections )-1 ) * NUM_CIRCLE_SAMPLES
    for i in range( NUM_CIRCLE_SAMPLES ):
        mesh.faces.append( ( len( mesh.vs )-1, offset + i, offset + ( (i+1) % NUM_CIRCLE_SAMPLES ) ) )
    
    ## TODO Q: Why do the caps flatten out?
    ## A1: It's because all other points are moving according to the smoothed distance field,
    ##     which pushes them ever-so-slightly along the spine.  I tried to push the end points
    ##     along the smoothed distance field, too, but the same thing happens unless it flows 2x
    ##     what it should.
    ## Now trace all points away from the spine.
    sweep_path = asarray( sweep_path )
    
    start_cap_vec = asarray( start_cap_pt ) - sweep_path[0]
    end_cap_vec = asarray( end_cap_pt ) - sweep_path[-1]
    
    start_cap_vec_mag = mag( start_cap_vec )
    end_cap_vec_mag = mag( end_cap_vec )
    
    start_cap_vec_dir = start_cap_vec
    if start_cap_vec_mag > kEps: start_cap_vec_dir = start_cap_vec_dir / start_cap_vec_mag
    end_cap_vec_dir = end_cap_vec
    if end_cap_vec_mag > kEps: end_cap_vec_dir = end_cap_vec_dir / end_cap_vec_mag
    
    mesh.vs[-2] = sweep_path[0] + trace_eps * start_cap_vec_dir
    mesh.vs[-1] = sweep_path[-1] + trace_eps * end_cap_vec_dir
    
    ## Slow python code
    '''
    from edge_distances import trace_pts_in_smooth_distance_field
    mesh.vs = trace_pts_in_smooth_distance_field(
        asarray( mesh.vs ).T,
        asarray( asarray(
            [ [radius] * len( cross_section ) for radius, cross_section in izip( sweep_radii, cross_sections ) ]
            ).ravel().tolist() + [ 2*mag( start_cap_vec ), 2*mag( end_cap_vec ) ] ),
        zip( sweep_path[:-1], sweep_path[1:] ),
        dt = 3e-2
        ).T
    ## We're not done adding points to the mesh, so convert vertices back to a list.
    mesh.vs = list( mesh.vs )
    '''
    
    ## Faster C++ code.
    from edge_distances_cpp import trace_point
    kEndCapMultiplier = 1
    '''
    distances = asarray( asarray(
            [ [radius] * len( cross_section ) for radius, cross_section in izip( sweep_radii, cross_sections ) ]
            ).ravel().tolist() + [ kEndCapMultiplier*start_cap_vec_mag, kEndCapMultiplier*end_cap_vec_mag ] )
    '''
    distances = onesf( len(sweep_radii) * len(circle_samples) + 2 )
    distances[ : len( distances )-2 ] = ( asarray( sweep_radii )[...,newaxis].repeat( len(circle_samples), 1 ) ).ravel()
    distances[-2] = kEndCapMultiplier*start_cap_vec_mag
    distances[-1] = kEndCapMultiplier*end_cap_vec_mag
    
    '''
    path_indices = asarray( asarray(
            [ [i] * len( cross_section ) for i, cross_section in enumerate( cross_sections ) ],
            dtype = int
            ).ravel().tolist() + [ 0, len( sweep_path )-1 ], dtype = int )
    assert len( mesh.vs ) == len( path_indices )
    kSpan = 100
    '''
    #print "+++ trace_point()'s"
    for i in range(len( mesh.vs )):
        #print path_indices[i]
        #mesh.vs[i] = trace_point( mesh.vs[i], sweep_path[ max( 0, path_indices[i] - kSpan ) : path_indices[i] + kSpan+1 ], distances[i], dt = 3e-2 )
        mesh.vs[i] = trace_point( mesh.vs[i], sweep_path, distances[i], dt = 3e-2 )
    #print "--- trace_point()'s"
    
    mesh.topology_changed()
    
    return mesh

def smooth_orthogonal_vectors_from_constraints( sweep_path, sparse_orthogonal_vectors, sweep_path_frames = None, path_diffuser = None ):
    '''
    Given a 3d path (sequence of N 3d points) 'sweep_path', and
    a sequence of sparse orthogonal vectors 'sparse_orthogonal_vectors' along the path,
    returns an N-by-3 array of dense, smoothly varying orthogonal vectors along the path.
    
    Each sparse orthogonal constraints is a tuple
    ( index of point in 'sweep_path', 3d vector orthogonal to sweep_path ).
    
    If optional parameter 'sweep_path_frames' is specified, it must be a sequence of N tuples,
    each of which contains the 3d x,y,z vectors corresponding to the frame along 'sweep_path'.
    By default, a minimally rotating frame is computed.
    
    If optional parameter 'path_diffuser' is specified, it must be the result of
    diffuse_values_along_sweep_spine_factored( sweep_path ).
    
    lightly tested/stepped-through
    '''
    
    ## We need three points (and no two tangents to be equal) so that we have curve normals.
    assert len( sweep_path ) >= 3
    ## We want all point indices to be valid and positive-indexed.
    assert False not in [ p_i >= 0 and p_i < len( sweep_path ) for p_i, xvec in sparse_orthogonal_vectors ]
    ## Frame vectors are 3d.
    assert False not in [ 3 == len( framevec ) for p_i, framevec in sparse_orthogonal_vectors ]
    
    ## Given sparse orthogonal vectors, our task is to interpolate them at every point along
    ## 'sweep_path'.
    ## We can slerp them, or we can convert them to a single parameter along
    ## the line's frame and call sweep_spine.diffuse_values_along_sweep_spine() on those.
    ## (The single parameter will be an angle, so I don't think this will actually work.
    ## On the other hand, if the frame vectors come from the symmetry sheet, they will all be
    ## within pi of each other relative to a fixed world-space vector, but perhaps not relative
    ## to a vector on the line's frame.
    ## Actually, I can just add multiples of 2*pi to the angles so that each angle is less than pi
    ## from its predecessor.)
    ## 
    ## 1 Generate a frame along the spine.
    ## 2 Convert frame vectors into angles on the frame, each less than pi from its predecessor.
    ## 3 Call sweep_spine.diffuse_values_along_sweep_spine on these angles.
    ## 4 Convert the angles back into frame vectors.
    
    ## 1
    if sweep_path_frames is not None:
        sweep_frames = [ tuple( map( asarray, xyz ) ) for xyz in sweep_path_frames ]
    else:
        sweep_frames = path_frames( sweep_path )
    assert len( sweep_path ) == len( sweep_frames )
    
    ## 2
    sparse_orthogonal_vecs = [ dir( asarray( framevec ) ) for p_i, framevec in sparse_orthogonal_vectors ]
    assert len( sparse_orthogonal_vecs ) == len( sparse_orthogonal_vectors )
    sparse_frame_angles = [
        atan2( dot( framevec, y ), dot( framevec, x ) )
        for framevec, (x,y,z) in zip(
            sparse_orthogonal_vecs,
            [ sweep_frames[ p_i ] for p_i, framevec in sparse_orthogonal_vectors ]
            )
        ]
    for i in range( 1, len( sparse_frame_angles ) ):
        last_angle = sparse_frame_angles[i-1]
        cur_angle = sparse_frame_angles[i]
        
        ## Add k*2*pi such that 0 <= last_angle - cur_angle < 2*pi.
        ## Then, optionally add one more 2*pi so that | last_angle - cur_angle | < pi.
        
        ## Hmm, there are a lot of ways to get this wrong.  I'll do something safe.
        cur_angle += 2*pi * floor( (last_angle - cur_angle) / (2*pi) )
        assert cur_angle <= last_angle
        while abs( last_angle - cur_angle ) > pi:
            cur_angle += 2*pi
            
            assert cur_angle < last_angle + 2*pi
        assert abs( last_angle - cur_angle ) <= pi
        
        sparse_frame_angles[i] = cur_angle
    
    ## 3
    if path_diffuser is None:
        import sweep_spine
        path_diffuser = sweep_spine.diffuse_values_along_sweep_spine_factored( sweep_path )
    
    assert len( sparse_orthogonal_vectors ) == len( sparse_frame_angles )
    dense_frame_angles = path_diffuser(
        [
            ( p_i, [ angle ] )
            for (p_i, framevec), angle
            in zip( sparse_orthogonal_vectors, sparse_frame_angles )
        ]
        ).ravel()
    
    ## 4
    assert len( dense_frame_angles ) == len( sweep_frames )
    dense_orthogonal_vecs = [
        cos( angle ) * x + sin( angle ) * y
        for angle, (x,y,z) in zip( dense_frame_angles, sweep_frames )
        ]
    return dense_orthogonal_vecs

def mesh_from_sparse_sweep_cross_sections( sweep_path, frame_xvecs, sparse_sweep_cross_sections, caps = 'spherical-clipped', path_diffuser = None ):
    '''
    Given a 3d path (sequence of 3d points) 'sweep_path',
    a sequence of 3d vectors corresponding to frame x-vectors along the path 'frame_xvecs', and
    a sparse set of cross sections 'sparse_sweep_cross_sections' along the path,
    returns a TriMesh object created by connecting each cross section.
    
    Each sparse cross section is defined by the tuple ( index of point in 'sweep_path', 2d line loop ),
    where the 2d line loop is defined to lie on the
        frame x-vector, cross( path tangent at point in 'sweep_path', frame x-vector ) plane
    and all 2d line loops are correspondingly sampled.
    
    Optional parameter 'caps' determines the end caps, and may be any value accepted by
    path_and_radii_with_caps().
    
    If optional parameter 'path_diffuser' is specified, it must be the result of
    diffuse_values_along_sweep_spine_factored( sweep_path ).
    
    lightly tested/stepped-through
    '''
    
    ## We need three points (and no two tangents to be equal) so that we have curve normals.
    assert len( sweep_path ) >= 3
    ## The path is 3d.
    assert False not in [ 3 == len( pt ) for pt in sweep_path ]
    ## We should have the same number of frame x-vectors as points in the path.
    assert len( sweep_path ) == len( frame_xvecs )
    ## Frame vectors are 3d.
    assert False not in [ 3 == len( framevec ) for framevec in frame_xvecs ]
    ## We want all cross section point indices to be valid and positive-indexed.
    assert False not in [ p_i >= 0 and p_i < len( sweep_path ) for p_i, ll in sparse_sweep_cross_sections ]
    ## Cross sections must have the number of samples.
    assert 1 == len( set([ len( xsection ) for p_i, xsection in sparse_sweep_cross_sections ]) )
    ## Cross sections are line loops, so they must have at least 3 points.
    assert len( sparse_sweep_cross_sections[0][1] ) >= 3
    
    ## Given sparse cross sections, our task is to interpolate them at every point along
    ## 'sweep_path'.
    ## We can lerp the cross sections' points, or we can call
    ## sweep_spine.diffuse_values_along_sweep_spine() for each sample in the cross sections.
    ## (Calling sweep_spine.diffuse_values_along_sweep_spine() means more code re-use.)
    ## 
    ## 1 For each cross section parameter (index), call
    ##   sweep_spine.diffuse_values_along_sweep_spine() on the samples in all cross sections.
    ## alt 1 Build a mesh with circles at non-specified cross sections, and then do a global
    ##       laplacian-minimizing solve with the end-caps and sparse cross sections fixed.
    ## 2 Return result of dense function mesh_from_sweep_cross_sections()
    
    ## 1
    if path_diffuser is None:
        import sweep_spine
        path_diffuser = sweep_spine.diffuse_values_along_sweep_spine_factored( sweep_path )
    
    num_cross_section_samples = len( sparse_sweep_cross_sections[0][1] )
    dense_xsections = path_diffuser(
        [ ( p_i, asarray( xsection ).ravel() ) for p_i, xsection in sparse_sweep_cross_sections ]
        )
    dense_xsections.reshape( ( dense_xsections.shape[0], num_cross_section_samples, 2 ) )
    
    ## 2
    assert len( frame_xvecs ) == len( dense_xsections )
    mesh = mesh_from_sweep_cross_sections( sweep_path, list(zip( frame_xvecs, dense_xsections )), caps )
    return mesh

def mesh_from_sweep_cross_sections( sweep_path, sweep_cross_sections, caps = (('spherical-clipped',),('spherical-clipped',)), no_pinch = False ):
    '''
    Given a 3d path (sequence of 3d points) 'sweep_path', and
    a corresponding sequence of cross sections 'sweep_cross_sections' along the path,
    returns a TriMesh object created by connecting each cross section.
    Each cross section is defined by the tuple ( frame x-vector, 2d line loop ),
    where the 2d line loop is defined to lie on the
        'frame x-vector', cross( path tangent, frame x-vector ) plane
    and the 2d line loops are correspondingly sampled.
    
    Optional parameter 'caps' determines the end caps, and may be any value accepted by
    path_and_radii_with_caps().
    
    Optional parameter 'no_pinch' determines whether or not to try to keep the shape from
    self-intersecting.
    
    lightly tested
    '''
    
    ## We need three points (and no two tangents to be equal) so that we have curve normals.
    assert len( sweep_path ) >= 3
    ## We want one cross section per sweep path sample.  (This is a low level function.)
    assert len( sweep_path ) == len( sweep_cross_sections )
    ## Frame vectors are 3d.
    assert False not in [ 3 == len( framevec ) for framevec, xsection in sweep_cross_sections ]
    ## Cross sections must have the number of samples.
    assert 1 == len( set([ len( xsection ) for framevec, xsection in sweep_cross_sections ]) )
    ## Cross sections are line loops, so they must have at least 3 points.
    assert len( sweep_cross_sections[0][1] ) >= 3
    
    ## 1 Call stitch_cross_sections().
    ## 2 Attach hemisphere caps to the beginning and end of the stitched mesh.
    
    ## 2
    ## Actually, attaching hemispherical start and end caps is best accomplished by
    ## adding points to the beginning and end of the sweep_path in the tangent direction
    ## with decreasing radii.
    sweep_path, sweep_cross_sections, start_cap_pt, end_cap_pt = path_and_cross_sections_with_caps(
        sweep_path, sweep_cross_sections, caps, cap_segments = NUM_SPHERE_CAP_SLICES
        )
    
    ## 1
    ## Use the forward difference as the tangent (and repeat the last one)
    sweep_path = asarray( sweep_path )
    sweep_tangents = sweep_path[1:] - sweep_path[:-1]
    sweep_tangents = sweep_tangents / sqrt( ( sweep_tangents**2 ).sum(1) )[:,newaxis]
    sweep_tangents = append( sweep_tangents, sweep_tangents[-1:], 0 )
    
    ## Build 3d frames.
    ## TODO Q: Should I use the cross product of the cross product for the x-direction, instead
    ##         of trusting that framexvec is actually orthogonal to the tangent?
    ## A1: Allow people to input whatever they want.  Anyone using the
    ##     mesh_from_sparse_sweep_cross_sections() interface is guaranteed this.
    #print 'dot( xsection x, tangent ):', [ dot( framexvec, tangent ) for ( framexvec, xsection ), tangent in izip( sweep_cross_sections, sweep_tangents ) ]
    sweep_frames = [
        ( dir( asarrayf( framexvec ) ), dir( cross( tangent, framexvec ) ) )
        for ( framexvec, xsection ), tangent
        in zip( sweep_cross_sections, sweep_tangents )
        ]
    
    assert len( sweep_path ) == len( sweep_tangents )
    assert len( sweep_path ) == len( sweep_frames )
    
    '''
    print "sweep_path:"
    pprint( sweep_path )
    
    print "sweep_tangents:"
    pprint( sweep_tangents )
    
    print "sweep_frames:"
    pprint( sweep_frames )
    '''
    
    cross_sections_3d = []
    for pt, ( frame0, frame1 ), ( xxx, xsection ) in zip( sweep_path, sweep_frames, sweep_cross_sections ):
        
        xsection = asarrayf( xsection )
        cross_sections_3d.append(
            pt + xsection[:,:1] * frame0 + xsection[:,1:] * frame1
            )
    
    assert len( cross_sections_3d ) == len( sweep_path )
    
    
    if no_pinch:
        trace_eps = 1e-3
        
        from edge_distances_cpp import trace_point
        
        #print "+++ trace_point()'s"
        for spine_point, xsection in zip( sweep_path, cross_sections_3d ):
            for i in range(len( xsection )):
                off = xsection[i] - spine_point
                off_mag = mag( off )
                xsection[i] = trace_point( spine_point + trace_eps * ( (1./off_mag) * off ), sweep_path, off_mag, dt = kNoPinchTraceDt )
        #print "--- trace_point()'s"
    
    
    ## 2
    mesh = stitch_cross_sections( cross_sections_3d )
    
    
    ## 3
    ## I know that we are assuming that every cross section has the same number of samples,
    ## but I'm not building the assumption in when it is not needed.
    ## This stitching-to-the-end-cap code at least can handle differently samples cross sections.
    start_xsection_samples = len( sweep_cross_sections[0][1] )
    end_xsection_samples = len( sweep_cross_sections[-1][1] )
    end_cap_rim_offset = len( mesh.vs ) - len( sweep_cross_sections[-1][1] )
    
    ## The start cap.  Attach all points in the first cross section to start_cap_pt.
    mesh.vs.append( tuple( start_cap_pt ) )
    if kGenUVs: mesh.uvs.append( ( 0.,0. ) )
    for i in range( start_xsection_samples ):
        mesh.faces.append( ( i, len( mesh.vs )-1, (i+1) % start_xsection_samples ) )
    
    ## The end cap.  Attach all points in the last cross section to end_cap_pt.
    mesh.vs.append( tuple( end_cap_pt ) )
    if kGenUVs: mesh.uvs.append( ( 0.,1. ) )
    for i in range( len( sweep_cross_sections[-1][1] ) ):
        mesh.faces.append( ( len( mesh.vs )-1, end_cap_rim_offset + i, end_cap_rim_offset + ( (i+1) % end_xsection_samples ) ) )
    
    mesh.topology_changed()
    
    return mesh

def circle_from_point_and_gradient( x, y, y_prime ):
    '''
    Given a point (x,y) on a circle, with dy/dx of 'y_prime', returns
    the tuple (c, r) for the corresponding circle centered
    at (c,0) --- on the x-axis --- with radius r.
    '''
    
    c = x + y*y_prime
    r = y * sqrt( 1 + y_prime**2 )
    
    return c, r

def test():
    #sweep_path = [ (0,0,0), (1,0,0) ]
    #sweep_radii = [ 1., 1. ]
    sweep_path = [ (0,0,0), (.5,.1,0), (1,0,0), (1.5,0,0), (2,0,0) ]
    sweep_radii = [ 1., 1., 1., 1., 1. ]
    
    mesh = mesh_from_sweep_circles( sweep_path, sweep_radii )
    #debugger()
    mesh.write_OBJ( "sweep_surface-test.obj" )

def main():
    test()

if __name__ == "__main__": main()
