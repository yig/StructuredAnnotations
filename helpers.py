from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *
from itertools import chain

def PAINT():
    print('PAINT: executing previously unexecuted code path:')
    import traceback
    traceback.print_stack()
    
    #from pydb import debugger
    #debugger()

class Struct( object ): pass

def members( obj ):
    from builtins import dir as dirr
    return [ (el, getattr( obj, el )) for el in dirr(obj) if not el.startswith( '__' ) ]

def expect( value ):
    if not value:
        print('Expect failed!')
        import traceback
        traceback.print_stack()
        
        #import pydb
        #pydb.debugger()
    
    #if not value:
    #    import sys
    #    print >> sys.stderr, 'Expect failed!'

from geometry import line_line_intersection, line_segments_cross

def get_contours_bbox( contours ):
    '''
    Returns the bounding box of the list of line loops 'contours' as a 4-tuple
    ( min x, max x, min y, max y ).
    Each line loop is a list of 2d points.
    '''
    
    xs = [ pt[0] for contour in contours for pt in contour ]
    ys = [ pt[1] for contour in contours for pt in contour ]
    maxx = max( xs )
    minx = min( xs )
    maxy = max( ys )
    miny = min( ys )
    
    dim = ( minx, maxx, miny, maxy )
    return dim

def contour_contains( contour, point ):
    '''
    Takes a simple, closed 'contour' (list of 2d points) and a 2d 'point' and
    returns True if the 'point' is contained within the contour and False otherwise.
    '''
    
    ## If we have less than three points, we can't enclose any area, so the point isn't inside.
    if len( contour ) < 3: return False
    
    ## assert 2d points
    assert False not in [ len( cpt ) == 2 for cpt in list(contour) + [point] ]
    
    from math import sin, cos, sqrt, pi
    
    bbox = get_contours_bbox( [contour] )
    diameter = 2 * sqrt( max( abs( bbox[1] ), abs( bbox[0] ) )**2 + max( abs( bbox[3] ), abs( bbox[2] ) )**2 )
    
    ## We'll rotate the other endpoint about a a circle larger than the contour until
    ## there are no epsilon issues with the line segment intersections
    crossing_eps = 1e-5
    angle_incr = 1e-2
    angle = 0.
    
    while True:
        end_pt = [ diameter * cos( angle ), diameter * sin( angle ) ]
        
        crossings = 0
        eps_issues = False
        for i in range(len( contour )):
            edge = ( contour[i], contour[ (i+1) % len(contour) ] )
            
            t_p, t_e = line_line_intersection( ( point, end_pt ), edge )
            
            ## There is a crossing if both t_p and t_e are in [0,1]
            if t_p >= 0. and t_p <= 1. and t_e >= 0. and t_e <= 1.:
                crossings += 1
                
                ## If they intersect, the far end point should be way beyond any edge in the contour
                ## (no epsilon issues).
                assert t_p < 1. - crossing_eps
                
                ## There is an epsilon issue if t_e is in [0,eps) or (1-eps,1].
                ## NOTE: We can't check t_p because 'point' is fixed.
                if abs( t_e - 0. ) < crossing_eps or abs( t_e - 1. ) < crossing_eps:
                    print('Epsilon issue!')
                    eps_issues = True
        
        ## If there were no epsilon issues, return.
        ## If we rotated all the way around, return.
        if not eps_issues or angle >= 2*pi:
            ## An odd number of crossings means the point is inside the contour.
            return (crossings % 2) == 1
        
        angle += angle_incr

def ray_linestrip_intersection_ts( ray, linestrip, ray_is_line = False ):
    '''
    Given a 2d 'ray' ( 2d point, 2d direction vector ), returns the list of all intersections
    of the ray with the line strip 'linestrip' (sequence of 2d points) as an unsorted list of
    tuples of the form (
        t along 'ray',
        index i of a point in 'linestrip',
        t_edge along the edge between linestrip[i] and linestrip[i+1] such that lerp( linestrip[i], linestrip[i+1], t_edge ) equals ray[0] + t * ray[1] )
        ).
    If optional parameter 'ray_is_line' is True, then returned 't' values may be negative.
    
    untested
    '''
    
    assert len( ray[0] ) == 2
    assert len( ray[1] ) == 2
    ray_edge = ( ray[0], asarray( ray[0] ) + asarray( ray[1] ) )
    
    isect_t_eps = 1e-10
    
    ts = []
    
    for i in range( len( linestrip )-1 ):
        try:
            t_ray, t_edge = line_line_intersection( ray_edge, ( linestrip[i], linestrip[i+1] ) )
            if t_edge > -isect_t_eps and t_edge < 1 + isect_t_eps:
                ts.append( ( t_ray, i, t_edge ) )
        
        ## Skip edges coincident with the ray (not parallel, but exactly the same line).
        except RuntimeError: pass
    
    if not ray_is_line: ts = [ t_etc for t_etc in ts if t_etc[0] >= 0. ]
    return ts

def ray_lineloop_intersection_ts( ray, lineloop, ray_is_line = False ):
    '''
    The same as ray_linestrip_intersection_ts(), except the ray is intersected with the line loop
    (closed curve) 'lineloop'.
    '''
    return ray_linestrip_intersection_ts( ray, list(lineloop) + [lineloop[0]], ray_is_line )

def angle2D( one, two ):
    '''
    Returns the angle in radians [-pi,pi] turned clockwise between vector one and vector
    two in the sequence { 0, one, two }.
    In other words, if one were to turn clockwise by the returned angle,
    it would have the same direction as two:
        two == rotate( one, -angle2D( one, two ) )
    '''
    from myarray import asarrayf, dir
    from numpy import dot
    from math import acos
    
    one = dir( asarrayf( one ) )
    two = dir( asarrayf( two ) )
    
    angle = acos( max( -1., min( 1., dot( one,two ) ) ) )
    if cross_z( one,two ) < 0.:
        angle = -angle
    return angle

def cross_z( a,b ):
    '''
    Returns the z component of the cross product of vectors a,b.
    
    NOTE: tested
    '''
    return a[0]*b[1] - a[1]*b[0]

def solve_for_fixed_points( points, src0_src1, dst0_dst1 ):
    '''
    Given:
        src0_src1: a pair of 2d points (src0, src1)
        dst0_dst1: a pair of 2d points (dst0, dst1)
    
    Finds the 2D transform (uniform scale, rotate, translate)
    that maps src0 to dst0 and src1 to dst1.
    Returns 'points', a list of 2d points, with this transform applied.
    
    If any parameters have more than two coordinates, only the first two coordinates are used,
    and 2d coordinates are returned.
    '''
    
    ## Unpack the parameters
    (src0, src1) = src0_src1
    (dst0, dst1) = dst0_dst1
    
    from numpy import array
    from numpy.linalg import solve
    
    coefs = array((
        ( src0[0], -src0[1], 1, 0 ),
        ( src1 [0], -src1 [1], 1, 0 ),
        ( src0[1],  src0[0], 0, 1 ),
        ( src1 [1],  src1 [0], 0, 1 )
        ))
    rhs = array( (
        dst0[0],
        dst1[0],
        dst0[1],
        dst1[1]
        ) )
    
    sol = solve( coefs, rhs )
    
    xformed = [
        (
        sol[2] + sol[0] * pt[0] - sol[1] * pt[1],
        sol[3] + sol[1] * pt[0] + sol[0] * pt[1]
        )
        for pt in points
        ]
    
    return xformed

def x_intercept( x0_y0, x1_y1 ):
    '''
    Given:
        x0_y0: a pair of points (x0,y0)
        x1_y1: a pair of points (x1,y1)
    Returns x such y(x) = 0 for the line defined by
    the pair of points (x0,y0) and (x1,y1).
    '''
    
    ## Unpack the parameters.
    (x0,y0) = x0_y0
    (x1,y1) = x1_y1
    
    if y1 == y0: return x0 + 1e-5
    
    ## The two points of the line are (x0, y0) and (x1, y1).
    ##  (x0, y0) + t * ( x1 - x0, y1 - y0 ) = ( x_s, 0 )
    ## <=>
    ##  y0 + t * (y1 - y0) = 0
    ## <=>
    ##  t = -y0 / ( y1 - y0 )
    ## =>
    ##  x_s = x0 + ( -y0 / (y1 - y0) ) * (x1 - x0)
    return x0 + ( -y0 / (y1 - y0) ) * (x1 - x0)


def squared_distance_point_to_points( pt, points ):
    '''
    Returns the squared distance from each point in points to 'pt'.
    '''
    assert len( pt ) == 2
    
    pt = asarrayf( pt )
    points = asarrayf( points )
    
    result = ( ( points - pt )**2 ).sum(-1)
    return result

def index_of_closest_in_points_to_pt( points, pt ):
    '''
    Takes a list of 2d points 'points' and a 2d point 'pt'
    and returns the index of the closest point in 'points' to 'pt' and
    the distance to it as a tuple (closest pt index, distance).
    
    untested
    '''
    
    ## get the distances all at once (vectorized)
    dists = squared_distance_point_to_points( pt, points )
    assert len( dists ) == len( points )
    dists = list(zip( dists, list(range(len( points ))) ))
    #dists = [ ( mag2( cand - pt ), i ) for i, cand in enumerate( points ) ]
    ## take the minimum element
    result = min( dists )
    ## return it sqrt()'d
    return result[1], sqrt( result[0] )

def get_path_t_for_index( pts, index ):
    '''
    Returns 't' such that get_point_and_tangent_on_path( t, pts )
    returns the point in 'pts' at 'index'.
    
    NOTE: untested
    '''
    
    pts = asarrayf( pts )
    
    ## Make an array with the length of the entire line up to
    ## pts[i].
    lengths = [0.]
    for i in range( 1, len(pts) ):
        lengths.append( lengths[-1] + mag( pts[i] - pts[i-1] ) )
    
    percents = [max( 0., min( 1., x / lengths[-1] ) ) for x in lengths]
    return percents[ index ]

def get_point_on_path( t, pts, info = None ):
    '''
    Returns point (x,y) 't' along the given pts, as measured
    by path lengths.
    
    To cache information about 'pts' for faster subsequent calls (with different 't's),
    use the input/output parameter 'info.'  Pass an empty list on the first call, which
    will be modified to contain the cached data.  Pass the same parameter for subsequent calls
    to take advantage off the cached data.
    '''
    
    assert len( pts ) >= 2
    
    if None == info: info = []
    if len(info) == 0:
        info.append( Struct() )
        info[0].lengths = None
        info[0].percents = None
    
    from numpy import asarray, array, dot
    from numpy import float64 as Float
    pts = asarray( pts, Float )
    
    def get_lengths( pts ):
        ## Make an array with the length of the entire line up to
        ## pts[i].
        lengths = append( [0.], cumsum( sqrt( ( ( pts[1:] - pts[:-1] )**2 ).sum(1) ) ) )
        return lengths
    
    def get_percents( lengths ):
        percents = maximum( 0., minimum( 1., lengths / lengths[-1] ) )
        return percents
    
    if info[0].lengths is None: info[0].lengths = get_lengths( pts )
    lengths = info[0].lengths
    assert len( lengths ) == len( pts )
    
    if info[0].percents is None: info[0].percents = get_percents( lengths )
    percents = info[0].percents
    assert len( percents ) == len( pts )
    
    #print 't:',t,'percents:', percents
    
    ## clip t to be inside the path
    t = max( percents[0], min( percents[-1], t ) )
    
    for i in range( 1, len(percents) ):
        ## if 't' is in the segment preceding i
        if t <= percents[i]:
            
            try: tval = ( t - percents[i-1] ) / ( percents[i] - percents[i-1] )
            except ZeroDivisionError: tval = percents[i-1]
            pt_on_path = lerp( pts[i-1], pts[i], tval )
            return pt_on_path
            
        #~if t <= percent
    #~for i in percents
    
    raise RuntimeError("We shouldn't be here.")

def lerp( a,b,t ):
    '''
    Linear interpolation from a (when t == 0) to b (when t == 1).
    '''
    return a + (b-a) * t

def lerpVal( a,b,val ):
    '''
    Returns 't' such that lerp( a,b,t ) == val
    '''
    return (val - a) / (b - a)

def remove_contour_self_intersections( contour, extra = None, too_short_edge_length = 0. ):
    '''
    Prevents self-intersections in a contour (closed path of points).
    That is, the contour should be modified to have self-intersections removed
    (the shorter loop will be removed).
    Takes as input a contour 'contour' (closed path of points) and an optional parameter 'extra'
    which, if not None must have the same number of elements as 'contour'.
    Returns 'contour' with self-intersections removed, or (contour, extra) with elements
    corresponding to the sharp corners removed if 'extra' is not None.
    '''
    assert extra is None or len( contour ) == len( extra )
    
    contour = list( contour )
    if extra is not None: extra = list( extra )
    
    ## To do this we find two edges of the contour which intersect, and consider the two possible
    ## contours created by cutting the contour at the intersection point and choosing the longer
    ## of the two simpler (one less intersection) contours.  We repeat until we have a simple
    ## contour (no self intersections).
    ## Q: Should I remove very short edges?
    ## A: Yes, zero length edges are a problem.
    ## TODO Q: With multiple self-intersections, is there a proper order to consider them in?
    ##         For example, I think with a lower-case 'a' where instead of a sharp up-down there is
    ##         a very narrow loop that causes several self-intersections and only one
    ##         order for removing them results in the outline we want.
    ## TODO: I could make this work with open paths (line strip instead of line loop),
    ##       if I iterate e0/f0 until len( contour )-1.
    def cut_one_self_intersection( contour, extra ):
        assert extra is None or len( contour ) == len( extra )
        
        for e0 in range( len( contour ) ):
            e1 = (e0 + 1) % len( contour )
            ## Edge 'e' is from contour[e0] to contour[e1].
            
            for f0 in range( e0 + 2, len( contour ) ):
                ## Don't compare adjacent edges or the same edge against itself.
                ## Actually, the way my loop counters work this should never be true
                ## unless e0 is 0 and f0 is the last (looping) edge.
                assert ( e0 == 0 and f0 == len( contour )-1 ) or not (abs( e0 - f0 ) <= 1 or abs( e0 - f0 ) == len( contour )-1)
                if abs( e0 - f0 ) <= 1 or abs( e0 - f0 ) == len( contour )-1: continue
                
                f1 = (f0 + 1) % len( contour )
                ## Edge 'f' is from contour[f0] to contour[f1].
                
                try:
                    t_e, t_f = line_line_intersection(
                        (contour[e0], contour[e1]),
                        (contour[f0], contour[f1])
                        )
                except RuntimeError:
                    ## Skip coincident lines.
                    continue
                ## If the line segments intersect
                #isect_t_eps = kEps
                isect_t_eps = 1e-10
                if t_e > -isect_t_eps and t_f > -isect_t_eps and t_e < 1 + isect_t_eps and t_f < 1 + isect_t_eps:
                    isect = tuple( lerp( asarray( contour[e0] ), asarray( contour[e1] ), t_e ) )
                    
                    ## This should be true because of the way we loop through edges.
                    assert f0 > e0
                    
                    ## Construct the inner loop of the contour.
                    contour_inner = contour[:e0+1] + [isect] + contour[f0+1:]
                    if extra is not None: extra_inner = extra[:e0+1] + [[]] + extra[f0+1:]
                    else: extra_inner = None
                    
                    ## Construct the outer loop of the contour.
                    contour_outer = [isect] + contour[e0+1:f0+1]
                    if extra is not None: extra_outer = [[]] + extra[e0+1:f0+1]
                    else: extra_outer = None
                    
                    ## How to choose: a) the longer piece of contour
                    ##                or
                    ##                b) the longer piece of generating stroke
                    ##                or
                    ##                c) count the number of points in the two contours
                    ##                (but this assumes something about them, such as the generating
                    ##                points we traced from were evenly spaced --- which they weren't.)
                    ## NOTE: I think this only makes a difference for pretty broken skeletons,
                    ##       such as if I use contour length then result15/Palatino-Italic 'd' is not what I want.
                    
                    ## a)
                    if line_loop_length( contour_inner ) > line_loop_length( contour_outer ):
                    ## c)
                    #if len( contour_inner ) > len( contour_outer ):
                    ## b)
                    ## NOTE: this assumes extra is the trace points (or at least that extra is a
                    ##       list of lists with the first element the 2d skeleton point).
                    #if line_loop_length( [tr[0] for tr in extra_inner if len(tr) > 0] ) > line_loop_length( [tr[0] for tr in extra_outer if len(tr) > 0] ):
                        return contour_inner, extra_inner, True
                    else:
                        return contour_outer, extra_outer, True
        
        return contour, extra, False
    
    if extra is None: contour = remove_too_short_edges( too_short_edge_length, contour, extra )
    else: contour, extra = remove_too_short_edges( too_short_edge_length, contour, extra )
    while True:
        ## Cut until we have nothing left to do.
        contour, extra, did_cut = cut_one_self_intersection( contour, extra )
        if not did_cut: break
        print('cut self-intersection')
        ## also remove zero-length lines
        if extra is None: contour = remove_too_short_edges( too_short_edge_length, contour, extra )
        else: contour, extra = remove_too_short_edges( too_short_edge_length, contour, extra )
    
    if None == extra: return contour
    else: return contour, extra

def remove_too_short_edges( too_short_edge, contour, extra = None ):
    '''
    Removes edges in the contour if they are too short.
    Takes as input a distance threshold for determining too short edges
    'too_short_edge', a contour 'contour' (closed path of points),
    and an optional parameter 'extra' which, if not None,
    must have the same number of elements as 'contour'.
    Returns 'contour' with short edges removed, or (contour, extra) with elements corresponding
    to the short edges removed if 'extra' is not None.
    '''
    assert extra is None or len( contour ) == len( extra )
    
    contour = list( contour )
    if extra is not None: extra = list( extra )
    
    while len( contour ) > 2:
        
        cut_once = False
        for i in range(len( contour )):
            i_next = (i + 1) % len(contour)
            v_next = asarray( contour[i_next] ) - asarray( contour[i] )
            
            if mag2( v_next ) <= too_short_edge**2:
                print('edge length %f < too_short_edge %f' % (mag( v_next ), too_short_edge))
                
                del contour[ i ]
                if extra is not None: del extra[ i ]
                
                cut_once = True
                break
        
        ## If we looped through without removing any vertices, we're done.
        if not cut_once: break
    
    if extra is None: return contour
    else: return contour, extra

def resample_line_loop( curve, N = 100 ):
    '''
    Given a line loop (list of 2d points) 'curve', resample it evenly by arc length N times.
    '''
    
    ## We can't do anything for a line loop with less than 3 points.
    if len( curve ) < 3: return curve
    
    ## We need at least three samples to make a loop.
    N = max( 3, N )
    
    return resample_line_strip( list(curve) + [curve[0]], N+1 )[:-1]

def resample_line_strip( curve, N = 100 ):
    '''
    Given a line strip (list of 2d points) 'curve', resample it evenly by arc length N times.
    '''
    
    ## We can't do anything for a line strip with only a single point.
    if len( curve ) < 2: return curve
    
    ## We need at least two samples to make a line.
    N = max( 2, N )
    
    info = []
    return [ get_point_on_path( float(i)/(N-1), curve, info ) for i in range(N) ]

def smooth_open_stroke1( curve, smooth = 10, resample_N = 100 ):
    '''
    Given an open curve (simple, open list of (x,y) points), returns a smoother
    (and optionally evenly sampled) curve.
    '''
    
    
    ## 1 Resample the path
    ## 2 Smooth the path
    
    
    ## 1
    if resample_N is not None: curve = resample_line_strip( curve, resample_N )
    
    ## 2
    ## Smooth it by moving each point to the midpoint of neighboring points.
    curve = [ asarray( pt ) for pt in curve ]
    for i in range(smooth):
        
        curve2 = list(curve)
        for ci in range( 1, len( curve ) - 1 ):
            ## This jumps all the way to the "smoothest" position, the average of the neighbors.
            ## This is equivalent to a forward timestep of 1.
            ## Unfortunately, neighbors are also jumping, so it's sometimes (occasionally)
            ## not stable.  I could instead lerp only part of the way, but it's rarely a problem:
            #curve2[ci] = lerp( .5 * ( curve[ ci - 1 ] + curve[ ci + 1 ] ), curve[ ci ], .5 )
            curve2[ci] = .5 * ( curve[ ci - 1 ] + curve[ ci + 1 ] )
        
        curve = curve2
        
        ## Resample again
        ## I don't think this is needed
        #if resample_N is not None: curve = resample_line_loop( curve, resample_N )
    
    return curve

def smooth_open_stroke2( curve, smooth = 10, resample_N = 100 ):
    '''
    Given an open curve (simple, open list of (x,y) points), returns a smoother
    (and optionally evenly sampled) curve.
    '''
    
    
    ## 1 Resample the path (round up to power-of-two).
    ## 2 Smooth the path.
    
    
    ## 1
    if resample_N is None:
        resample_N = len( curve )
    
    '''
    ## Round up to a power of 2.
    ## UPDATE: I don't think this is necessary.
    n = 1
    while n < resample_N:
        n *= 2
    '''
    n = resample_N
    ## We need at least 6 points.
    resample_N = max( 6, n )
    curve = resample_line_strip( curve, resample_N )
    ## We also need an array-like curve.
    curve = asarray( curve )
    
    ## 2
    ## Smooth it by subsampling and then upsampling a la
    ## "Local Filters of B-spline Wavelets" (Samavati and Bartels 2004)
    
    def decompose( F ):
        n = len( F )
        ## I don't think this is necessary.
        #assert ( n % 2 ) == 0
        r = n//2 - 3
        assert r >= 0
        
        ## subsample
        C = zerosf( ( 4+r, F.shape[1] ) )
        C[0] = F[0]
        C[1] = -1./2.*F[0] + F[1] + 3./4.*F[2] - 1./4.*F[3]
        j = 2
        for i in range( 2, n-5, 2 ):
            C[j] = -1./4.*F[i] + 3./4.*F[i+1] + 3./4.*F[i+2] - 1./4.*F[i+3]
            j = j + 1
        C[j] = -1./4.*F[-4] + 3./4.*F[-3] + F[-2] - 1./2.*F[-1]
        C[j+1] = F[-1]
        assert j+1 == len(C)-1
        
        return C
    
    def recompose( C ):
        r = len( C ) - 4
        assert r >= 0
        n = 2 * ( r + 3 )
        ## I don't think this is necessary.
        #assert ( n % 2 ) == 0
        
        ## upsample
        F = zerosf( ( n, C.shape[1] ) )
        F[0] = C[0]
        F[1] = 1./2.*C[0] + 1./2.*C[1]
        j = 2
        for i in range( 1, 4+r - 2 ):
            F[j] = 3./4.*C[i] + 1./4.*C[i+1]
            F[j+1] = 1./4.*C[i] + 3./4.*C[i+1]
            j = j + 2
        F[j] = 1./2.*C[-2] + 1./2.*C[-1]
        F[j+1] = C[-1]
        assert j+1 == len(F)-1
        
        return F
    
    try:
        
        C = curve
        #print "smooth:", smooth
        #smooth = int(ceil(log2( smooth )))
        #print "log2 smooth:", smooth
        for i in range( smooth ):
            if len( C ) < 6:
                smooth = i
                break
            C = decompose( C )
        F = C
        for i in range( smooth ):
            F = recompose( F )
    
    except:
        try: from pydb import debugger
        except ImportError:
            import pdb
            def debugger(): pdb.set_trace()
        debugger()
    
    return F
## TODO: Use scipy.interpolate.fitpack.splprep / splev
def smooth_open_stroke3( curve, smooth = 10, resample_N = 100 ):
    '''
    Given an open curve (simple, open list of (x,y) points), returns a smoother
    (and optionally evenly sampled) curve.
    '''
    
    if resample_N is None: resample_N = len( curve )
    
    from scipy.interpolate.fitpack import splprep
    from scipy.interpolate.fitpack import splev
    
    curve = asarray( curve )
    
    m = len( curve )
    min_s, max_s = ( m-sqrt(2*m), m+sqrt(2*m) )
    
    #from pydb import debugger
    #debugger()
    
    try:
        spline, u = splprep( curve.T, s = lerp( min_s, max_s, smooth/10. ), quiet = True )
        smooth_curve = asarray( splev( linspace( 0., 1., resample_N ), spline ) ).T
        return smooth_curve
    except TypeError:
        return smooth_open_stroke1( curve, smooth, resample_N )

smooth_open_stroke = smooth_open_stroke1
#smooth_open_stroke = smooth_open_stroke2
#smooth_open_stroke = smooth_open_stroke3

def smooth_closed_stroke( curve, smooth = 1, resample_N = 100, simplify = True ):
    '''
    Given a curve (simple, closed list of (x,y) points), returns a smoother
    (and optionally evenly sampled) curve.
    '''
    
    
    ## 1 Clip the path if it self-intersects.  (Choose the longer of the two sides.)
    ## 2 Resample the path
    ## 3 Smooth the path
    
    
    ## 1
    if resample_N is not None: curve = resample_line_loop( curve, resample_N )
    if simplify: curve = remove_contour_self_intersections( curve )
    '''
    ## We start out indexing at 2, because the first two segments can't intersect each other.
    for ci in xrange( 2, len( curve ) - 1 ):
        ## If the segment curve[ ci : ci+2 ] intersects any previous segments
        ## (up to curve[ ci-2 : ci ]), clip the curve there.
        if True in helpers.line_segments_cross( 
    '''
    
    ## 2
    if resample_N is not None: curve = resample_line_loop( curve, resample_N )
    
    ## 3
    ## TODO: Use scipy.interpolate.fitpack.splprep / splev
    
    ## 3
    ## Don't allow 100% smoothing; it will return all zeros.
    return ft_smooth_closed_uniformly_sampled_curve( curve, min( smooth, 90 ) )
    
    ## 3
    ## Smooth it by moving each point to the midpoint of neighboring points.
    curve = [ asarray( pt ) for pt in curve ]
    for i in range(smooth):
        
        ## To make my indexing easier, put the first point at the end and the last point
        ## at the beginning.
        curve_wrap = [curve[-1]] + list(curve) + [curve[0]]
        curve2 = list( curve_wrap )
        
        for ci in range( 1, len( curve_wrap ) - 1 ):
            curve2[ci] = .5 * ( curve_wrap[ ci - 1 ] + curve_wrap[ ci + 1 ] )
        
        curve = curve2[1:-1]
        
        ## Resample again
        ## I don't think this is needed
        #if resample_N is not None: curve = resample_line_loop( curve, resample_N )
    
    return curve
smooth_curve = smooth_closed_stroke

def ft_smooth_closed_uniformly_sampled_curve( curve, percent_deviation = 10 ):
    '''
    Given a uniformly sampled, closed N-dimensional curve curve[0], curve[1], ..., curve[n-1],
    returns a smoother version of 'curve' that deviates from 'curve' by less than
    'percent_deviation' percent.
    
    used
    '''
    
    ## Assumes the curve is uniformly sampled and closed.
    
    curve = asarray( curve )
    for i in range( curve.shape[1] ):
        curve[:,i] = ft_smooth_periodic_function( curve[:,i], percent_deviation )
    
    return curve

def ft_smooth_periodic_function( x, percent_deviation = 10 ):
    '''
    Given a uniformly sampled periodic function x[0], x[1],, ..., x[n-1],
    returns a smoother version of 'x' that deviates from 'x' by less than
    'percent_deviation' percent.
    
    NOTE: This will return all zeros if percent_deviation is 100.
    
    used
    '''
    
    ## TODO Q: Should I change the inputs to have some definite notion of smoothness, so that
    ##         calling the function on its output with the same smoothing parameter has no effect.
    ## A1: Yes, to prevent over-smoothing.
    ## A2: No, I want the user to be able to choose a "smooth curve" operation that smoothes a
    ##     certain amount each time, without having to keep track of the smoothness parameter.
    
    percent_deviation = max( 0., min( 100., percent_deviation ) )
    
    X = fft.rfft( x )
    X_mags = sqrt( X.real**2 + X.imag**2 ).cumsum()
    
    target_mag = ((100 - percent_deviation)/100.) * X_mags[-1]
    assert target_mag >= 0.
    
    from bisect import bisect_left
    i = bisect_left( X_mags, target_mag )
    assert X_mags[i] >= target_mag
    ## TODO Q: Should we not touch the lowest frequency (i == 0) or perhaps shift it towards 1?
    assert i == 0 or X_mags[i-1] < target_mag
    
    X[i+1:] = 0.
    X[i] *= ( target_mag - ( X_mags[i-1] if i > 0 else 0. ) ) / X_mags[i]
    
    return fft.irfft( X, len(x) )

def line_strips_intersect( linestrip1, linestrip2 ):
    '''
    Returns True if the two line strips 'linestrip1' and 'linestrip2' intersect,
    where each is a sequence of 2d points.
    Returns False otherwise.
    
    untested
    '''
    
    for i in range( len( linestrip1 )-1 ):
        for j in range( len( linestrip2 )-1 ):
            if line_segments_cross( ( linestrip1[i], linestrip1[i+1] ), ( linestrip2[j], linestrip2[j+1] ) ):
                return True
    
    return False

def polyline_self_intersects( polyline, polyline_is_loop = False ):
    '''
    Returns True if the polyline 'polyline' intersects itself,
    where the polyline is a sequence of 2d points.
    Returns False otherwise.
    If optional parameter 'polyline_is_loop' is True, then 'polyline'
    is interpreted as a closed line loop with an implicit edge between
    the last point and the first point.
    The default value of 'polyline_is_loop' is False.
    
    untested
    '''
    
    end_buffer = 0 if polyline_is_loop else 1
    
    for e0 in range( len( polyline ) - end_buffer ):
        e1 = ( e0 + 1 ) % len( polyline )
        assert e0 != e1
        assert polyline_is_loop or e1 > e0
        
        ## The inner loops only tests edges whose first index is
        ## greater than 'e0'; edges whose first index is less than 'e0'
        ## will have been tested earlier.
        ## Two line segments that share a vertex will always self-intersect;
        ## start 2 away from 'e0'.
        for f0 in range( e0 + 2, len( polyline ) - end_buffer ):
            f1 = ( f0 + 1 ) % len( polyline )
            assert f0 != f1
            assert polyline_is_loop or f1 > f0
            
            if line_segments_cross(
                ( polyline[e0], polyline[e1] ),
                ( polyline[f0], polyline[f1] )
                ):
                return True
    
    return False

def line_loop_length( pts ):
    '''
    Returns the total length of the line loop 'pts' (list of 2d points).
    NOTE: 'pts' is interpreted as a closed path (there is an edge between pts[-1] and pts[0]).
    '''
    if len( pts ) == 0: return 0.
    else: return line_strip_length( list( pts ) + [pts[0]] )

def line_strip_length( pts ):
    '''
    Returns the total length of the line strip 'pts' (list of 2d points).
    NOTE: 'pts' is interpreted as an open path, not a closed path (there is no edge between pts[-1] and pts[0]).
    '''
    seg_lengths = line_strip_seg_lengths( pts )
    return sum( seg_lengths )

def line_strip_seg_lengths( pts ):
    '''
    Returns the length of each edge in the line strip 'pts' (list of 2d points).
    The returned list will have length len( pts ) - 1.
    NOTE: 'pts' is interpreted as an open path, not a closed path (there is no edge between pts[-1] and pts[0]).
    '''
    pts = asarray( pts )
    seg_lengths = sqrt( ( ( pts[1:] - pts[:-1] )**2 ).sum(1) )
    return seg_lengths

def closest_point_index( pt, list_of_pts ):
    '''
    Returns the index of the point in 'list_of_pts' closest to the point passed as 'pt'.
    '''
    pt = asarray( pt ).reshape( ( 1, len(pt) ) )
    list_of_pts = asarray( list_of_pts )
    
    assert prod( list_of_pts.shape ) >= prod( pt.shape )
    assert ( prod( list_of_pts.shape ) % prod( pt.shape ) ) == 0
    
    return ((list_of_pts - pt)**2).sum(1).argmin()

def farthest_point_index( pt, list_of_pts ):
    '''
    Returns the index of the point in 'list_of_pts' farthest from the point passed as 'pt'.
    '''
    pt = asarray( pt ).reshape( ( 1, len(pt) ) )
    list_of_pts = asarray( list_of_pts )
    
    assert prod( list_of_pts.shape ) >= prod( pt.shape )
    assert ( prod( list_of_pts.shape ) % prod( pt.shape ) ) == 0
    
    return ((list_of_pts - pt)**2).sum(1).argmax()

def tangents_and_normals_along_line_strip_2d( line_strip ):
    '''
    Given an open curve as a sequence of 2d points 'line_strip',
    returns a tuple containing two sequences,
        unit tangent vectors and
        unit (right-turn) normals,
    Each sequence contains one 2d vector for each edge in 'line_strip'.
    
    NOTE: If points in line_strip have dimension higher than two, tangents will still be correct,
          but normals will be a 90-degree clockwise rotation of the tangents about the z axis.
    
    used
    '''
    
    ## Any line strip must have at least two points
    assert len( line_strip ) >= 2
    
    ## Convert the line strip to an array, so we can vectorize our computations.
    line_strip = asarrayf( line_strip )
    ## Normalize each edge.
    unit_edges = normalize_vectors( line_strip[1:] - line_strip[:-1] )
    
    normals = array( unit_edges )
    ## This is a 90-degree right-turn (clockwise rotation) about the z axis.
    normals[:,0] = unit_edges[:,1]
    normals[:,1] = -unit_edges[:,0]
    
    return unit_edges, normals

def normalize_vectors( vecs, eps = 1e-8 ):
    '''
    Returns the array-like sequence of vectors, normalized.  If the vector has magnitude less
    than optional parameter 'eps', it will be returned unchanged.
    '''
    
    vecs = asarray( vecs )
    
    vec_norms = sqrt( ( vecs**2 ).sum(1) )
    vec_norms[ vec_norms < eps ] = 1.
    return vecs / vec_norms[:,newaxis]

def line_loop_from_tilted_circle( center, radius, normal, num_samples = 40 ):
    '''
    Given a tilted circle defined by a 3d center point, a radius, and a normal,
    returns a line loop sampling the circle.
    The number of samples is specified by the optional parameter 'num_samples'.
    '''
    
    center = asarray( center ).ravel()
    
    from vecutil import get_rotation_matrix, get_axis_angle_rotation
    rot = asarray( get_rotation_matrix( *get_axis_angle_rotation( ( 0,0,1 ), normal ) ) )
    
    '''
    return [
        center + dot( rot, array( ( radius*cos( theta ), radius*sin( theta ), 0. ) ) )
        for theta in linspace( 0., 2*pi, num_samples )[:-1]
        ]
    '''
    
    theta = linspace( 0., 2*pi, num_samples )[:-1]
    return center + radius * dot( rot, vstack( ( cos( theta ), sin( theta ), zeros( theta.shape ) ) ) ).T

def topo_less( values, i, j ):
    '''
    Topologically perfect less-than value comparison.
    If indices differ (i != j), then one value is always less than the other
    (values[i] < values[j] or values[j] < values[i]).
    Uses a tie-breaking scheme based on comparing the indices themselves, so
    equality only occurs when i == j.
    '''
    if values[i] < values[j]: return True
    if values[i] == values[j] and i < j: return True
    return False

def histogram( vals ):
    '''
    Given a sequence of values 'vals', return a dictionary mapping a value to the number
    of times it appears in 'vals'.
    '''
    
    hist = {}
    for val in vals:
        hist.setdefault( val, 0 )
        hist[ val ] += 1
    return hist

def isiterable( a ):
    '''
    Returns True is 'a' is iterable, False otherwise.
    '''
    ## http://bytes.com/forum/thread514838.html
    try:
        it = iter(a)
        return True
    except TypeError:
        return False

def color_spiral( i, N ):
    '''
    Takes parameters 'i' and 'N' where 'i' is a number from 0 to 'N'.
    Returns a 3-tuple ( red, green, blue ) of well-stratified colors, each value in the range [0,1].
    '''
    #from math import sin, sqrt
    #col = ( sin( 1 + i ), sin( 1 + 7 * i ), sin( 1 + 11 * i ) )
    #from random import random
    #col = ( random(), random(), random() )
    from colorsys import hsv_to_rgb
    t = float( i ) / N
    col = hsv_to_rgb(
        ( 2 * sqrt(t) * t ) % 1,
        sqrt( 1. - .8 * t ),
        1
        )
    
    col = tuple( map( abs, list( col ) ) )
    
    return col

def shortest_path( vs, v2neighbors, source, dest, **kwargs ):
    '''
    Given a list of vertex objects 'vs',
    a function 'v2neighbors' mapping a vertex object to a sequence of its neighboring vertex objects,
    source vertex object 'source', and
    destination vertex object 'dest',
    returns the list of vertex objects [ 'source', ..., 'dest' ] such that there is an edge
    between each successive pair of vertex objects and the path returned is the shortest path
    according to the (by default) number of edges between the vertex objects.
    If no such path exists, returns None.
    
    Optional named parameters:
        metric, a function that takes two adjacent vertex objects and returns the distance between them.
                defaults to a function always returning 1.
    
    tested
    '''
    
    ## Adapted from annotations_rotund.shortest_path_on_mesh()
    
    kwargs.setdefault( 'metric', lambda v1, v2: 1 )
    metric = kwargs[ 'metric' ]
    
    ## We follow the wikipedia article on Dijkstra's algorithm:
    ##   http://en.wikipedia.org/wiki/Dijkstra's_algorithm
    dist = {}
    previous = {}
    
    dist[ source ] = 0.
    Q = set( vs )
    while len( Q ) > 0:
        distQ = [ ( dist[v], v ) for v in Q if v in dist ]
        if len( distQ ) == 0: break
        u = min( distQ )[1]
        Q.remove( u )
        if u == dest: break
        for v in v2neighbors( u ):
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

def tEquidistant( p1, dir1, p2, dir2 ):
    '''
    Returns smallest positive 't' such that
        distance( p1 + t * dir1, p2 + t * dir2 ) == t
    If no such positive t exists, returns 0.
    If dir1 and dir2 are parallel, or if one of them is zero and
    the other is perpendicular to (p1-p2), returns distance( p1, p2 ).
    NOTE: If both dir1 and dir2 are 0, distance( p1, p2 ) is returned.
    '''
    
    eps = 1e-5
    
    #print 'p1:', p1, 'dir1:', dir1, 'p2:', p2, 'dir2:', dir2
    
    dp = asarray( p1 ) - asarray( p2 )
    
    #'''
    ## if either dir is 0, return t equal to the distance between p1 and p2.
    ## if dot(dir1,dir1) < kEps or dot(dir2,dir2) < kEps or abs(dot(dir1,dir2) - 1.) < kEps:
    ## XXX no, how about:
    ## if dir1 and dir2 are parallel or dirN is 0 and dirM is
    ## perpendicular to dp, return t equal to dp.
    if( abs(dot(dir1,dir2) - 1.) < eps
        or ( dot(dir1,dir1) < eps and dot(dir2,dp) < eps )
        or ( dot(dir2,dir2) < eps and dot(dir1,dp) < eps )
        ):
        return sqrt( dot(dp,dp) )
    #'''
    
    ## assuming dir1, dir2 are normalized
    dd = asarray( dir1 ) - asarray( dir2 )
    #print 'dd:', dd
    
    ## We want t such that
    ## t = | ( p1 + t * dir1 ) - ( p1 + t * dir2 ) |
    ##   = | ( p1 - p2 ) + t * ( dir1 - dir2 ) |
    ##   = | dp + t * dd |
    ## <=>
    ## t^2 = ( dp + t * dd ) ( dp + t * dd )
    ##     = dot( dp, dp ) + t^2 * dot( dd, dd ) + 2 * t * dot( dp, dd )
    ## <=>
    ## 0 = t^2 * ( 1 - dot( dd, dd ) + t * -2 * dot( dp, dd ) - dot( dp, dp )
    
    a =  1. - dot( dd, dd )
    b = -2. * dot( dp, dd )
    c = -1. * dot( dp, dp )
    
    #print 'a:',a,'b:',b,'c:',c
    ## if a is zero, we have a linear equation
    if abs(a) < eps:
        return -c / b
    
    ## quadratic equation
    t_discrim = sqrt( b*b - 4*a*c )
    t_numer1 = -b + t_discrim
    t_numer2 = -b - t_discrim
    t_denom = 2*a
    
    t1 = t_numer1 / t_denom
    t2 = t_numer2 / t_denom
    
    
    # return max( t1, t2 )
    
    ## we want the smallest positive t
    # return min( [t for t in (t1, t2) if t > 0.] )
    
    ts = [t for t in (t1, t2) if t > 0.]
    #print ts
    if len( ts ) > 0:
        return min( ts )
    else:
        return 0.

def point_is_inside_rect( pos, rect ):
    '''
    Given a 2d point 'pos' and a sequence of two 2d positions specifying the corners of a rectangle
    'rect',
    returns True if the point is inside the rectangle, and False otherwise.
    NOTE: 'rect' may be a sequence of four values; in other words, the sequence of two 2d values can
          be flattened.
    '''
    
    pos = asarray( pos )
    rect = asarray( rect ).reshape( (2,2) )
    return ( pos < rect.max( axis = 0 ) ).all() and ( pos > rect.min( axis = 0 ) ).all()

def asin_periodic( x ):
    '''
    A periodic version of asin; when |x| > 1, x is shifted back into range and the number of shifts
    determines the number of extra radians to add to the angle (positive or negative).
    The result angle is not limited to any bounds, and is roughly proportional to x.
    '''
    from math import asin, pi
    
    #return pi * max( -1., min( 1., x ) )
    
    x_shift = ( ( x + 1 ) % 2 ) - 1
    return asin( max( -1., min( 1., x_shift ) ) ) + pi * floor( ( x + 1 ) * .5 )

def lines_to_line_strips( lines ):
    '''
    Given a sequence of line segments (pairs of points), returns a list of line strips (sequences
    of connected points).
    '''
    
    eps = 1e-5
    
    Q = list( lines )
    
    strips = []
    while len( Q ) > 0:
        strip = list(map( asarray, Q.pop() ))
        ## Iterate last to first because we will delete elements as we go.
        for i in range(len(Q))[::-1]:
            line = Q[i]
            
            ## Attach to the end?
            if mag2( strip[-1] - line[0] ) < eps**2:
                strip.append( asarray( line[1] ) )
                del Q[i]
            elif mag2( strip[-1] - line[1] ) < eps**2:
                strip.append( asarray( line[0] ) )
                del Q[i]
            
            ## Attach to the beginning?
            elif mag2( strip[0] - line[0] ) < eps**2:
                strip.insert( 0, asarray( line[1] ) )
                del Q[i]
            elif mag2( strip[0] - line[1] ) < eps**2:
                strip.insert( 0, asarray( line[0] ) )
                del Q[i]
        
        strips.append( strip )
    
    '''
    ## To debug drawing glitches.
    from pprint import pprint
    print "strips:"
    pprint( strips )
    max_strips = [ abs( asarray( strip ) ).max( axis = 0 ) for strip in strips ]
    print "max_strips:"
    pprint( max_strips )
    if True in [ max_strip.max() > 1. for max_strip in max_strips ]:
        from pydb import debugger
        debugger()
    '''
    return strips

def test_shortest_path():
    import helpers
    vs = list(range(10))
    print("vs:", vs)
    ## [0,1,2,3]
    print(helpers.shortest_path( vs, lambda v: [ n for n in (v-1,v+1) if n in vs ], 0, 3 ))
    ## [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    print(helpers.shortest_path( vs, lambda v: [ n for n in (v-1,v+1) if n in vs ], 0, 9 ))
    ## [3, 4, 5, 6, 7, 8, 9]
    print(helpers.shortest_path( vs, lambda v: [ n for n in (v-1,v+1) if n in vs ], 3, 9 ))
    ## [9, 8, 7, 6, 5, 4, 3]
    print(helpers.shortest_path( vs, lambda v: [ n for n in (v-1,v+1) if n in vs ], 9, 3 ))
    ## [9, 0, 1, 2, 3]
    print(helpers.shortest_path( vs, lambda v: [ n for n in ((v+len(vs)-1) % len(vs), (v+1) % len(vs)) if n in vs ], 9, 3 ))
    ## [3, 2, 1, 0, 9]
    print(helpers.shortest_path( vs, lambda v: [ n for n in ((v+len(vs)-1) % len(vs), (v+1) % len(vs)) if n in vs ], 3, 9 ))
    ## [1, 0, 9]
    print(helpers.shortest_path( vs, lambda v: [ n for n in ((v+len(vs)-1) % len(vs), (v+1) % len(vs)) if n in vs ], 1, 9 ))
    ## [4, 3, 2, 1, 0, 9]
    print(helpers.shortest_path( vs, lambda v: [ n for n in ((v+len(vs)-1) % len(vs), (v+1) % len(vs)) if n in vs ], 4, 9 ))
    ## [5, 6, 7, 8, 9]
    print(helpers.shortest_path( vs, lambda v: [ n for n in ((v+len(vs)-1) % len(vs), (v+1) % len(vs)) if n in vs ], 5, 9 ))
    vs = list(range(11))
    print("vs:", vs)
    ## [5, 6, 7, 8, 9]
    print(helpers.shortest_path( vs, lambda v: [ n for n in ((v+len(vs)-1) % len(vs), (v+1) % len(vs)) if n in vs ], 5, 9 ))
    ## [4, 5, 6, 7, 8, 9]
    print(helpers.shortest_path( vs, lambda v: [ n for n in ((v+len(vs)-1) % len(vs), (v+1) % len(vs)) if n in vs ], 4, 9 ))
    ## [3, 2, 1, 0, 10, 9]
    print(helpers.shortest_path( vs, lambda v: [ n for n in ((v+len(vs)-1) % len(vs), (v+1) % len(vs)) if n in vs ], 3, 9 ))
