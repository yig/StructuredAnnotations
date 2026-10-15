from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.

def membrane_smooth_curve_range_preserving_tangents( curve, smooth_start_vertex_index, smooth_end_vertex_index ):
    '''
    Performs in-place smoothing of curve, such that the range
    [ smooth_start_vertex_index, smooth_end_vertex_index ]
    is replaced with the membrane solution.
    
    Note that
    curve[ smooth_start_vertex_index ] and
    curve[ smooth_start_vertex_index + 1 ] and
    curve[ smooth_end_vertex_index ] and
    curve[ smooth_end_vertex_index + 1 ]
    will all be constrained, so that tangents are preserved.
    
    If there would be no live vertices, an error message is printed out and the curve is returned
    unchanged.
    '''
    
    assert int( smooth_start_vertex_index ) == smooth_start_vertex_index
    assert int( smooth_end_vertex_index ) == smooth_end_vertex_index
    assert smooth_end_vertex_index < len( curve )
    assert smooth_start_vertex_index >= 0
    
    if smooth_start_vertex_index + 1 < smooth_end_vertex_index - 1:
        from sweep_spine import diffuse_values_along_sweep_spine
        domain = array( curve[ smooth_start_vertex_index : smooth_end_vertex_index + 1 ] )
        domain_n = smooth_end_vertex_index + 1 - smooth_start_vertex_index
        constraints = [
            ( 0, domain[ 0 ] ),
            ( 1, domain[ 1 ] ),
            ( domain_n-2, domain[ -2 ] ),
            ( domain_n-1, domain[ -1 ] )
            ]
        #print "smoothing domain size:", domain_n
        #print "domain (before):"
        #print domain
        #print "constraints:"
        #print constraints
        curve[ smooth_start_vertex_index : smooth_end_vertex_index + 1 ] = diffuse_values_along_sweep_spine(
            domain,
            constraints
            )
        #print "domain (after):"
        #print curve[ smooth_start_vertex_index : smooth_end_vertex_index + 1 ]
        #print "domain | before - after |^2:", ( ( curve[ smooth_start_vertex_index : smooth_end_vertex_index + 1 ] - domain )**2 ).sum()
    else:
        print("Can't smooth; smoothing domain is too small.")

def curve_from_base_curve_and_oversketch_curves():
    pass

def curve_from_base_curve_and_oversketch_curve( curve, oversketch_curve, scale = 1., ts = [] ):
    '''
    Given a curve 'curve' (sequence of 2d points) and
    a second, oversketching curve 'oversketch_curve',
    an optional scale parameter (default is 1.; the curves are assumed to be displayed in a region [-1/scale..1/scale, -1/scale..1/scale]),
    an optional list of t's along the curve (each t is an arc length parameter along 'curve'),
    returns a new curve by modifying 'curve' to follow 'oversketch_curve' where appropriate
    ('oversketch_curve' can modify a range of vertices, or it can extend 'curve').
    The optional input/output parameter 'ts' will be updated so that each of its t values is
    correct for the new curve.
    '''
    
    curve = asarray( curve )
    oversketch_curve = asarray( oversketch_curve )
    
    assert curve.shape[0] >= 2
    assert oversketch_curve.shape[0] >= 2
    
    ## 3 cases:
    ##   1 oversketch_curve replaces a range of curve
    ##   2 oversketch_curve replaces a range of curve, extending curve beyond its start or end
    ##   3 oversketch_curve replaces curve
    
    ## How close the oversketch curve should be to the curve to attach them.
    kDistanceThreshold = .05 / scale
    ## Whether or not to smooth.
    kSmooth = True
    ## How large of an area of the curve and oversketch curve to smooth.
    kSmoothCurvePercent = 20.
    
    from annotations_sweep import edge_index2t
    from vecutil import closest_distsqr_and_point_and_edge_index_on_line_strip_to_point
    
    
    ## case 1
    ## Does oversketch_curve start and end near curve?
    ## If so, replace that range with oversketch_curve, and smooth heavily the vertices near the
    ## start and end of the range (or, better, set the vertices to a membrane solution (laplacian)).
    ## t values before and after the oversketch_curve range should be preserved, and those
    ## within should be rescaled.
    start_info = closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
        curve, oversketch_curve[0]
        )
    end_info = closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
        curve, oversketch_curve[-1]
        )
    if start_info[0] < kDistanceThreshold**2 and end_info[0] < kDistanceThreshold**2:
        ## The oversketch curve shouldn't be a loop starting and ending at the same curve edge.
        if start_info[2] == end_info[2]: return curve
        
        ## Reverse oversketch_curve if needed.
        if start_info[2] > end_info[2]:
            start_info, end_info = end_info, start_info
            oversketch_curve = oversketch_curve[ ::-1 ]
        
        ## Build the new curve.
        result_curve = list( curve[ : start_info[2] + 1 ] )
        result_curve = result_curve + list( oversketch_curve )
        result_curve = result_curve + list( curve[ end_info[2] + 1 : ] )
        result_curve = asarray( result_curve )
        
        
        ## Smooth near the attachment points.
        if kSmooth:
            five_percent_of_curve = int(round( (kSmoothCurvePercent/100.) * len( curve ) ))
            five_percent_of_oversketch_curve = int(round( (kSmoothCurvePercent/100.) * len( oversketch_curve ) ))
            ## Clamp to 50% of the oversketch curve, so our smoothing of the start and end
            ## attachments don't conflict (smooth the same area with different boundary constraints).
            result_oversketch_middle_edge_index = start_info[2]+1 + len( oversketch_curve )//2
            
            ## The first attachment point.
            ## -1 because we want to constrain the tangent.
            smooth_start_vertex_index = start_info[2] - 1 - five_percent_of_curve
            ## +1 because we want to constrain the tangent.
            smooth_end_vertex_index = start_info[2]+1 + 1 + five_percent_of_oversketch_curve
            
            smooth_start_vertex_index = max( 0, smooth_start_vertex_index )
            smooth_end_vertex_index = min( result_oversketch_middle_edge_index+1, smooth_end_vertex_index )
            assert smooth_start_vertex_index >= 0
            assert smooth_end_vertex_index < len( result_curve )
            
            membrane_smooth_curve_range_preserving_tangents( result_curve, smooth_start_vertex_index, smooth_end_vertex_index )
            
            
            ## The second attachment point.
            ## -1 because we want to constrain the tangent.
            smooth_start_vertex_index = start_info[2]+1 + len( oversketch_curve )-1 - 1 - five_percent_of_oversketch_curve
            ## +1 because we want to constrain the tangent.
            smooth_end_vertex_index = start_info[2]+1 + len( oversketch_curve ) + 1 + five_percent_of_curve
            
            smooth_start_vertex_index = max( result_oversketch_middle_edge_index, smooth_start_vertex_index )
            smooth_end_vertex_index = min( len( result_curve )-1, smooth_end_vertex_index )
            assert smooth_start_vertex_index >= 0
            assert smooth_end_vertex_index < len( result_curve )
            
            membrane_smooth_curve_range_preserving_tangents( result_curve, smooth_start_vertex_index, smooth_end_vertex_index )
        
        
        ## Re-map ts along the curve.
        ## The edge attaching oversketch curve to 'curve' nearest the start.
        remap_old_t_start = edge_index2t( start_info[2], curve )
        remap_new_t_start = edge_index2t( start_info[2], result_curve )
        ## The edge attaching oversketch curve to 'curve' nearest the end.
        remap_old_t_end = edge_index2t( end_info[2], curve )
        remap_new_t_end = edge_index2t( start_info[2] + 1 + len( oversketch_curve ) - 1, result_curve )
        
        remap_old_t_span = remap_old_t_end - remap_old_t_start
        remap_new_t_span = remap_new_t_end - remap_new_t_start
        assert remap_old_t_span > 0
        assert remap_new_t_span > 0
        
        ## We have to watch out for remap_old_t_start == 0 and remap_old_t_end == 1, since
        ## we may divide by remap_old_t_start or ( 1 - remap_old_t_end ).
        ## Make sure to get our inequality signs right.
        ## UPDATE: This still breaks sometimes.  I'll catch ZeroDivisionError's and do
        ##         the right thing.
        new_ts = []
        for t in ts:
            assert t >= 0.
            assert t <= 1.
            if t < remap_old_t_start:
                try:
                    new_ts.append(
                        ( t / remap_old_t_start ) * remap_new_t_start
                        )
                except ZeroDivisionError:
                    new_ts.append( 0. )
            elif remap_old_t_start <= t and t <= remap_old_t_end:
                try:
                    new_ts.append(
                        ( ( t - remap_old_t_start ) / remap_old_t_span ) * remap_new_t_span + remap_new_t_start
                        )
                except ZeroDivisionError:
                    new_ts.append( remap_new_t_start )
            else:
                try:
                    new_ts.append(
                        ( t - remap_old_t_end ) / ( 1. - remap_old_t_end ) * ( 1. - remap_new_t_end ) + remap_new_t_end
                        )
                except ZeroDivisionError:
                    new_ts.append( remap_new_t_end )
        
        ## Clamp ts.
        new_ts = [ max( 0., min( 1., t ) ) for t in new_ts ]
        
        assert len( ts ) == len( new_ts )
        ts[:] = new_ts
        return result_curve
    
    
    ## case 2
    ## Does oversketch_curve start xor end near curve?
    ## If so, replace the range from the intersection (or near intersection) up to the end with
    ## oversketch_curve.  Just as in case 2, smooth or replace with membrane solution the area
    ## around the intersection.
    close_infos = [ info for info in ( start_info, end_info ) if info[0] < kDistanceThreshold**2 ]
    if len( close_infos ) > 0:
        assert len( close_infos ) == 1
        info = close_infos[0]
        
        ## Reverse oversketch_curve if needed, so that its start point is closest to 'curve'.
        if info is end_info:
            oversketch_curve = oversketch_curve[ ::-1 ]
        ## Test if the tangents point in the same direction to tell which piece of 'curve' to keep
        ## (before or after the closest edge index).
        if dot( curve[ info[2] + 1 ] - curve[ info[2] ], oversketch_curve[1] - oversketch_curve[0] ) > 0.:
            result_curve = list( curve[ : info[2] + 1 ] ) + list( oversketch_curve )
            result_curve = asarray( result_curve )
            
            ## The value that controls the remapping is the t of the edge that joins the two curves.
            remap_t_old = edge_index2t( info[2], curve )
            remap_t_new = edge_index2t( info[2], result_curve )
            
            ## Smooth near the attachment points.
            if kSmooth:
                five_percent_of_curve = (kSmoothCurvePercent/100.) * len( curve )
                five_percent_of_oversketch_curve = (kSmoothCurvePercent/100.) * len( oversketch_curve )
                ## -1 because we want to constrain the tangent.
                smooth_start_vertex_index = int(round( info[2] - 1 - five_percent_of_curve ))
                ## +1 because we want to constrain the tangent.
                smooth_end_vertex_index = int(round( info[2] + 1 + 1 + five_percent_of_oversketch_curve ))
                
                ## Clamp to 50% of the oversketch curve.
                ## UPDATE: Actually, I see no reason to priviledge the oversketch curve; just obey kSmoothCurvePercent.
                smooth_start_vertex_index = max( 0, smooth_start_vertex_index )
                #smooth_end_vertex_index = min( info[2] + 1 + len( oversketch_curve )//2 + 1, smooth_end_vertex_index )
                smooth_end_vertex_index = min( len(result_curve) - 1, smooth_end_vertex_index )
                assert smooth_start_vertex_index >= 0
                assert smooth_end_vertex_index < len( result_curve )
                
                membrane_smooth_curve_range_preserving_tangents( result_curve, smooth_start_vertex_index, smooth_end_vertex_index )
        
        else:
            result_curve = list( oversketch_curve )[::-1] + list( curve[ info[2] + 1 : ] )
            result_curve = asarray( result_curve )
            
            ## The value that controls the remapping is the t of the edge that joins the two curves.
            remap_t_old = edge_index2t( info[2], curve )
            remap_t_new = edge_index2t( len( oversketch_curve ) - 1, result_curve )
            
            ## Smooth near the attachment points.
            if kSmooth:
                five_percent_of_curve = int(round( (kSmoothCurvePercent/100.) * len( curve ) ))
                five_percent_of_oversketch_curve = int(round( (kSmoothCurvePercent/100.) * len( oversketch_curve ) ))
                ## -1 because we want to constrain the tangent.
                smooth_start_vertex_index = len( oversketch_curve )-1 - 1 - five_percent_of_oversketch_curve
                ## +1 because we want to constrain the tangent.
                smooth_end_vertex_index = len( oversketch_curve ) + 1 + five_percent_of_curve
                
                ## Clamp to 50% of the oversketch curve.
                ## UPDATE: Actually, I see no reason to priviledge the oversketch curve; just obey kSmoothCurvePercent.
                #smooth_start_vertex_index = max( len( oversketch_curve )//2 - 1, smooth_start_vertex_index )
                smooth_start_vertex_index = max( 0, smooth_start_vertex_index )
                smooth_end_vertex_index = min( len( result_curve ) - 1, smooth_end_vertex_index )
                assert smooth_start_vertex_index >= 0
                assert smooth_end_vertex_index < len( result_curve )
                
                membrane_smooth_curve_range_preserving_tangents( result_curve, smooth_start_vertex_index, smooth_end_vertex_index )
        
        ## Re-map ts along the curve.
        new_ts = []
        for t in ts:
            if t < remap_t_old:
                try:
                    new_ts.append(
                        ( t / remap_t_old ) * remap_t_new
                        )
                except ZeroDivisionError:
                    new_ts.append( 0. )
            else:
                try:
                    new_ts.append(
                        ( ( t - remap_t_old ) / ( 1. - remap_t_old ) ) * ( 1. - remap_t_new ) + remap_t_new
                        )
                except ZeroDivisionError:
                    new_ts.append( remap_t_new )
        
        ## Clamp ts.
        new_ts = [ max( 0., min( 1., t ) ) for t in new_ts ]
        
        assert len( ts ) == len( new_ts )
        ts[:] = new_ts
        
        return result_curve
    
    
    ## case 3
    ## TODO: If curve's endpoints are close enough to oversketch_curve's endpoints, average the two curves.
    #if mag2( curve[0] - oversketch_curve[0] ) < kDistanceThreshold**2 and mag2( curve[-1] - oversketch_curve[-1] ) < kDistanceThreshold**2:
    return oversketch_curve


def curve_from_bundle_of_curves():
    pass
