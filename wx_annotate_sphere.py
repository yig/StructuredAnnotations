#!/usr/bin/env python

from math import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *

try: from pydb import debugger
except ImportError:
    import pdb
    def debugger(): pdb.set_trace()



import wx
from wxBasicGLCanvas import BasicGLCanvas
from ui_actions import *
from wx_util import Timers
from mesh_gui import TriMesh, Camera, set_headlight, draw_mesh_faces
from helpers import Struct, solve_for_fixed_points

from vecutil import *
from OpenGL.GL import *
from glhelpers import *

from annotations_ellipsoid import kDebugging, kEps

kUseTexture = True
kOnlyShowAxis0 = True
## Screen-space (pixels) picking distance threshold
kPickScreenThreshold = 20

from annotations_ellipsoid import Annotations, AnnotationData, kTestNewOverlap

def intersection_curve_ellipsoid_ellipsoid( ellipsoid1, ellipsoid2 ):
    try:
        return XXX_intersection_curve_ellipsoid_ellipsoid( ellipsoid1, ellipsoid2 )
    except:
        debugger()
        raise
def XXX_intersection_curve_ellipsoid_ellipsoid( ellipsoid1, ellipsoid2 ):
    '''
    Given two Ellipsoids, returns points on the intersection curve between them.
    TODO: Points are not organized into a loop, they are structured in a not-so-useful manner.
    
    untested
    '''
    
    #debugger()
    
    ellipse1_bbox = ellipse_bbox( ( ellipsoid1.center, ellipsoid1.axes ) )
    ellipse2_bbox = ellipse_bbox( ( ellipsoid2.center, ellipsoid2.axes ) )
    both_bbox = bbox_intersection( ellipse1_bbox, ellipse2_bbox )
    
    maxx = both_bbox[0][0]
    maxy = both_bbox[0][1]
    minx = both_bbox[1][0]
    miny = both_bbox[1][1]
    
    look = asarray( ( 0,0,-1 ) )
    oo_accuracy = 1e3
    
    ## Sample a uniform grid inside the bounding box
    N = max( 1, int( round( oo_accuracy ) ) )
    ## We want to sample space uniformly.
    ## This means N samples should cover the area (maxx-minx) * (maxy-miny) uniformly.
    ## Let NX be the number of samples in x.
    ## Then we want NX * (maxy-miny)/(maxx-minx) samples in y
    ## and N = NX * NX * (maxy-miny)/(maxx-minx).
    ## So
    NX = max( 3, sqrt( N / (maxy-miny) * (maxx-minx) ) )
    NY = max( 3, NX * (maxy-miny)/(maxx-minx) )
    ## Round to the nearest integer.
    NX = int( round( NX ) )
    NY = int( round( NY ) )
    #print 'NX, NY:', NX, NY
    
    ## Add a sample to the beginning and end, so that we don't have to worry about
    ## the first and last values being intersection minima.
    xs = linspace( minx - 1./(NX-1), maxx + 1./(NX-1), NX+2 )
    ys = linspace( miny - 1./(NY-1), maxy + 1./(NY-1), NY+2 )
    
    ## 1 Collect the front/front, front/back, back/front, and back/back intersections' zdiffs.
    ##   These are solid areas of sampled intersections: the set intersection of ray intersections
    ##   with the front and back of the two ellipsoids.
    ## TODO Q: Should I return all four line strips?
    ## 2 Pick the two collections whose zdifferences are the smallest --- these are points on the
    ##   actual intersection curves.
    ## 3 Order each collection of points into a line strip.
    ## TODO: If the intersection curves are always quadrics, we only need 6 samples on each.
    
    ## 1
    front_front_xyzdiff = []
    front_back_xyzdiff = []
    back_front_xyzdiff = []
    back_back_xyzdiff = []
    for y in ys:
        for x in xs:
            ts1 = ray_ellipse_intersection_ts(
                ( ellipsoid1.center, ellipsoid1.axes ),
                ( (x,y,0.), look ),
                ray_is_line = True
                )
            ts2 = ray_ellipse_intersection_ts(
                ( ellipsoid2.center, ellipsoid2.axes ),
                ( (x,y,0.), look ),
                ray_is_line = True
                )
            
            ts1.sort()
            ts2.sort()
            assert len( ts1 ) in (0,1,2)
            assert len( ts2 ) in (0,1,2)
            
            try: front_front_xyzdiff.append( ( x, y, abs( ts1[0] - ts2[0] ) ) )
            except IndexError: pass
            try: front_back_xyzdiff.append( ( x, y, abs( ts1[0] - ts2[1] ) ) )
            except IndexError: pass
            try: back_front_xyzdiff.append( ( x, y, abs( ts1[1] - ts2[0] ) ) )
            except IndexError: pass
            try: back_back_xyzdiff.append( ( x, y, abs( ts1[1] - ts2[1] ) ) )
            except IndexError: pass
    
    ## 2
    all_xyzdiffs = []
    all_xyzdiffs.append( front_front_xyzdiff )
    all_xyzdiffs.append( front_back_xyzdiff )
    all_xyzdiffs.append( back_front_xyzdiff )
    all_xyzdiffs.append( back_back_xyzdiff )
    ## Degenerate intersection: cull empty points.
    ## Do this before we sort, because empty lists have sum() 0.
    all_xyzdiffs = [ xyzdiffs for xyzdiffs in all_xyzdiffs if len( xyzdiffs ) > 0 ]
    
    score2points = [
        ( sum( [ zdiff for x,y,zdiff in xyzdiffs ] ) / len( xyzdiffs ), xyzdiffs )
        for xyzdiffs in all_xyzdiffs
        ]
    assert len( score2points ) <= 4
    score2points.sort()
    print( [ score for score, points in score2points ] )
    del score2points[2:]
    collections = [ [ (x,y) for x,y,zdiff in points ] for score, points in score2points ]
    
    ## 3
    def order_line_strip( points ):
        '''
        Given a list of points from a line, order them along the line.
        '''
        #return points
        
        ## Assume the points can be roughly parameterized by an ellipse.
        center, axes = lsq_ellipse( points )
        center = asarray( center )
        axes = asarray( axes )
        ellipse = center, axes
        
        ellipse_points = [ closest_point_on_ellipse_to_point( ellipse, pt ) for pt in points ]
        thetas = ellipse_thetas_for_points( ellipse, ellipse_points )
        return [ pt for theta, pt in sorted( zip( thetas, points ) ) ]
    
    line_strips = [ order_line_strip( points ) for points in collections ]
    #return line_strips
    result = []
    for line_strip in line_strips: result = result + list( line_strip )
    return result
    
    
    
    ## The zdiff threshold will be the larger of the x,y step sizes.
    #ZDIFF_THRESHOLD = max( xs[1] - xs[0], ys[1] - ys[0] )
    ZDIFF_THRESHOLD = 1.
    
    def local_minima_indices( values ):
        ## Topologically correct values comparisons:
        ## break ties by index.
        def less( lhs, rhs ):
            return values[lhs] < values[rhs] or ( values[lhs] == values[rhs] and lhs < rhs )
        
        result = []
        for i in range( 1, len(values) - 1 ):
            if less( i, i-1 ) and less( i, i+1 ):
                result.append(i)
        return result
    
    ## Collect (x,y) points where the min zdiff at (x,y) is a local minima across all y.
    ## TODO: Collect these into a line loop.  Right now pairs of points are across the loop from
    ##       each other.
    result = []
    num_ones = 0
    for x in xs:
        ## Get a row of intersections.
        y_minzdiffs = x_y_minzdiff[ x ]
        ## The values we want are the local minima, if those minima are small enough.
        local_mins = local_minima_indices( [ y_minzdiff[1] for y_minzdiff in y_minzdiffs ] )
        ## Make a list of tuples: ( min zdiff, index )
        local_mins = [ ( y_minzdiffs[ local_min ][1], local_min ) for local_min in local_mins ]
        ## It's come up that we get three local mins, but one is much larger than the others.
        ## We can cull based on a zdiff threshold. Or we can keep the two smallest zdiffs.
        ## Cull values larger than the second-smallest zdiff or ZDIFF_THRESHOLD.
        row_threshold = ZDIFF_THRESHOLD
        #if len( local_mins ) > 2: debugger()
        if len( local_mins ) > 2: row_threshold = min( sorted( local_mins )[2][0], row_threshold )
        local_mins = [ val for val in local_mins if val[0] < row_threshold ]
        ## Convert local_mins back into a list of indices.
        local_mins = [ index for ( minzdiff, index ) in local_mins ]
        
        ## Collect x,y pairs of points in local_mins.
        xy_isectvals = [ ( x, y_minzdiffs[ local_min ][0] ) for local_min in local_mins ]
        
        ## Each row should have 0, 1, or 2 points on the intersection curve.
        assert len( xy_isectvals ) in (0,1,2)
        ## In general, there should be two points.  There can be 1 point twice, at the top and bottom.
        if len( xy_isectvals ) == 1:
            #expect( num_ones < 2 )
            num_ones += 1
        
        ## In general there are two points.  They are opposite each other in the intersection
        ## curve loop.  Put one at the beginning and one at the end.
        if len( xy_isectvals ) == 2:
            result.insert( 0, xy_isectvals[0] )
            result.append( xy_isectvals[1] )
        ## If there is only one point, put it at the end.
        ## This should happen twice, one for the first point and one for the last point.
        ## The first point will be appended to an empty list.
        ## Because it's a loop, the last point added at the end could also have been added at
        ## the beginning.
        elif len( xy_isectvals ) == 1:
            result.append( xy_isectvals[0] )
        elif len( xy_isectvals ) == 0: pass
        else:
            ## We shouldn't be here, but for some reason we get here sometimes.  Let's visualize it.
            assert len( xy_isectvals ) > 2
            result.insert( 0, xy_isectvals[0] )
            result.extend( xy_isectvals[1:] )
    
    ## There can be 1 point on the intersection curve twice, at the top and bottom.
    #expect( num_ones == 2 )
    return result

class CircleDrawer( EventHandler ):
    def __init__( self, ann, event2world, refresh_cb, failed_cb, success_cb, debugging = False ):
        
        self.stroke = []
        
        self.ann = ann
        self.ann.stroke = self.stroke
        self.event2world = event2world
        self.refresh_callback = refresh_cb
        self.success_callback = success_cb
        self.failed_callback = failed_cb
        
        self.debugging = debugging
    
    def doMouseDown( self, evt ):
        del self.stroke[:]
        self.stroke.append( self.event2world( evt.GetPosition() ) )
        if self.debugging: print(self.stroke[-1])
        return self
    
    def doMouseDrag( self, evt ):
        pt = self.event2world( evt.GetPosition() )
        
        if pt != self.stroke[-1]:
            self.stroke.append( pt )
            if self.debugging: print(self.stroke[-1])
            
            if pt != self.stroke[0]:
                self.ann.SetCircle( self.stroke )
                self.refresh_callback()
        
        return self
    
    def doMouseUp( self, evt ):
        self.doMouseDrag( evt )
        self.ann.stroke = None
        kTooNarrowAxesRatio = 50
        if self.ann.ellipse is None:
            self.failed_callback( self.ann )
        elif max( list(map( mag, self.ann.ellipse[1:] )) ) / min( list(map( mag, self.ann.ellipse[1:] )) ) > kTooNarrowAxesRatio:
            print('Ellipse too narrow!  Aborting.')
            self.failed_callback( self.ann )
        else:
            self.success_callback( self.ann )
        '''
        else:
            mags = [ mag( axis ) for axis in self.data.ellipse[1] ]
            ## If the ellipse is too big (our on-screen world is 2x2) or too eccentric, don't show it.
            if max( mags ) > 3. or 10. < max( mags ) / min( mags ):
                self.failed_callback( self.data )
        '''
        self.refresh_callback()
        
        unforward_events( self )
        return None

class HandleDragger( EventHandler ):
    def __init__( self, annotations, handle, event2world, refresh_cb, debugging = False ):
        
        self.stroke = []
        
        self.annotations = annotations
        self.handle = handle
        self.event2world = event2world
        self.refresh_callback = refresh_cb
        
        self.debugging = debugging
    
    def doMouseDown( self, evt ):
        del self.stroke[:]
        self.stroke.append( tuple( self.event2world( evt.GetPosition() ) ) )
        if self.debugging: print(self.stroke[-1])
        return self
    
    def doMouseDrag( self, evt ):
        pt = tuple( self.event2world( evt.GetPosition() ) )
        
        if pt != self.stroke[-1]:
            self.stroke.append( pt )
            if self.debugging: print(self.stroke[-1])
            
            self.handle.translate( asarray( self.stroke[-1] ) - asarray( self.stroke[-2] ) )
            self.refresh_callback()
        
        return self
    
    def doMouseUp( self, evt ):
        self.doMouseDrag( evt )
        
        self.annotations.PlaceAnnotationsInDepth()
        self.refresh_callback()
        
        unforward_events( self )
        return None
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        ## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
        ## in wxPython.
        if evt.GetKeyCode() in (wx.WXK_DELETE,8):
            print('deleting an annotation')
            self.annotations.remove( self.handle.annotation )
            self.annotations.PlaceAnnotationsInDepth()
            self.refresh_callback()
            unforward_events( self )
            return None
        else:
            evt.Skip()
        return self

class AnnotateCanvasCommon( BasicGLCanvas ):
    def __init__( self, parent, annotations, pos = wx.DefaultPosition, size = wx.DefaultSize, shared = None ):
        BasicGLCanvas.__init__( self, parent, pos, size, debugging = True, shared = shared )
        
        self.annotations = annotations
        self.show_annotations = True
        
        if kUseTexture:
            import os
            #self.annotations.SetImage( '/mixed media/pictures/ColossalKatamari.png' )
            self.annotations.SetImage( os.path.join( 'scenes', 'cactus.png' ) )
        
        self.debugging = False
    
    def GetBestSize( self ):
        return (500,500)
    
    def PickAxisLength( self, ann, pos ):
        '''
        Given an AnnotationData object 'ann' and a 2d screen-space position 'pos',
        returns the index of the AnnotationData.ellipsoid's axis whose length is displayed
        under the position 'pos'.  Returns -1 if no axis is picked (or there is no ellipsoid).
        '''
        
        assert ann is not None
        
        #debugger()
        
        pos = asarray( pos )
        centers = [ self.world2screen( ann.GetClosestAxisMidpointInDepth( i, (0,0,-1) ) ) for i in range(3) ]
        #print 'pos:', tuple( pos )
        #print 'centers:', map( tuple, centers )
        dists2 = [ mag2( pos - asarray( center ) ) for center in centers ]
        closest_axis = argmin( dists2 )
        if kOnlyShowAxis0: closest_axis = 0
        if dists2[ closest_axis ] < kPickScreenThreshold**2:
            return closest_axis
        else:
            return -1
        
        '''
        for i in (0,1,2):
            # TODO: if RectFromCornerAndSize( self.world2screen( self.data.ellipsoid.center + .5 * self.data.ellipsoid.axes[i] ), draw_text_pixel_size( ... ) ).Contains( pos ): return i
            if mag2(
                asarrayf( pos )
                - self.world2screen( self.data.ellipsoid.center + .5 * self.data.ellipsoid.axes[i] * sign( self.data.ellipsoid.axes[i][2] ) )
                ) < kPickScreenThreshold**2:
                return i
        
        return -1
        '''
    
    def PickAnnotationIndex( self, pos ):
        '''
        Given a 2d screen-space position 'pos', returns the index of the AnnotationData
        under the position 'pos'.  Returns -1 if no AnnotationData is picked.
        '''
        
        if len( self.annotations ) == 0: return -1
        
        world_pos = self.screen2world( pos )
        
        '''
        self.annotations.sortz()
        for i in range(len( self.annotations ))[::-1]:
            ann = self.annotations[i]
            
            ## Do we intersect the originally drawn ellipse?
            if ellipse_contains_point(
                ( ann.ellipse[0], ann.ellipse[1:] ),
                world_pos
                ):
                return i
        '''
        
        ## We'll shoot a ray under the mouse and pick the closest annotation.
        r = ( zerosf( 3 ), (0,0,-1) )
        r[0][:2] = world_pos
        
        t2annotation = []
        for ann_index, ann in enumerate( self.annotations ):
            ts = ray_ellipse_intersection_ts( ( ann.ellipsoid.center, ann.ellipsoid.axes ), r, ray_is_line = True )
            if len( ts ) > 0:
                t2annotation.append( ( min(ts), ann_index ) )
        
        if len( t2annotation ) > 0:
            t2annotation.sort()
            return t2annotation[0][1]
        
        return -1
    
    def PickAnnotation( self, pos ):
        '''
        Given a 2d screen-space position 'pos', returns the AnnotationData object
        under the position 'pos'.  Returns None if no AnnotationData is picked.
        '''
        
        ann_index = self.PickAnnotationIndex( pos )
        if -1 == ann_index: return None
        else: return self.annotations[ ann_index ]
    
    def PickAnnotationAxisTip( self, pos, limit_axes = (0,1,2) ):
        '''
        Given a 2d screen-space position 'pos', picks an ellipse axis near the position 'pos'.
        Returns the AnnotationData object, the index of the axis, and whether the positive (versus
        negative) endpoint along the axis was picked.
        If 'limit_axes' is specified, only axes in 'limit_axes' will be considered.
        Returns (None, undefined, undefined) if no AnnotationData's axis is picked.
        '''
        
        assert False not in [ i in (0,1,2) for i in limit_axes ]
        assert len( limit_axes ) == len( set( limit_axes ) )
        
        pos = asarray( pos )[:2]
        
        dist2result = []
        self.annotations.sortz()
        for ann in reversed( self.annotations ):
            for ai in limit_axes:
                dist2result.append( (
                    mag2(
                        pos - asarray( self.world2screen( ann.ellipsoid.center + ann.ellipsoid.axes[ai] ) )[:2]
                        ),
                    ( ann, ai, True )
                    ) )
                dist2result.append( (
                    mag2(
                        pos - asarray( self.world2screen( ann.ellipsoid.center - ann.ellipsoid.axes[ai] ) )[:2]
                        ),
                    ( ann, ai, False )
                    ) )
        
        dist2result.sort()
        if len( dist2result ) == 0 or dist2result[0][0] >= kPickScreenThreshold**2:
            return None, None, None
        else:
            return dist2result[0][1]
    
    def prepare2D( self ):
        glDisable( GL_DEPTH_TEST )
        glDisable( GL_LIGHTING )
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        ## Center the square inside the viewport.
        glOrtho(*( list( self.world_bounds ) + [ -100, 100 ] ))
		
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()


class ModeProxy( EventHandler ):
    def __init__( self, canvas ):
        self.canvas = canvas
    
    ## This magic allows non-overloaded member names to magically pass through to self.canvas.
    def __getattr__( self, name ):
        if not name.startswith('_'): return getattr( self.canvas, name )
        else: raise AttributeError("I am not proxying attribute '%s'" % (name,))

class OldAnnotate2Mode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def doMouseDown( self, evt ):
        ## If we pick a handle, manipulate it.
        handle = self.PickHandle( self.event2screen( evt.GetPosition() ) )
        if handle is not None:
            forward_events( self, HandleDragger(
                self.annotations,
                handle,
                self.event2world,
                self.GetTopLevelParent().Refresh,
                debugging = self.debugging
                ) )
        ## Otherwise, so long as we didn't pick an annotation, draw a new one.
        elif self.PickAnnotation( self.event2screen( evt.GetPosition() ) ) is None:
            self.annotations.append( AnnotationData(
                z = self.annotations.MaxZ()
                ) )
            forward_events( self, CircleDrawer(
                self.annotations[-1],
                self.event2world,
                self.GetTopLevelParent().Refresh,
                lambda ann: self.annotations.remove( ann ),
                self.NewAnnotationCreated,
                debugging = self.debugging
                ) )
        ## If we're doing nothing, don't call self.doMouseDown() again.
        else: return self
        
        ## Now that we've forwarded events, call doMouseDown() again.  This time it will forward.
        return self.doMouseDown( evt )
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
        ## in wxPython.
        if keycode in (wx.WXK_DELETE,8):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            ## If the mouse is not over any annotations, check if we're over a handle
            ## (they protrude a little bit, and this is consistent with delete-while-dragging
            ## a handle).
            try:
                ann = self.PickAnnotation( screen_pos ) or self.PickHandle( screen_pos ).annotation
                print('deleting an annotation')
                self.annotations.remove( ann )
                
                self.annotations.PlaceAnnotationsInDepth()
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        elif keycode in (wx.WXK_DOWN,wx.WXK_UP):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## Get the annotation or handle's annotation under the mouse.
            try:
                ann = self.PickAnnotation( screen_pos ) or self.PickHandle( screen_pos ).annotation
                ## Modify its z
                z_delta = { wx.WXK_DOWN: -1, wx.WXK_UP: 1 }[ keycode ]
                ann.z += z_delta
                #if ann.parentage: ann.parentage.front_sign = z_delta
                
                self.annotations.PlaceAnnotationsInDepth()
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        else:
            evt.Skip()
        
        return self
    
    def doDraw( self ):
        self.DrawScene( draw_handles = True )
    
    def doMouseMove( self, evt ):
        ## I have mouse-over (hover) feedback.
        self.Refresh()
        return self

class AnnotationCreateMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def doMouseDown( self, evt ):
        ## Create a new annotation.
        self.annotations.append( AnnotationData(
            z = self.annotations.MaxZ()
            ) )
        ## Forward events to the class which creates it.
        forward_events( self, CircleDrawer(
            self.annotations[-1],
            self.event2world,
            self.GetTopLevelParent().Refresh,
            lambda ann: self.annotations.remove( ann ),
            self.NewAnnotationCreated,
            debugging = self.debugging
            ) )
        ## Now that we've forwarded events, call doMouseDown() again.  This time it will forward.
        return self.doMouseDown( evt )
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
        ## in wxPython.
        if keycode in (wx.WXK_DELETE,8):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            ## If the mouse is not over any annotations, check if we're over a handle
            ## (they protrude a little bit, and this is consistent with delete-while-dragging
            ## a handle).
            try:
                ann = self.PickAnnotation( screen_pos ) or self.PickHandle( screen_pos ).annotation
                print('deleting an annotation')
                self.annotations.remove( ann )
                
                self.annotations.PlaceAnnotationsInDepth()
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        else:
            evt.Skip()
        
        return self
    
    def doDraw( self ):
        self.DrawScene( draw_handles = False )

class HandleManipMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def doMouseDown( self, evt ):
        ## If we pick a handle, manipulate it.
        handle = self.PickHandle( self.event2screen( evt.GetPosition() ) )
        if handle is not None:
            forward_events( self, HandleDragger(
                self.annotations,
                handle,
                self.event2world,
                self.GetTopLevelParent().Refresh,
                debugging = self.debugging
                ) )
            ## Now that we've forwarded events, call doMouseDown() again.  This time it will forward.
            return self.doMouseDown( evt )
        
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
            ## If the mouse is not over any annotations, check if we're over a handle
            ## (they protrude a little bit, and this is consistent with delete-while-dragging
            ## a handle).
            try:
                ann = self.PickAnnotation( screen_pos ) or self.PickHandle( screen_pos ).annotation
                print('deleting an annotation')
                self.annotations.remove( ann )
                
                self.annotations.PlaceAnnotationsInDepth()
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        elif keycode in (wx.WXK_DOWN,wx.WXK_UP):
            
            ## evt.GetPosition() returns garbage.
            screen_pos = self.event2screen( self.GetMousePosition() )
            
            ## Get the annotation or handle's annotation under the mouse.
            try:
                ann = self.PickAnnotation( screen_pos ) or self.PickHandle( screen_pos ).annotation
                ## Modify its z
                z_delta = { wx.WXK_DOWN: -1, wx.WXK_UP: 1 }[ keycode ]
                ann.z += z_delta
                #if ann.parentage: ann.parentage.front_sign = z_delta
                
                self.annotations.PlaceAnnotationsInDepth()
                self.GetTopLevelParent().Refresh()
                return self
            except: pass
        
        else:
            evt.Skip()
        
        return self
    
    def doDraw( self ):
        self.DrawScene( draw_handles = True )
    
    def doMouseMove( self, evt ):
        ## I have mouse-over (hover) feedback.
        self.Refresh()
        return self

class ParentageMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        
        self.mouse_down_ann = None
        self.mouse_down_world_pos = None
    
    def doMouseDown( self, evt ):
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return self
        
        self.mouse_down_ann = self.PickAnnotation( self.event2screen( evt.GetPosition() ) )
        self.mouse_down_world_pos = self.event2world( evt.GetPosition() )
        return self
    
    def doMouseUp( self, evt ):
        mouse_up_ann = self.PickAnnotation( self.event2screen( evt.GetPosition() ) )
        
        ## If we picked an annotation on mouse down and mouse up and they're different and they
        ## overlap, parent the smaller one with the larger one.
        if( self.mouse_down_ann is not None
            and
            mouse_up_ann is not None
            and
            self.mouse_down_ann is not mouse_up_ann
            and
            ## They must overlap or else being the parent isn't so meaningful.
            self.mouse_down_ann.Overlaps( mouse_up_ann )
            and
            ## Make sure they are not already parented to each other.
            (
                ( self.mouse_down_ann.parentage is None or self.mouse_down_ann.parentage.parent is not mouse_up_ann )
                and
                ( mouse_up_ann.parentage is None or mouse_up_ann.parentage.parent is not self.mouse_down_ann )
            )
            ):
            ## Make the larger one the parent of the smaller one.
            size_ann = [
                ( ellipse_volume( self.mouse_down_ann.ellipse[1:] ), self.mouse_down_ann ),
                ( ellipse_volume( mouse_up_ann.ellipse[1:] ), mouse_up_ann )
                ]
            size_ann.sort()
            ann = size_ann[0][1]
            parent = size_ann[1][1]
            
            ann.parentage = Struct()
            ann.parentage.parent = parent
            ## This element is a 2-element tuple of integers \in {0,1}.
            ## The first element determines whether the attachment point
            ## is ann's front or back (0 means front, 1 means back).
            ## The second element determines whether the attachment point
            ## is ann.parent's front or back (0 means front, 1 means back).
            ## By default, make ann's front attach to its parent's front.
            ann.parentage.ann_parent_front_back_signature = (0,0)
            ann.parentage.center_offset = (
                center_of_mass_of_ellipses_intersection(
                    ( ann.ellipse[0], ann.ellipse[1:] ),
                    ( parent.ellipse[0], parent.ellipse[1:] )
                    )
                -
                ann.ellipse[0]
                )
            
            ## Place all annotations in depth
            self.annotations.PlaceAnnotationsInDepth()
            self.GetTopLevelParent().Refresh()
        
        elif self.mouse_down_ann is not None and mouse_up_ann is None:
            ## Destroy the parentage relationship of mouse_down_ann.
            self.mouse_down_ann.parentage = None
            
            ## Place all annotations in depth
            self.annotations.PlaceAnnotationsInDepth()
            self.GetTopLevelParent().Refresh()
        
        self.mouse_down_ann = None
        self.mouse_down_world_pos = None
        return self
    
    def doMouseDrag( self, evt ):
        ## We draw while dragging, so refresh.
        self.Refresh()
        return self
    
    def doContextMenu( self, evt ):
        evt_pos = self.GetMousePosition()
        print('context menu:', self.event2world( evt_pos ))
        
        ## This is the list of annotations, ordered from back to front.
        picked = []
        
        world_pos = self.event2world( evt_pos )
        self.annotations.sortz()
        for ann in self.annotations:
            ## Do we intersect the originally drawn ellipse?
            if ellipse_contains_point(
                ( ann.ellipse[0], ann.ellipse[1:] ),
                world_pos
                ):
                picked.append( ann )
        
        print('context menu picked:', picked)
        
        ## Show a contextual menu for the closest two overlapping annotations.
        del picked[ :-2 ]
        if len( picked ) < 2:
            evt.Skip()
            return self
        
        menu = wx.Menu()
        menu_ids = []
        
        ## Are they already connected?  If so, add a disconnect option.
        if( picked[0].parentage and picked[0].parentage.parent is picked[1]
            or
            picked[1].parentage and picked[1].parentage.parent is picked[0]
            ):
            
            def disconnect( evt ):
                self.annotations.RemoveEdge( picked )
                
                self.annotations.PlaceAnnotationsInDepth()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Disconnect" )
            self.Bind( wx.EVT_MENU, disconnect, id = item.GetId() )
            menu_ids.append( item.GetId() )
            menu.AppendSeparator()
        
        def connect_with_signature( signature ):
            if picked[0].parentage and picked[0].parentage.parent is picked[1]:
                picked[0].parentage.ann_parent_front_back_signature = signature
            elif picked[1].parentage and picked[1].parentage.parent is picked[0]:
                picked[1].parentage.ann_parent_front_back_signature = signature
            else:
                ## We attach the smaller one to the larger one.
                size_ann = [
                    ( ellipse_volume( picked[0].ellipse[1:] ), picked[0] ),
                    ( ellipse_volume( picked[1].ellipse[1:] ), picked[1] )
                    ]
                size_ann.sort()
                ann = size_ann[0][1]
                parent = size_ann[1][1]
                ## If the smaller one is already attached, attach the larger one to the smaller one.
                if ann.parentage is not None:
                    ann, parent = parent, ann
                ## If they're both already attached elsewhere, we can't attach them without creating
                ## a cycle!
                if ann.parentage is not None:
                    print('ERROR: Refusing to attach, a cycle would be created.')
                    return
                
                ann.parentage = Struct()
                ann.parentage.parent = parent
                ## This element is a 2-element tuple of integers \in {0,1}.
                ## The first element determines whether the attachment point
                ## is ann's front or back (0 means front, 1 means back).
                ## The second element determines whether the attachment point
                ## is ann.parent's front or back (0 means front, 1 means back).
                ## By default, make ann's front attach to its parent's front.
                ann.parentage.ann_parent_front_back_signature = signature
                ann.parentage.center_offset = (
                    center_of_mass_of_ellipses_intersection(
                        ( ann.ellipse[0], ann.ellipse[1:] ),
                        ( parent.ellipse[0], parent.ellipse[1:] )
                        )
                    -
                    ann.ellipse[0]
                    )
            
            self.annotations.PlaceAnnotationsInDepth()
            self.GetTopLevelParent().Refresh()
        
        def connect_front_front( evt ):
            print('connect_front_front')
            connect_with_signature( (0,0) )
        def connect_back_back( evt ):
            print('connect_back_back')
            connect_with_signature( (1,1) )
        
        item = menu.Append( wx.ID_ANY, "Connect, Front" )
        self.Bind( wx.EVT_MENU, connect_front_front, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        item = menu.Append( wx.ID_ANY, "Connect, Back" )
        self.Bind( wx.EVT_MENU, connect_back_back, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        
        self.PopupMenu( menu )
        menu.Destroy()
        for id in menu_ids: self.Unbind( wx.EVT_MENU, id = id )
        
        return self
    
    def doDraw( self ):
        def draw_dragging_line():
            if self.mouse_down_world_pos is None: return
            
            depth_test = glIsEnabled( GL_DEPTH_TEST )
            if depth_test: glDisable( GL_DEPTH_TEST )
            
            ## A thick, orange line.
            glLineWidth( 6 )
            glColor4f( 1, .5, 0, .8 )
            draw_linestrip( [ self.mouse_down_world_pos, self.event2world( self.GetMousePosition() ) ] )
            draw_linestrip( solve_for_fixed_points(
                [ (-1, -1), (0,0), (1,-1) ],
                ( (0,0), (0,10) ),
                ( self.event2world( self.GetMousePosition() ), asarray( self.event2world( self.GetMousePosition() ) ) + dir( asarray(self.event2world( self.GetMousePosition() )) - asarray(self.mouse_down_world_pos) ) )
                ) )
            
            if depth_test: glEnable( GL_DEPTH_TEST )
        
        self.DrawScene( draw_handles = True, callback = draw_dragging_line )

def PickHandle( annotations, world_pos, kind2world_radius ):
    '''
    Given a sequence of annotations, a 2d world-space position 'pos',
    and a function mapping AnnotationData.*Handle types to world-space radii,
    picks a Handle near the position 'pos'.
    Returns the picked Handle object, or None if no Handle axis is picked.
    '''
    
    z2handle = []
    for ann in annotations:
        handle_z = PickAnnotationHandleZ( ann, world_pos, kind2world_radius )
        if handle_z is not None: z2handle.append( ( handle_z[1], handle_z[0] ) )
    
    if len( z2handle ) == 0:
        return None
    else:
        return max( z2handle )[1]

def PickAnnotationHandle( ann, world_pos, kind2world_radius ):
    '''
    Given an AnnotationData object 'ann', a 2d world-space position 'pos',
    and a function mapping AnnotationData.*Handle types to world-space radii,
    picks a Handle near the position 'pos'.
    Returns the picked Handle object, or None if no Handle axis is picked.
    '''
    
    handle_z = PickAnnotationHandleZ( ann, world_pos, kind2world_radius )
    if handle_z is None:
        return handle_z
    else:
        return handle_z[0]

def PickAnnotationHandleZ( ann, world_pos, kind2world_radius ):
    '''
    Given an AnnotationData object 'ann', a 2d world-space position 'pos',
    and a function mapping AnnotationData.*Handle types to world-space radii,
    picks a Handle near the position 'pos'.
    Returns the tuple ( picked Handle object, picked Handle's z ),
    or None if no Handle axis is picked.
    '''
    
    pos = asarrayf( world_pos ).ravel()[:2]
    kind2radius = kind2world_radius
    
    z2result = []
    whole_ellipse_handle_z = None
    for handle in ann.Handles:
        if handle.kind == AnnotationData.WholeEllipseHandle:
            r = ( append( pos, 0. ), asarray( (0,0,-1) ) )
            ts = ray_ellipse_intersection_ts( ( handle.annotation.ellipsoid.center, handle.annotation.ellipsoid.axes ), r, ray_is_line = True )
            ## If we don't intersect, return our guess
            if len( ts ) > 0:
                whole_ellipse_handle_z = ( handle, r[0][2] + min(ts) * r[1][2] )
            continue
        
        handle_pos = handle.get_pos()
        if( handle.kind == AnnotationData.Axis0Handle and
            handle_pos[2] + kEps < handle.annotation.ellipsoid.center[2]
            ):
            continue
        
        dist2 = mag2( pos - asarray( handle_pos )[:2] )
        if dist2 > kind2radius( handle.kind )**2: continue
        
        z2result.append( ( handle_pos[2], handle ) )
    
    ## If we didn't pick any ordinary handles, return the whole ellipse handle
    ## (which may be None).
    if len( z2result ) == 0: return whole_ellipse_handle_z
    
    ## This chooses the handle under the mouse closest to the camera.
    max_z2result = max( z2result )
    return max_z2result[1], max_z2result[0]

def PickTAnnotations( annotations, world_pos, world_dir = (0,0,-1) ):
    '''
    Given a sequence of AnnotationData 'annotations', a 2+d world-space position 'world_pos',
    and optional world-space vector 'world_dir' (defaults to -z),
    returns a list of unsorted (t, annotation) pairs corresponding to intersections of the
    world-ray with annotations.
    Note that each annotation will appear at most once, with the smallest 't'.
    Note that t's may be negative.
    '''
    
    ## We'll shoot a ray under the mouse and pick the closest annotation.
    r = ( zerosf( 3 ), world_dir )
    r[0][:2] = world_pos
    
    t2annotations = []
    for ann in annotations:
        ts = ray_ellipse_intersection_ts( ( ann.ellipsoid.center, ann.ellipsoid.axes ), r, ray_is_line = True )
        if len( ts ) > 0:
            t2annotations.append( ( min(ts), ann ) )
    
    return t2annotations

class AnnotateCanvas2( AnnotateCanvasCommon ):
    
    kAxis0HandleRadius = .03
    kAxis1HandleRadius = .02
    kCrossSectionHandleRadius = .02
    ## Dynamically generated since we use pixel-sized text.
    # kDiscreteZHandleRadius =
    kParentIntersectionHandle = .02
    
    def __init__( self, *args, **kwargs ):
        AnnotateCanvasCommon.__init__( self, *args, **kwargs )
        
        #self.SetMode( OldAnnotate2Mode( self ) )
        self.SetMode( AnnotationCreateMode( self ) )
        #self.SetMode( HandleManipMode( self ) )
    
    def SetMode( self, mode ):
        ## Oops: I can't use forward_events because it's broken and only works to one level.
        ## My modes need to be able to forward events, too.
        #if self.mode is not None: unforward_events( self.mode )
        self.mode = mode
        #forward_events( self, self.mode )
    
    kOldAnnotate2Mode = 1
    kAnnotationCreateMode = 2
    kHandleManipMode = 3
    kParentageMode = 4
    
    @staticmethod
    def GetAllModeIDs():
        return [
            AnnotateCanvas2.kOldAnnotate2Mode,
            AnnotateCanvas2.kAnnotationCreateMode,
            AnnotateCanvas2.kHandleManipMode,
            AnnotateCanvas2.kParentageMode
            ]
    
    @staticmethod
    def GetLongNameForModeID( ID ):
        return {
            AnnotateCanvas2.kOldAnnotate2Mode: 'Old Annotation 2 Mode',
            AnnotateCanvas2.kAnnotationCreateMode: 'Create shapes',
            AnnotateCanvas2.kHandleManipMode: 'Manipulate annotations',
            AnnotateCanvas2.kParentageMode: 'Connect shapes'
            }[ ID ]
    
    @staticmethod
    def GetShortNameForModeID( ID ):
        return {
            AnnotateCanvas2.kOldAnnotate2Mode: 'Annotation2',
            AnnotateCanvas2.kAnnotationCreateMode: 'Create',
            AnnotateCanvas2.kHandleManipMode: 'Manipulate',
            AnnotateCanvas2.kParentageMode: 'Connect'
            }[ ID ]
    
    def SetModeByID( self, ID ):
        mode_type = {
            AnnotateCanvas2.kOldAnnotate2Mode: OldAnnotate2Mode,
            AnnotateCanvas2.kAnnotationCreateMode: AnnotationCreateMode,
            AnnotateCanvas2.kHandleManipMode: HandleManipMode,
            AnnotateCanvas2.kParentageMode: ParentageMode
            }[ ID ]
        print('Switching to mode:', self.GetLongNameForModeID( ID ))
        self.SetMode( mode_type( self ) )
    
    def PickHandle( self, screen_pos ):
        '''
        Given a sequence of annotations and a 2d screen-space position 'pos',
        picks a Handle near the position 'pos'.
        Returns the picked Handle object, or None if no Handle axis is picked.
        '''
        
        world_pos = self.screen2world( screen_pos )
        kind2world_radius = {
            AnnotationData.Axis0Handle: self.kAxis0HandleRadius,
            AnnotationData.Axis1Handle: self.kAxis1HandleRadius,
            AnnotationData.CrossSectionHandle: self.kCrossSectionHandleRadius,
            AnnotationData.ParentIntersectionHandle: self.kParentIntersectionHandle,
            #AnnotationData.DiscreteZHandle: self.kAxis0HandleRadius #10 / 500. #mag( asarray( self.screen2world( (20,20) ) ) )/2
            #AnnotationData.DiscreteZHandle: max( asarray( draw_text_world_size( self.world2screen, self.screen2world, 'Z: ' + str( ann.z ) ) ) )/2
            }
        
        return PickHandle( self.annotations, world_pos, lambda key: kind2world_radius[key] )
    
    def NewAnnotationCreated( self, ann ):
        ## Calling this prevents misbehavior the first time an axis0 handle is dragged.
        ann.SetCrossSectionFromStroke( 1, True, [ ann.ellipsoid.center[:2] + ann.ellipsoid.axes[2][:2] ], (0,0,-1) )
        
        ## Set the z of the annotation to be just in front of the shape it most overlaps
        overlaps = [
            (
                ellipses_overlap_area(
                    ( ann.ellipse[0], ann.ellipse[1:] ),
                    ( other_ann.ellipse[0], other_ann.ellipse[1:] )
                    ),
                other_ann
            )
            for other_ann in self.annotations
            if ann is not other_ann
            ]
        
        # print 'overlaps:', overlaps
        ann.parentage = None
        try:
            max_overlap = max( overlaps )
            #print 'max_overlap:', max_overlap
            if max_overlap[0] > 0:
                ann.z = max_overlap[1].z + 1
                ## Store a struct of 'parent', 'ann_parent_front_back_signature'
                ## (a tuple indicating whether ann's front or back is attached
                ## to its parent's front or back), and 'center_offset' (the offset
                ## from ann's center to a point on the intersection curve).
                ## 'ann_parent_front_back_signature' is a 2-element tuple of integers \in {0,1}.
                ## The first element determines whether the attachment point
                ## is ann's front or back (0 means front, 1 means back).
                ## The second element determines whether the attachment point
                ## is ann.parent's front or back (0 means front, 1 means back).
                ## By default, make ann's front attach to its parent's front.
                ## TODO: There are four combinations: ann.front is on parent's front,
                ##       ann.back is on parent's front, ann.front is on parent's back,
                ##       ann.back is on parent's back.
                ##       On the other hand, ann.front on parent.front covers nearly all of
                ##       ann.back on parent.front, so maybe I can ignore that case.
                ##       Also, perhaps I can leave precise parent.back placement to symmetry constraints.
                ##       Alternatively, I can show the curve (dotted if obscured), and make the
                ##       curve draggable.  If the back curve is dragged, then the back touches the front.
                parent = max_overlap[1]
                ann.parentage = Struct()
                ann.parentage.parent = parent
                ann.parentage.ann_parent_front_back_signature = (0,0)
                ann.parentage.center_offset = (
                    center_of_mass_of_ellipses_intersection(
                        ( ann.ellipse[0], ann.ellipse[1:] ),
                        ( parent.ellipse[0], parent.ellipse[1:] )
                        )
                    -
                    ann.ellipse[0]
                    )
        except ValueError: pass
        
        ## Place all annotations in depth
        self.annotations.PlaceAnnotationsInDepth()
    
    
    def DrawAnnotation_old( self, ann, picked_handle = None, draw_handles = True ):
        glDisable( GL_LIGHTING )
        
        ## Project the annotation so that all of its z values are equal to ann.z
        ## (project onto the z = ann.z plane).
        glMatrixMode( GL_MODELVIEW )
        glPushMatrix()
        glMultMatrixf( (1,0,0,0, 0,1,0,0, 0,0,0,0, 0,0,ann.z,1) )
        
        ## If the annotation still has a stroke, it's currently being drawn.
        ## Draw the stroke and return.
        if ann.stroke is not None:
            glLineWidth( 8 )
            glColor4ub( 174, 174, 0, 100 )
            draw_linestrip( ann.stroke )
            
            glPopMatrix()
            return
        
        assert ann.ellipse is not None and ann.ellipsoid is not None
        
        ## Draw circles spanning the axes.
        #elif ann.ellipsoid is not None:
        widths = [ 2, 2, 1 ]
        cols = [ (50,50,100,200), (174,0,174,200), (50,50,100,100) ]
        start_angles = [ 0, 0, 0 ]
        end_angles = [ pi, 2*pi, 2*pi ]
        start_angles[2], end_angles[2] = (-pi/2, pi/2) if ann.ellipsoid.axes[0][2] > -kEps else (pi/2, 3*pi/2)
        
        ## Draw the 2d ellipse.
        ## We can't draw rely on drawing the ellipse spanning perpendiculars to axis 2 (below)
        ## because the skeleton axis is chosen to keep the 2d ellipse constant.
        ## Draw the outline in black.
        glLineWidth( 2 )
        glColor4ub( 0,0,0,200 )
        draw_ellipse( *ann.ellipse )
        
        ## 0 is the cross section
        ## 1 is the major axis (looks like the skeleton)
        ## 2 is the silhouette (only if the cross section isn't tilted, otherwise it's skewed)
        if kDebugging:
            ellipse_axes = (0,1,2)
        else:
            ellipse_axes = (0,2)
        for i in ellipse_axes:
            glLineWidth( widths[i] )
            glColor4ub( *cols[i] )
            draw_ellipse(
                ann.ellipsoid.center,
                ann.ellipsoid.axes[ (i+1) % 3 ],
                ann.ellipsoid.axes[ (i+2) % 3 ],
                start_radians = start_angles[i],
                end_radians = end_angles[i]
                )
        
        
        ## Draw the handles.
        handles = ann.Handles if draw_handles else []
        for handle in handles:
            fill_alpha = 1 if handle == picked_handle else .5
            
            if handle.kind == ann.Axis0Handle:
                if handle.get_pos()[2] + kEps < handle.annotation.ellipsoid.center[2]: continue
                
                ## Blue handles at either end of the ellipse's major axis.
                glColor4f( 0,0,.5, fill_alpha )
                draw_circle_interior( handle.get_pos(), AnnotateCanvas2.kAxis0HandleRadius )
                glLineWidth( 2 )
                glColor3f( 0,0,.5 )
                draw_circle( handle.get_pos(), AnnotateCanvas2.kAxis0HandleRadius )
            elif handle.kind == ann.Axis1Handle:
                ## Light red handles at either end of the ellipse's minor axis.
                glColor4f( .5,.25,.25, fill_alpha )
                draw_circle_interior( handle.get_pos(), AnnotateCanvas2.kAxis1HandleRadius )
                glLineWidth( 2 )
                glColor3f( .5,.25,.25 )
                draw_circle( handle.get_pos(), AnnotateCanvas2.kAxis1HandleRadius )
            elif handle.kind == ann.CrossSectionHandle:
                ## A handle for the cross section.
                glColor4f( .5,.25,.5, fill_alpha )
                draw_circle_interior( handle.get_pos(), AnnotateCanvas2.kCrossSectionHandleRadius )
                glLineWidth( 2 )
                glColor3f( .5,.25,.5 )
                draw_circle( handle.get_pos(), AnnotateCanvas2.kCrossSectionHandleRadius )
            elif handle.kind == ann.ParentIntersectionHandle:
                ## A handle for the intersection between this annotation and its parent.
                glColor4f( 1,.5,0, fill_alpha )
                draw_circle_interior( handle.get_pos(), AnnotateCanvas2.kParentIntersectionHandle )
                glLineWidth( 2 )
                glColor3f( 1,.5,0 )
                draw_circle( handle.get_pos(), AnnotateCanvas2.kParentIntersectionHandle )
            elif handle.kind == ann.DiscreteZHandle:
                ## A handle for the discrete z depth.
                ## TODO: An arrow which indicates that the annotation at the start is in front of
                ##       the annotation at the back.
                glColor4f( 0,0,0, fill_alpha )
                txt = 'Z: ' + str( ann.z )
                world_size = append( draw_text_world_size( self.world2screen, self.screen2world, txt ), 0. )
                draw_text(
                    asarray( handle.get_pos() ) - .5 * asarray( world_size ),
                    txt
                    )
            elif handle.kind == ann.WholeEllipseHandle:
                pass
            else:
                print('Unknown handle kind:', handle.kind)
        
        if ann.parentage is not None: #and kDebugging:
            from surface_intersections import ellipsoid_ellipsoid_intersection_lines
            line_segments = ellipsoid_ellipsoid_intersection_lines(
                ann.ellipsoid,
                ann.parentage.parent.ellipsoid
                )
            glColor4f( 0, .5, 1., .75 )
            glLineWidth( 2 )
            glBegin( GL_LINES )
            for line in line_segments:
                glVertex3f( *line[0] )
                glVertex3f( *line[1] )
            glEnd()
            
            '''
            ## Shoot a bunch of rays through this annotation and its parent.
            pts = intersection_curve_ellipsoid_ellipsoid(
                ann.ellipsoid,
                ann.parentage.parent.ellipsoid
                )
            glColor4f( 0, .5, 1., .75 )
            #glLineWidth( 2 )
            #draw_lineloop( pts )
            glPointSize( 2 )
            glBegin( GL_POINTS )
            for pt in pts: glVertex2f( *pt )
            glEnd()
            '''
        
        glMatrixMode( GL_MODELVIEW )
        glPopMatrix()
    
    def DrawScene_old( self, draw_handles = True, callback = None ):
        glClearColor( .9,.9,.9,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        self.prepare2D()
        
        ## Draw the gltex if there is one.
        if self.annotations.gltex is not None:
            glColor4f( 1,1,1, .5 )
            draw_gltexture( self.annotations.gltex, self.annotations.tex_pos, self.annotations.tex_size )
        
        if self.show_annotations:
            glEnable( GL_DEPTH_TEST )
            
            ## Sort by z so we can draw the annotations in depth order
            self.annotations.sortz()
            
            ## Draw the annotations as depth masks, front to back.
            ## Draw them white, since they are also the background of the annotations.
            glColor4f( 1,1,1, .5 )
            for ann in self.annotations[::-1]:
                ## Skip in-progress strokes
                if ann.stroke is not None: continue
                
                assert ann.ellipse is not None
                
                ## Draw a white background for the interior, which also sets the depth buffer
                ## (we set the color at the beginning of this loop).
                draw_ellipse_interior( append( ann.ellipse[0], ann.z - .1 ), append( ann.ellipse[1], 0 ), append( ann.ellipse[2], 0 ) )
            
            ## Now turn off drawing into the depth mask, but continue to compare against it.
            glDepthMask( GL_FALSE )
            
            ## Draw the annotations' shadows and white backgrounds, back to front
            for ann in self.annotations:
                if ann.stroke is not None: continue
                
                ## Draw the shadow.
                shadow_scale = 1.
                shadow_offset = .05
                smaller_radius1 = mag( ann.ellipse[1] )
                smaller_radius2 = mag( ann.ellipse[2] )
                larger_radius1 = shadow_scale * smaller_radius1 + shadow_offset
                larger_radius2 = shadow_scale * smaller_radius2 + shadow_offset
                
                ## The shadow lives below ann.z (and its depth mask) but above the next-lower layer.
                glColor( 1,0,0,1 )
                draw_ellipse_crust(
                    ( ann.ellipse[0][0], ann.ellipse[0][1], ann.z - .5 ),
                    ann.ellipse[1], ann.ellipse[2],
                    larger_radius1, larger_radius2,
                    smaller_radius1, smaller_radius2,
                    larger_color = (0,0,0, 0),
                    smaller_color = (0,0,0, .5)
                    )
                
                ## Draw a white background for the interior.
                ## UPDATE: Actually, we do this above, when we draw into the z buffer.
                #glColor4f( 1,1,1, .5 )
                #draw_ellipse_interior( append( ann.ellipse[0], ann.z ), append( ann.ellipse[1], 0 ), append( ann.ellipse[2], 0 ) )
            
            ## Draw the annotations, back to front
            picked_handle = None
            if draw_handles:
                picked_handle = self.PickHandle( self.event2screen( self.GetMousePosition() ) )
                #print 'picked:', picked_handle
            for ann in self.annotations: self.DrawAnnotation( ann, picked_handle = picked_handle, draw_handles = draw_handles )
            
            ## Turn the depth mask back on (normal OpenGL behavior)
            glDepthMask( GL_TRUE )
        
        if callback is not None: callback()
        
        self.SwapBuffers()
    
    def DrawAnnotation( self, ann, picked_handle = None, draw_handles = True ):
        glDisable( GL_LIGHTING )
        
        ## If the annotation still has a stroke, it's currently being drawn.
        ## Draw the stroke and return.
        if ann.stroke is not None:
            glDisable( GL_DEPTH_TEST )
            
            glLineWidth( 8 )
            glColor4ub( 174, 174, 0, 100 )
            draw_linestrip( ann.stroke )
            
            glEnable( GL_DEPTH_TEST )
            return
        
        assert ann.ellipse is not None and ann.ellipsoid is not None
        
        ## Draw the silhouette.  We have to do this before we draw the object itself.
        ## Also, we can't draw rely on drawing the ellipse spanning perpendiculars to axis 2 (below)
        ## because the skeleton axis is chosen to keep the 2d ellipse constant.
        ## Draw the silhouette in black.
        glColor4ub( 0,0,0, 200 )
        glLineWidth( 4 )
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
        glEnable( GL_CULL_FACE )
        glCullFace( GL_FRONT )
        ann.ellipsoid.Draw()
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL )
        glCullFace( GL_BACK )
        glDisable( GL_CULL_FACE )
        
        ## Draw the 3d ellipse solid white.
        ## TODO: Draw the ellipsoid transparent white.
        glColor3f( 1,1,1 )
        ## This offset guarantees that all polygon faces have a z-buffer value at least 1 greater.
        glPolygonOffset( 1.0, 1.0 )
        glEnable( GL_POLYGON_OFFSET_FILL )
        ann.ellipsoid.Draw()
        glDisable( GL_POLYGON_OFFSET_FILL )
        
        ## Draw circles spanning the axes.
        #elif ann.ellipsoid is not None:
        widths = [ 2, 2, 1 ]
        cols = [ (50,50,100,200), (174,0,174,200), (50,50,100,100) ]
        start_angles = [ 0, 0, 0 ]
        end_angles = [ pi, 2*pi, 2*pi ]
        start_angles[2], end_angles[2] = (-pi/2, pi/2) if ann.ellipsoid.axes[0][2] > -kEps else (pi/2, 3*pi/2)
        
        ## 0 is the cross section
        ## 1 is the major axis (looks like the skeleton)
        ## 2 is the silhouette (only if the cross section isn't tilted, otherwise it's skewed)
        if kDebugging:
            ellipse_axes = (0,1,2)
        else:
            ellipse_axes = (0,2)
        for i in ellipse_axes:
            glLineWidth( widths[i] )
            glColor4ub( *cols[i] )
            draw_ellipse(
                ann.ellipsoid.center,
                ann.ellipsoid.axes[ (i+1) % 3 ],
                ann.ellipsoid.axes[ (i+2) % 3 ],
                start_radians = start_angles[i],
                end_radians = end_angles[i]
                )
        
        
        ## Draw the handles.
        ## Offset them a little towards the camera in z so they appear above the shapes.
        glMatrixMode( GL_MODELVIEW )
        glPushMatrix()
        glTranslatef( 0,0,.01 )
        handles = ann.Handles if draw_handles else []
        for handle in handles:
            fill_alpha = 1 if handle == picked_handle else .5
            
            if handle.kind == ann.Axis0Handle:
                if handle.get_pos()[2] + kEps < handle.annotation.ellipsoid.center[2]: continue
                
                ## Blue handles at either end of the ellipse's major axis.
                glColor4f( 0,0,.5, fill_alpha )
                draw_circle_interior( handle.get_pos(), AnnotateCanvas2.kAxis0HandleRadius )
                glLineWidth( 2 )
                glColor3f( 0,0,.5 )
                draw_circle( handle.get_pos(), AnnotateCanvas2.kAxis0HandleRadius )
            elif handle.kind == ann.Axis1Handle:
                ## Light red handles at either end of the ellipse's minor axis.
                glColor4f( .5,.25,.25, fill_alpha )
                draw_circle_interior( handle.get_pos(), AnnotateCanvas2.kAxis1HandleRadius )
                glLineWidth( 2 )
                glColor3f( .5,.25,.25 )
                draw_circle( handle.get_pos(), AnnotateCanvas2.kAxis1HandleRadius )
            elif handle.kind == ann.CrossSectionHandle:
                ## A handle for the cross section.
                glColor4f( .5,.25,.5, fill_alpha )
                draw_circle_interior( handle.get_pos(), AnnotateCanvas2.kCrossSectionHandleRadius )
                glLineWidth( 2 )
                glColor3f( .5,.25,.5 )
                draw_circle( handle.get_pos(), AnnotateCanvas2.kCrossSectionHandleRadius )
            elif handle.kind == ann.ParentIntersectionHandle:
                ## A handle for the intersection between this annotation and its parent.
                glColor4f( 1,.5,0, fill_alpha )
                draw_circle_interior( handle.get_pos(), AnnotateCanvas2.kParentIntersectionHandle )
                glLineWidth( 2 )
                glColor3f( 1,.5,0 )
                draw_circle( handle.get_pos(), AnnotateCanvas2.kParentIntersectionHandle )
            elif handle.kind == ann.DiscreteZHandle:
                ## A handle for the discrete z depth.
                ## TODO: An arrow which indicates that the annotation at the start is in front of
                ##       the annotation at the back.
                glColor4f( 0,0,0, fill_alpha )
                txt = 'Z: ' + str( ann.z )
                world_size = append( draw_text_world_size( self.world2screen, self.screen2world, txt ), 0. )
                draw_text(
                    asarray( handle.get_pos() ) - .5 * asarray( world_size ),
                    txt
                    )
            elif handle.kind == ann.WholeEllipseHandle:
                pass
            else:
                print('Unknown handle kind:', handle.kind)
        
        if ann.parentage is not None: #and kDebugging:
            from surface_intersections import ellipsoid_ellipsoid_intersection_lines
            line_segments = ellipsoid_ellipsoid_intersection_lines(
                ann.ellipsoid,
                ann.parentage.parent.ellipsoid
                )
            
            ## Since we're intersecting with something approximate, project onto the actual surface.
            ## TODO: move this inside ellipsoid_ellipsoid_intersection_lines()
            line_segments = [
                ( closest_point_on_ellipse_to_point( ( ann.ellipsoid.center, ann.ellipsoid.axes ), seg[0] ),
                closest_point_on_ellipse_to_point( ( ann.ellipsoid.center, ann.ellipsoid.axes ), seg[1] ) )
                for seg in line_segments
                ]
            
            glColor4f( 0, .5, 1., .75 )
            glLineWidth( 2 )
            glBegin( GL_LINES )
            for line in line_segments:
                glVertex3f( *line[0] )
                glVertex3f( *line[1] )
            glEnd()
            
            '''
            ## Shoot a bunch of rays through this annotation and its parent.
            pts = intersection_curve_ellipsoid_ellipsoid(
                ann.ellipsoid,
                ann.parentage.parent.ellipsoid
                )
            glColor4f( 0, .5, 1., .75 )
            #glLineWidth( 2 )
            #draw_lineloop( pts )
            glPointSize( 2 )
            glBegin( GL_POINTS )
            for pt in pts: glVertex2f( *pt )
            glEnd()
            '''
        
        glMatrixMode( GL_MODELVIEW )
        glPopMatrix()
    
    def DrawScene( self, draw_handles = True, callback = None ):
        glClearColor( .9,.9,.9,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        self.prepare2D()
        
        ## Draw the gltex if there is one.
        if self.annotations.gltex is not None:
            glColor4f( 1,1,1, .5 )
            draw_gltexture( self.annotations.gltex, self.annotations.tex_pos, self.annotations.tex_size )
        
        if self.show_annotations:
            glEnable( GL_DEPTH_TEST )
            
            ## Sort by z so we can draw the annotations in depth order (important for transparency).
            self.annotations.sortz()
            
            ## Draw the annotations and their handles, back to front.
            picked_handle = None
            if draw_handles:
                picked_handle = self.PickHandle( self.event2screen( self.GetMousePosition() ) )
                #print 'picked:', picked_handle
            for ann in self.annotations: self.DrawAnnotation( ann, picked_handle = picked_handle, draw_handles = draw_handles )
        
        if callback is not None: callback()
        
        self.SwapBuffers()
    
    def doDraw( self ):
        return self.mode.doDraw()
    
    def doMouseDown( self, evt ):
        return self.mode.doMouseDown( evt )
    
    def doMouseUp( self, evt ):
        return self.mode.doMouseUp( evt )
    
    def doMouseMove( self, evt ):
        return self.mode.doMouseMove( evt )
    
    def doMouseDrag( self, evt ):
        return self.mode.doMouseDrag( evt )
    
    def doKeyDown( self, evt ):
        return self.mode.doKeyDown( evt )
    
    def doKeyUp( self, evt ):
        return self.mode.doKeyUp( evt )
    
    def doContextMenu( self, evt ):
        return self.mode.doContextMenu( evt )


class ResultCanvas( BasicGLCanvas ):
    
    EYE_DISTANCE = 100.
    
    def __init__( self, parent, annotations, pos = wx.DefaultPosition, size = wx.DefaultSize, shared = None ):
        BasicGLCanvas.__init__( self, parent, pos, size, debugging = True, shared = shared )
        
        self.camera = Camera()
        ## Center the square inside the viewport.
        self.camera.proj = Camera.Ortho( *self.world_bounds )
        self.camera.eye[2] = self.EYE_DISTANCE
        self.camera.near_far[1] = 1000.
        
        self.annotations = annotations
        self.image_depth = -.5
        
        debugging = False
        self.rotator = RotateHandler( self.camera, self.Refresh, self.event2world, debugging = debugging )
        forward_events( self, self.rotator )
    
    def GetBestSize( self ):
        return (500,500)
    
    def doDraw( self ):
        glClearColor( .5,.5,.5,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        self.prepare3D()
        
        #self.drawTexture()
        ## Clear the depth buffer.
        ## Actually, this looks weird.
        # glClear( GL_DEPTH_BUFFER_BIT )
        
        for ann in self.annotations:
            if ann.ellipse is not None:
                glDisable( GL_LIGHTING )
                
                glLineWidth( 3 )
                glColor4ub( 174, 174, 0, 200 )
                draw_ellipse( *ann.ellipse )
            
            if ann.ellipsoid is not None:
                glEnable( GL_LIGHTING )
                
                glEnable( GL_COLOR_MATERIAL )
                glColor3ub( 200, 200, 255 )
                ann.ellipsoid.Draw()
            
            if kDebugging:
                if hasattr( ann, '_decs' ):
                    glDisable( GL_LIGHTING )
                    for dec in ann._decs:
                        if dec[0] == 'line':
                            [ None, None, None, glColor3f, glColor4f ][ len( dec[1] ) ]( *dec[1] )
                            glBegin( GL_LINES )
                            for pt in dec[2:]:
                                [ None, None, glVertex2d, glVertex3d, glVertex4d ][ len( pt ) ]( *pt )
                            glEnd()
        
            if kDebugging:
                if ann.ellipsoid is not None:
                    glDisable( GL_LIGHTING )
                    glDisable( GL_DEPTH_TEST )
                    
                    cols = [ (50,50,100,200), (174,0,174,200), (50,50,100,100) ]
                    
                    for i in range(3):
                        glLineWidth( 2 )
                        glColor4ub( *cols[i] )
                        draw_ellipse(
                            ann.ellipsoid.center,
                            ann.ellipsoid.axes[ (i+1) % 3 ],
                            ann.ellipsoid.axes[ (i+2) % 3 ]
                            )
                    
                    glEnable( GL_DEPTH_TEST )
        
        ## Draw the gltex if there is one.
        if self.annotations.gltex is not None:
            self.drawTexture(
                self.annotations.gltex,
                self.annotations.tex_pos,
                self.annotations.tex_size
                )
        
        self.SwapBuffers()
    
    def drawTexture( self, gltex, tex_pos, tex_size ):
        if self.image_depth is None: return
        
        glDisable( GL_LIGHTING )
        
        glMatrixMode( GL_MODELVIEW )
        glPushMatrix()
        glTranslate( 0,0,self.image_depth )
        
        glColor4f( 1,1,1, .5 )
        draw_gltexture( gltex, tex_pos, tex_size )
        
        glPopMatrix()
    
    def prepare2D( self ):
        glDisable( GL_DEPTH_TEST )
        glDisable( GL_LIGHTING )
        
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(-1, 1, -1, 1, -100, 100)
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
    
    def prepare3D( self ):
        glEnable( GL_DEPTH_TEST )
        glEnable( GL_LIGHTING )
        
        self.center_camera_on_data()
        ## Center the square inside the viewport.
        self.camera.proj = Camera.Ortho( *self.world_bounds )
        self.camera.Apply()
        ## TODO: Something with the light.  Right now, we're getting OpenGL's default light.
        
        set_headlight( self.camera )
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE )
        glShadeModel( GL_SMOOTH )
    
    def reset_rotation( self ):
        #self.camera.center = arrayf( ( 0,0,0 ) )
        self.camera.eye = arrayf( ( 0,0,self.EYE_DISTANCE ) )
        self.camera.up = arrayf( ( 0,1,0 ) )
    
    def event2world( self, pos ):
        ## Use this opportunity to make sure the camera is focused on the object's COM.
        self.center_camera_on_data()
        return BasicGLCanvas.event2world( self, pos )
    
    def center_camera_on_data( self ):
        ## By default, the origin.
        self.camera.center = zeros(3)
        ## Keep the up vector pointed in positive y, so we don't get twisted.
        self.camera.up[:] = (0,1,0)
        
        ## UPDATE: I don't want this anymore.
        return
        
        ## Collect all ellipsoid centers.
        centers = [ ann.ellipsoid.center for ann in self.annotations if ann.ellipsoid is not None ]
        centers = asarray( centers )
        
        if len( centers ):
            ## Average the centers of the ellipses.
            #self.camera.center = mean( centers, 0 )
            
            ## Take the middle of the bounding box of the centers.
            self.camera.center = .5 * ( centers.min(0) + centers.max(0) )


class AnnotateFrame( wx.Frame ):
    
    RESET_BUTTON = 10
    ITERATE_BUTTON = 20
    ITERATE_INF_BUTTON = 30
    
    def __init__( self, parent, annotations ):
        wx.Frame.__init__( self, parent, wx.ID_ANY, 'Annotation Modeler' )
        
        self.annotations = annotations
        
        ## Add a toolbar for the modes.  Do this before adding the canvases.
        TBFLAGS = ( wx.TB_HORIZONTAL
            | wx.NO_BORDER
            | wx.TB_FLAT
            #| wx.TB_TEXT
            #| wx.TB_HORZ_LAYOUT
            )
        tb = self.CreateToolBar( TBFLAGS )
        tbsize = ( 64, 64 )
        tb.SetToolBitmapSize( tbsize )
        
        panel = wx.Panel( self, wx.ID_ANY, size=(1000,530) )
        self.sizer = wx.BoxSizer( wx.VERTICAL )
        
        ## Add the side-by-side gl views
        glviews = wx.BoxSizer( wx.HORIZONTAL )
        ## Create the AnnotateCanvas.
        canvas = AnnotateCanvas2( panel, self.annotations )
        canvas.Show()
        self.annotate_canvas = canvas
        canvas.SetSize( canvas.GetBestSize() )
        glviews.Add( canvas, 1, wx.EXPAND )
        ## Create the ResultCanvas.
        canvas = ResultCanvas( panel, self.annotations, shared = self.annotate_canvas.GetContext() )
        self.result_canvas = canvas
        canvas.SetSize( canvas.GetBestSize() )
        glviews.Add( canvas, 1, wx.EXPAND )
        self.sizer.Add( glviews, 1, wx.EXPAND )
        
        ## Add the buttons
        buttons = wx.BoxSizer( wx.HORIZONTAL )
        ## Reset button
        btn = wx.Button( panel, self.RESET_BUTTON, "Reset" )
        self.Bind( wx.EVT_BUTTON, self.OnClick, btn )
        btn.SetDefault()
        btn.SetSize( btn.GetBestSize() )
        buttons.Add( btn, 0, wx.ALIGN_CENTER )
        
        ## Iterate button
        #btn = wx.Button( panel, self.ITERATE_BUTTON, "Iterate Once" )
        #self.Bind( wx.EVT_BUTTON, self.OnClick, btn )
        #btn.SetDefault()
        #btn.SetSize( btn.GetBestSize() )
        #buttons.Add( btn, 0, wx.ALIGN_CENTER )
        ## Iterate button
        #btn = wx.Button( panel, self.ITERATE_INF_BUTTON, "Iterate inf" )
        #self.Bind( wx.EVT_BUTTON, self.OnClick, btn )
        #btn.SetDefault()
        #btn.SetSize( btn.GetBestSize() )
        #buttons.Add( btn, 0, wx.ALIGN_CENTER )
        
        ## Add the buttons sizer to the parent sizer
        self.sizer.Add( buttons, 0, wx.ALIGN_BOTTOM )
        
        panel.SetSizer( self.sizer )
        self.sizer.Fit( panel )
        self.Fit()
        
        import os
        self.last_wd = os.getcwd()
        self.save_filename = None
        
        menubar = wx.MenuBar()
        
        filem = wx.Menu()
        item = filem.Append( wx.ID_ANY, "Open Scene...\tctrl-o" )
        self.Bind( wx.EVT_MENU, self.OnFileMenuOpenScene, item )
        item = filem.Append( wx.ID_ANY, "Save Scene...\tctrl-s" )
        self.Bind( wx.EVT_MENU, self.OnFileMenuSaveScene, item )
        item = filem.Append( wx.ID_ANY, "Save Scene As...\tctrl-shift-s" )
        self.Bind( wx.EVT_MENU, self.OnFileMenuSaveSceneAs, item )
        filem.AppendSeparator()
        item = filem.Append( wx.ID_ANY, "Load Background PNG...\tctrl-l" )
        self.Bind( wx.EVT_MENU, self.OnFileMenuLoadBackground, item )
        filem.AppendSeparator()
        item = filem.Append( wx.ID_ANY, "Save Annotation Screenshot..." ) #\tctrl-3" )
        self.Bind( wx.EVT_MENU, lambda evt: self.SaveScreenShot( self.annotate_canvas, 'annotate canvas' ), item )
        item = filem.Append( wx.ID_ANY, "Save Result Screenshot..." ) #\tctrl-4" )
        self.Bind( wx.EVT_MENU, lambda evt: self.SaveScreenShot( self.result_canvas, 'result canvas' ), item )
        menubar.Append( filem, "&File" )
        
        modem = wx.Menu()
        mode_ids = AnnotateCanvas2.GetAllModeIDs()
        for mode_id in mode_ids:
            ## TODO: make these into a check item group and have only one checked.
            item = modem.Append( mode_id, AnnotateCanvas2.GetLongNameForModeID( mode_id ) + "\tctrl-%s" % (mode_id,) )
            self.Bind( wx.EVT_MENU, lambda evt: self.annotate_canvas.SetModeByID( evt.Id ), item )
        menubar.Append( modem, "&Mode" )
        
        viewm = wx.Menu()
        item = viewm.Append( wx.ID_ANY, "Reset 3D Result View\tctrl-r" )
        self.Bind( wx.EVT_MENU, self.OnReset3DResultView, item )
        item = viewm.AppendCheckItem( wx.ID_ANY, "Result View Image Behind\tctrl-b" )
        if self.result_canvas.image_depth < 0: item.Check()
        self.Bind( wx.EVT_MENU, self.OnToggleResultViewImageDepth, item )
        viewm.AppendSeparator()
        item = viewm.AppendCheckItem( wx.ID_ANY, "Show Annotations\tctrl-a" )
        if self.annotate_canvas.show_annotations: item.Check()
        self.Bind( wx.EVT_MENU, self.OnShowAnnotations, item )
        menubar.Append( viewm, "&View" )
        
        testm = wx.Menu()
        item = testm.AppendCheckItem( wx.ID_ANY, "Test new overlap method\tctrl-t" )
        if kTestNewOverlap: item.Check()
        self.Bind( wx.EVT_MENU, self.OnToggleTestNewOverlap, item )
        menubar.Append( testm, "&Test" )
        
        self.SetMenuBar( menubar )
        
        
        ## Add modes to the toolbar
        ## TODO: make these into a check item group and have only one checked.
        mode_ids = AnnotateCanvas2.GetAllModeIDs()
        for mode_id in mode_ids:
            short_name = AnnotateCanvas2.GetShortNameForModeID( mode_id )
            long_name = AnnotateCanvas2.GetLongNameForModeID( mode_id )
            
            img = wx.Image( os.path.join( ".", "modes", short_name.lower() + ".png" ), wx.BITMAP_TYPE_ANY )
            img = img.Scale( tbsize[0], tbsize[1], quality = wx.IMAGE_QUALITY_HIGH )
            bitmap = img.ConvertToBitmap()
            tb.AddRadioLabelTool( mode_id, short_name, bitmap, shortHelp=short_name, longHelp=long_name )
            self.Bind( wx.EVT_TOOL, lambda evt: self.annotate_canvas.SetModeByID( evt.Id ), id=mode_id )
        
        tb.Realize()
        
        ## Do I need to call this again?
        self.Fit()
    
    def SaveScreenShot( self, canvas, suggested_name ):
        dialog = wx.FileDialog(
            self, "Save screenshot...", self.last_wd,
            suggested_name + '.png',
            "PNG|*.png|JPEG|*.jpg",
            wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
            )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            
            self.save_filename = str( dialog.GetPath() )
            canvas.SaveToFile( self.save_filename )
    
    def OnToggleResultViewImageDepth( self, evt ):
        self.result_canvas.image_depth *= -1.
        self.result_canvas.Refresh()
    
    def OnToggleTestNewOverlap( self, evt ):
        global kTestNewOverlap
        kTestNewOverlap = not kTestNewOverlap
    
    def OnReset3DResultView( self, evt ):
        self.result_canvas.reset_rotation()
        self.result_canvas.Refresh()
    
    def OnShowAnnotations( self, evt ):
        self.annotate_canvas.show_annotations = not self.annotate_canvas.show_annotations
        self.annotate_canvas.Refresh()
    
    def OnFileMenuOpenScene( self, evt ):
        dialog = wx.FileDialog( self, "Open scene file...", self.last_wd, '', "Annotation Scene|*.annscene", style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            
            save_filename = str( dialog.GetPath() )
            try:
                self.annotations.Load( save_filename )
            except:
                raise
            else:
                self.save_filename = save_filename
            finally:
                self.GetTopLevelParent().Refresh()
    
    def OnFileMenuLoadBackground( self, evt ):
        dialog = wx.FileDialog( self, "Load background PNG...", self.last_wd, '', "PNG|*.png", style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            
            self.annotations.SetImage( str( dialog.GetPath() ) )
            self.GetTopLevelParent().Refresh()
    
    def OnFileMenuSaveScene( self, evt ):
        if self.save_filename is not None:
            self.annotations.Save( self.save_filename, overwrite = True )
            return
        
        suggested_name = 'annotated scene'
        import os
        if self.annotations.image_path: suggested_name = os.path.split( os.path.splitext( self.annotations.image_path )[0] )[1]
        
        dialog = wx.FileDialog(
            self, "Save scene file...", self.last_wd,
            suggested_name + '.annscene',
            "Annotation Scene|*.annscene",
            wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
            )
        if dialog.ShowModal() == wx.ID_OK:
            self.last_wd = os.path.split( dialog.GetPath() )[0]
            
            self.save_filename = str( dialog.GetPath() )
            self.annotations.Save( self.save_filename, overwrite = True )
    
    def OnFileMenuSaveSceneAs( self, evt ):
        self.save_filename = None
        self.OnFileMenuSaveScene( evt )
    
    def OnClick( self, event ):
        if event.GetId() == self.RESET_BUTTON:
            print("reset!")
            self.annotations.Reset( keep_texture = True )
            ## Refresh all
            self.GetTopLevelParent().Refresh()
        
        elif event.GetId() == self.ITERATE_BUTTON:
            print("iterate!")
            self.IterateOnce()
        
        elif event.GetId() == self.ITERATE_INF_BUTTON:
            if hasattr( self, 'iterate_inf_timer_id' ):
                Timers.CancelTimer( self.iterate_inf_timer_id )
                del self.iterate_inf_timer_id
            else:
                self.iterate_inf_timer_id = Timers.RequestTimer( self.IterateOnce, 10, once = False )
    
    def IterateOnce( self ):
        if len( self.annotations ) > 0 and self.annotations[-1].Iterate( 1 ):
            ## Refresh all
            self.GetTopLevelParent().Refresh()

def parse_options( argv ):
    global kDebugging
    
    try:
        argv.remove( '-d' )
        kDebugging = True
    except: pass
    try:
        argv.remove( '--debug' )
        kDebugging = True
    except: pass

def main():
    import sys
    argv = list( sys.argv[1:] )
    parse_options( argv )
    
    import wx
    # app = wx.App() ## <-- on OS X tracebacks and stdout/stderr appear in a dialog
    app = wx.PySimpleApp()
    
    ## load the PyCrust shell
    ## NOTE: This generates the following errors and is very slow switching between the app window
    ##       and the PyCrust window.
    ## Sun Jun 15 20:19:36 yotronagon.local Python[32825] <Error>: CGContextRestoreGState: invalid context
    if kDebugging:
        from wxBasicGLCanvas import load_PyCrust
        crust = load_PyCrust( app )
    
    annotations = Annotations()
    
    frame = AnnotateFrame( None, annotations )
    frame.Show(True)
    
    for fname in argv:
        if fname.lower().endswith( '.annscene' ):
            annotations.Load( fname )
            frame.Refresh()
        elif fname.lower().endswith( '.png' ):
            annotations.SetImage( fname )
            frame.Refresh()
    
    ## So I can easily access the annotations from 'app' in PyCrust
    app.annotations = annotations
    app.MainLoop()
    app.Destroy()

if __name__ == '__main__': main()
