#!/usr/bin/env python3

#!/usr/bin/arch -i386 /usr/bin/python2.6
## I have all of my python modules installed for only one specific
## python version which is now explicitly asked-for with the #!.
##!/usr/bin/env python

from wx_annotate_shell import *
## I have to import wx_annotate_shell and refer to wx_annotate_shell.kDebugging
## because the imported kDebugging symbol is not a reference to the dict entry
## with the same name in wx_annotate_shell.
import wx_annotate_shell

import wx_annotate_symmetry_xsection
import wx_annotate_sphere
import annotations_ellipsoid
import annotations_sweep
import annotations_symbols
from wx_annotate_symmetry_xsection import SymmetryXSectionShape
from annotations_sweep import SweepShape
from annotations_rotund import Connection
from annotations_glwidgets import draw_widget

from wx_annotate_symmetry_xsection import DrawDebugging, PlaceInDepth

## kHackDrawEdges = False
# from mesh_gui import draw_mesh_edges

from helpers import *

kRealTimeHacks = True

## The delete key on my MacBook Pro sends key code 8, which doesn't correspond to anything
## in wxPython.
kDeleteKeys = (wx.WXK_DELETE,8)

## When picking, ignore edges farther than this distance (in world coordinates).
kEdgeDistanceThreshold = .03
def EdgeDistanceThreshold( canvas ):
    return kEdgeDistanceThreshold / canvas.GetZoom()
## When picking, ignore points farther than this distance (in world coordinates).
#kPointDistanceThreshold = .03
## For the paper:
kPointDistanceThreshold = .04
def PointDistanceThreshold( canvas ):
    return kPointDistanceThreshold / canvas.GetZoom()

## This is the default world-space line width.
#kWorldDefaultLineWidth = .01
## For the paper:
kWorldDefaultLineWidth = .02
def WorldDefaultLineWidth( canvas ):
    return kWorldDefaultLineWidth / canvas.GetZoom()

try:
    import FTGL
    font = FTGL.TextureFont( 'OpenSans-ExtraBold.ttf' )
    font.FaceSize( 24, 72 )
except Exception:
    class FakeFont( object ):
        def __init__( self ):
            self.descender = 0.
        def Render( self, text ): pass
        def BBox( self, text ): return ( 0, 0, 0, 100*len( text ), 100, 10 )
    font = FakeFont()

def draw_mode_label_text( text ):
    glColor4ub( 0,0,0, 200 )
    ## Inset by the descender horizontally, too.  It's a nice padding.
    origin = (font.descender, font.descender)
    
    glMatrixMode( GL_PROJECTION )
    glPushMatrix()
    glLoadIdentity()
    
    ## Like parameters to glOrtho: xmin, xmax, ymin, ymax
    world_bounds = [ -1, 1, -1, 1 ]
    width, height = glGetIntegerv( GL_VIEWPORT )[2:]
    if width >= height:
        aspect = width / float( height )
        world_bounds[ 0 ] = -aspect
        world_bounds[ 1 ] = aspect
    else:
        aspect = height / float( width )
        world_bounds[ 2 ] = -aspect
        world_bounds[ 3 ] = aspect
    glOrtho(*( world_bounds + [ -100, 100 ] ))
    
    glMatrixMode( GL_MODELVIEW )
    glPushMatrix()
    glLoadIdentity()
    
    glTranslatef( world_bounds[0], world_bounds[2], 0 )
    scale = .4 * 1./72.
    glScalef( scale, scale, scale )
    glTranslatef( -origin[0], -origin[1], 0. )
    
    glEnable( GL_TEXTURE_2D )
    font.Render( text )
    glDisable( GL_TEXTURE_2D )
    
    glMatrixMode( GL_MODELVIEW )
    glPopMatrix()
    glMatrixMode( GL_PROJECTION )
    glPopMatrix()


## These must be integers, because glColor3ub() complains if it's passed floating point sequences.
Color = {
    'mesh': ( 200, 200, 255 ),
    #'silhouette': ( 234, 145, 248 ),
    'silhouette': ( 0, 0, 0 ),
    'symmetry_xsection': ( 174//2, 237//2, 111//2 ),
    'circular_xsection': ( 213, 99, 191 ),
    'diameter': ( 213//2, 99//2, 191//2 ),
    'connection': ( 174//2, 237//2, 111//2 ),
    'spine': ( 0, 0, 0 ),
    'symmetry_sheet': ( 103, 83, 99 ),
    'orange': ( 232, 146, 26 ),
    'shape_fill': ( 255,255,255, int(round(.9*255)) )
    }
def ub2f( color ):
    '''
    Given an unsigned-byte color values (sequence of numbers from 0 to 255),
    returns a sequence of floating point color values (numbers from 0 to 1).
    '''
    return asarray( color, dtype = float ) / 255.

class StrokeDrawer( ModeProxy ):
    def __init__( self, canvas, stroke, auto_unforward_events = True ):
        ModeProxy.__init__( self, canvas )
        
        assert len( stroke ) == 0
        
        self.stroke = stroke
        
        self.auto_unforward_events = auto_unforward_events
    
    def event2stroke( self, event_pos ):
        return self.event2world( event_pos )
    
    def doMouseDown( self, evt ):
        del self.stroke[:]
        self.stroke.append( self.event2stroke( evt.GetPosition() ) )
        
        self.doPointAdded()
    
    def doMouseDrag( self, evt ):
        pt = self.event2stroke( evt.GetPosition() )
        
        if pt != self.stroke[-1]:
            self.stroke.append( pt )
            #if wx_annotate_shell.kDebugging: print pt
            
            self.doPointAdded()
    
    def doMouseUp( self, evt ):
        self.doMouseDrag( evt )
        
        self.doStrokeFinished()
        
        if self.auto_unforward_events:
            self.stroke = None
            unforward_events( self )
        else:
            del self.stroke[:]
    
    def doDraw( self ):
        if self.stroke is None: return
        
        self.doSetGLLineStyle()
        #draw_linestrip( self.stroke )
        draw_thick_linestrip2d( self.stroke, self.doGetLineWidthWorld() )
    
    def doPointAdded( self ): pass
    def doStrokeFinished( self ): pass
    def doSetGLLineStyle( self ):
        glColor3ub( 0,0,0 )
        self.glLineWidthWorld( self.doGetLineWidthWorld() )
    def doGetLineWidthWorld( self ): return 2*WorldDefaultLineWidth( self )


def PickShapeSurfaces( annotations, world_pos, world_dir = (0,0,-1) ):
    '''
    Given a sequence 'annotations' of objects with a '.ray_surface_intersection_ts()' method,
    a 2+d world-space position 'world_pos',
    and optional 3d world-space vector 'world_dir' (defaults to -z),
    returns a list of unsorted (t, annotation) tuples
    corresponding to intersections of the world-ray with annotations.
    The .ray_surface_intersection_ts() method should take a 3d ray as the first parameter
    and optional 'ray_is_line' parameter as the second, and return a sequence of t's along the
    ray which intersect the shape.  The optional 'ray_is_line' parameter determines whether
    negative t's should be returned.
    
    Note that each annotation will appear at most once, with the smallest 't'.
    Note that t's may be negative.
    Note that the returned list is not sorted in 't', but is stable in 'annotations', meaning
    that the order of returned elements will match their order in 'annotations.'
    '''
    
    ray = ( zerosf( 3 ), asarrayf( world_dir ) )
    ray[0][:2] = world_pos[:2]
    
    t_anns = []
    ## Iterate over all annotations with a ray_surface_intersection_ts() method.
    for shape in [ ann for ann in annotations if hasattr( ann, 'ray_surface_intersection_ts' ) ]:
        ts = None
        try:
            ts = shape.ray_surface_intersection_ts( ray, ray_is_line = True )
        ## We get a RuntimeError if no mesh has been generated yet in some shapes; ignore them for now.
        except RuntimeError:
            print("PickShapeSurfaces() got a RuntimeError because the mesh hasn't been generated yet.")
        
        if ts is not None and len( ts ) > 0:
            t_anns.append( ( min(ts), shape ) )
    
    return t_anns

def ConnectionBetween( anns, shape1, shape2 ):
    '''
    Given a list_type_queryable of annotations 'anns' and two shapes,
    returns the Connection between them in 'anns' or None, if no such Connection exists.
    
    used
    '''
    
    cons = ConnectionsBetween( anns, shape1, shape2 )
    assert len( cons ) <= 1
    return cons[0] if len( cons ) > 0 else None

def ConnectionsBetween( anns, shape1, shape2 ):
    '''
    Given a list_type_queryable of annotations 'anns' and two shapes,
    returns the list of Connections between them in 'anns'.
    
    used
    '''
    
    result = []
    ## A shape isn't connected to itself.
    if shape1 == shape2: return result
    
    for con in anns.of_type( Connection ):
        if shape1 in con.regions and shape2 in con.regions:
            result.append( con )
    
    return result

def RemoveAnnotation( anns, ann_to_remove ):
    '''
    Given a list of annotations 'anns' and a specific annotation to remove 'ann_to_remove',
    removes 'ann_to_remove' from 'anns' and performs any necessary cleanup,
    such as also removing meta-annotations referencing by 'ann_to_remove'.
    '''
    
    anns.remove( ann_to_remove )
    GarbageCollectMetaAnnotations( anns )

def GarbageCollectMetaAnnotations( anns ):
    '''
    Annotations which reference other annotations, such as connection curves
    and parallel plane constraints, should be deleted when the annotations they reference
    are removed from the list of annotations.
    Given a list of annotations 'anns', this routine garbage collects such annotations.
    '''
    
    remove_me = []
    
    ## Garbage collect connections.
    for connection in anns.of_type( Connection ):
        if False in [ region in anns for region in connection.regions ]:
            remove_me.append( connection )
    
    ## Garbage collect tilt annotations.
    for same_tilt in anns.of_type( annotations_symbols.SameTilt ):
        same_tilt.garbage_collect( anns )
        if same_tilt.is_garbage( anns ):
            remove_me.append( same_tilt )
    
    ## Garbage collect size annotations.
    for same_size in anns.of_type( annotations_symbols.SameSize ):
        same_size.garbage_collect( anns )
        if same_size.is_garbage( anns ):
            remove_me.append( same_size )
    
    ## Garbage collect radius annotations.
    for same_radius in anns.of_type( annotations_symbols.SameRadius ):
        same_radius.garbage_collect( anns )
        if same_radius.is_garbage( anns ):
            remove_me.append( same_radius )
    
    ## Garbage collect mirror across annotations.
    for mirror_ann in anns.of_type( annotations_symbols.MirrorAcrossParent ):
        mirror_ann.garbage_collect( anns )
        if mirror_ann.is_garbage( anns ):
            remove_me.append( mirror_ann )
    
    ## Garbage collect mirror align annotations.
    for mirror_ann in anns.of_type( annotations_symbols.MirrorAlignAcrossParent ):
        mirror_ann.garbage_collect( anns )
        if mirror_ann.is_garbage( anns ):
            remove_me.append( mirror_ann )
    
    for ann in remove_me: anns.remove( ann )

def ApplyGeometryConstraints( annotations ):
    '''
    Modifies sequence of annotations 'anns' by applying all constraints (placement in depth, same
    tilt angle, etc).
    '''
    
    #print '+++ ApplyGeometryConstraints()'
    ## Force generation of mesh vertex normals, which will be copied to anns_beautified.
    ## NOTE: This doesn't speed things up if beautification does anything, because the meshes
    ##       will all be regenerated.
    #for ann in annotations.anns:
    #    if hasattr( ann, 'mesh' ): ann.mesh.vertex_normals
    from copy import deepcopy
    #print '+++ deepcopy( annotations.anns )'
    annotations.anns_beautified = deepcopy( annotations.anns )
    #print '--- deepcopy( annotations.anns )'
    annotations_symbols.ApplyConstraints( annotations.anns_beautified )
    ## Moved inside of annotations_symbols.ApplyConstraints() so that it can happen before mirroring.
    #PlaceInDepth( annotations.anns_beautified )
    #print '--- ApplyGeometryConstraints()'

class AnnotationsBeautify( Annotations ):
    def Beautify( self ):
        ApplyGeometryConstraints( self )

def GetCurrentShapeState( anns, shape ):
    '''
    Finds the shape in the list of annotations 'anns' with the same uuid as 'shape'.
    If the user performs an undo or redo, 'shape' will no longer be in the main list
    of annotations; this function will find it.
    
    If no shape with the same uuid exists in 'anns', returns None.  This can happen if the shape
    is created or destroyed.
    '''
    if shape not in anns:
        same = [
            candidate for candidate in anns.of_type( type( shape ) )
            if shape.uuid == candidate.uuid
            ]
        if len( same ) == 0: return None
        assert len( same ) == 1
        shape = same[0]
    return shape

def neighbors( anns, shape ):
    '''
    Given a type-list of annotations 'anns' and a shape 'region',
    returns all shapes that share a Connection with 'region',
    '''
    
    result = set()
    for con in anns.of_type( Connection ):
        if shape in con.regions:
            result.update( con.regions )
    result.remove( shape )
    return result

def run_dfs( connections, pre_visit = None, post_visit = None, roots = None ):
    '''
    Given a sequence of Connection objects 'connections',
    run a depth-first search on all regions,
    calling pre_visit( incoming_connection, region ) and post_visit( incoming_connection, region ) for every region in every
    connection.region.
    If roots is not None, the search will begin on each region in 'roots'; otherwise, it will
    be called on all regions in all Connections until they have all been visited.
    If 'pre_visit' is None, it will not be called.
    If 'post_visit' is None, it will not be called.
    
    untested
    '''
    
    visited_regions = set()
    
    connected_regions = set()
    for connection in connections:
            connected_regions.update( connection.regions )
    
    ## A map from a region to neighboring (connected) regions.
    region2connection_neighbors = {}
    ## Initialize all regions' neighbors to empty sets.
    for region in connected_regions: region2connection_neighbors.setdefault( region, set() )
    if roots is not None:
        for region in roots: region2connection_neighbors.setdefault( region, set() )
    ## Fill all regions' neighbor sets.
    for connection in connections:
        for region in connection.regions:
            region2connection_neighbors[ region ].update( [ ( connection, r ) for r in connection.regions ] )
    ## Right now, regions may be their own neighbors.
    ## If so, remove them from their own neighbor sets.
    for region, connection_neighbors in region2connection_neighbors.items():
        to_remove = [ cn for cn in connection_neighbors if cn[1] is region ]
        for cn in to_remove: connection_neighbors.remove( cn )
    
    def dfs( incoming_connection, region ):
        if region in visited_regions: return
        
        visited_regions.add( region )
        
        ## pre-visit region
        if pre_visit is not None: pre_visit( incoming_connection, region )
        
        for connection, child in region2connection_neighbors[ region ]:
            dfs( connection, child )
        
        ## post-visit region
        if post_visit is not None: post_visit( incoming_connection, region )
    
    if roots is not None:
        for root in roots:
            dfs( None, root )
    else:
        for root in connected_regions:
            dfs( None, root )

#===================
class EllipseAnnotationCreateMode( ModeProxy ):
    class EllipseSilhouetteDrawer( StrokeDrawer ):
        def __init__( self, canvas ):
            ## Create a new annotation
            self.shape = annotations_ellipsoid.AnnotationData()
            self.shape.stroke = []
            StrokeDrawer.__init__( self, canvas, self.shape.stroke )
        
        def doPointAdded( self ):
            self.Refresh()
        
        def doStrokeFinished( self ):
            self.shape.SetCircle( smooth_open_stroke( self.shape.stroke ) )
            self.shape.stroke = None
            
            self.Refresh()
            
            kTooNarrowAxesRatio = 50
            ## This happens if we couldn't fit an ellipse to the stroke.
            if self.shape.ellipse is None:
                #RemoveAnnotation( self.annotations.anns, self.shape )
                pass
            elif max( list(map( mag, self.shape.ellipse[1:] )) ) / min( list(map( mag, self.shape.ellipse[1:] )) ) > kTooNarrowAxesRatio:
                print('Ellipse too narrow!  Aborting.')
                #RemoveAnnotation( self.annotations.anns, self.shape )
                pass
            else:
                ## Calling this prevents misbehavior the first time an axis0 handle is dragged.
                ## TODO Q: Should I move this inside SetCircle()?
                self.shape.SetCrossSectionFromStroke( 1, True, [ self.shape.ellipsoid.center[:2] + self.shape.ellipsoid.axes[2][:2] ], (0,0,-1) )
                self.SaveStateForUndo()
                self.annotations.anns.append( self.shape )
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
                
                ## Switch to arrow mode.
                #self.SetModeByID( self.GetModeIDForShortName( "Select" ) )
                self.SwitchToolbarTools( self.GetModeIDForShortName( "Select" ) )
                self.mode.selection.add( self.shape )
            
            self.shape = None
        
        def doSetGLLineStyle( self ):
            glColor4ub( *( Color[ 'silhouette' ] + (200,) ) )
            self.glLineWidthWorld( self.doGetLineWidthWorld() )
    
    
    class HandleDragger( ModeProxy ):
        def __init__( self, canvas, handle ):
            ModeProxy.__init__( self, canvas )
            self.SaveStateForUndo()
            
            self.stroke = []
            self.handle = handle
        
        def doMouseDown( self, evt ):
            self.stroke = [ tuple( self.event2world( evt.GetPosition() ) ) ]
            if wx_annotate_shell.kDebugging: print(self.stroke[-1])
        
        def doMouseDrag( self, evt ):
            pt = tuple( self.event2world( evt.GetPosition() ) )
            if pt != self.stroke[-1]:
                self.stroke.append( pt )
                #if wx_annotate_shell.kDebugging: print self.stroke[-1]
                
                #print "del:", asarray( self.stroke[-1] ) - asarray( self.stroke[-2] )
                self.handle.translate( asarray( self.stroke[-1] ) - asarray( self.stroke[-2] ) )
                ## Q: Can we handle placing in depth in real-time, or is it too slow?
                ## A: No, Beautify() is too slow.
                #self.annotations.Beautify()
                #self.GetTopLevelParent().Refresh()
                self.Refresh()
        
        def doMouseUp( self, evt ):
            self.doMouseDrag( evt )
            
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            unforward_events( self )
        
        def doDraw( self ):
            self.DrawEllipseHandles( [self.handle.annotation] )
    
    def PickHandle( self, world_pos ):
        picked_handle = wx_annotate_sphere.PickHandle(
            self.annotations.anns.of_type( annotations_ellipsoid.AnnotationData ),
            world_pos,
            lambda kind: PointDistanceThreshold( self )
            )
        return picked_handle
    
    def doMouseDown( self, evt ):
        ## If we pick a handle, manipulate it.  Otherwise, create a new ellipse annotation.
        
        picked_handle = self.PickHandle( self.event2world( evt.GetPosition() ) )
        if picked_handle:
            print("Picked handle:", picked_handle.kind)
            forward_events( self, self.HandleDragger( self, picked_handle ) )
            self.doMouseDown( evt )
        else:
            forward_events( self, self.EllipseSilhouetteDrawer( self ) )
            self.doMouseDown( evt )
    
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        if keycode in kDeleteKeys:
            
            ## evt.GetPosition() returns garbage.
            world_pos = self.event2world( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            ## If the mouse is not over any annotations, check if we're over a handle.
            ann = None
            
            ## TODO: Convenience function for this annotation picking.
            #ann = self.PickAnnotationOfType( screen_pos, annotations_ellipsoid.AnnotationData )
            t_anns = wx_annotate_sphere.PickTAnnotations(
                self.annotations.anns.of_type( annotations_ellipsoid.AnnotationData ),
                world_pos
                )
            t_anns.sort()
            if len( t_anns ) > 0:
                ann = t_anns[0][1]
            
            if ann is not None:
                self.SaveStateForUndo()
                print('deleting an annotation')
                RemoveAnnotation( self.annotations.anns, ann )
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
        
        else:
            evt.Skip()
    
    
    def doMouseMove( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doDraw( self ):
        def draw_child():
            child = event_child( self )
            if child is not None:
                child.doDraw()
                return True
        
        def draw_callback():
            if draw_child(): return
            
            self.prepare2D()
            
            handle = self.PickHandle( self.event2world( self.GetMousePosition() ) )
            if handle:
                draw_handles = [handle.annotation]
                ## XXX TODO: This is ugly.
                #if wx_annotate_shell.kPaperFigureMode: draw_handles = None
                self.DrawEllipseHandles( draw_handles )
                
                if handle.kind == annotations_ellipsoid.AnnotationData.WholeEllipseHandle:
                    draw_mode_label_text( 'Move' )
                ## TODO: Either don't draw the unpickable, obscured axis 0 (long axis) handle,
                ##       or make it pickable.
                elif handle.kind == annotations_ellipsoid.AnnotationData.Axis0Handle:
                    draw_mode_label_text( 'Rotate/Scale' )
                elif handle.kind == annotations_ellipsoid.AnnotationData.Axis1Handle:
                    draw_mode_label_text( 'Change circle radius' )
                elif handle.kind == annotations_ellipsoid.AnnotationData.CrossSectionHandle:
                    draw_mode_label_text( 'Tilt circular cross section' )
            else:
                draw_mode_label_text( 'Draw ellipse' )
            
        
        self.DrawScene( callback = draw_callback )
#===================

#===================
class SymmetryXSectionPart(object):
    ## SilhouetteCurveInterior is the whole shape
    SilhouetteCurveInterior = 0
    SymmetryCurveInterior = 1
    SilhouetteCurve = 2
    SymmetryCurve = 3
    CircularCrossSection = 4
def PickSymmetryXSectionHandle( shape, world_pos, world_dir = (0,0,-1) ):
    '''
    Given a SymmetryXSection object 'shape' and a world-space position 'world_pos' and optional
    world-space direction 'world_dir' (defaults to -z),
    returns the tuple of one or more elements ( picked SymmetryXSectionPart, per-part info... )
    corresponding to the picked part, or None if no part was picked.
    '''
    world_ray = ( zerosf( 3 ), asarrayf( world_dir ) )
    world_ray[0][:2] = world_pos
    
    ## For each circular cross section, are we inside the ellipse?
    for i, (center, radius, normal) in enumerate( shape.circular_xsections_3d ):
        #print "shape.circular_xsections_3d[i]:", shape.circular_xsections_3d[i]
        plane_isect = ray_plane_intersection( world_ray, ( center, normal ), True )
        if plane_isect is not None and mag2( asarray( center ) - asarray( plane_isect ) ) <= radius**2:
            return ( SymmetryXSectionPart.CircularCrossSection, i )
    
    ## Are we close to the symmetry curve?
    if shape.symmetry_xsection:
        distsqr_pt_edge = closest_distsqr_and_point_and_edge_index_on_line_loop_to_point(
            shape.symmetry_xsection, world_pos[:2]
            )
        if distsqr_pt_edge[0] < EdgeDistanceThreshold( self )**2:
            return ( SymmetryXSectionPart.SymmetryCurve, )
    
    ## Are we close to the silhouette?
    if shape.silhouette:
        distsqr_pt_edge = closest_distsqr_and_point_and_edge_index_on_line_loop_to_point(
            shape.silhouette, world_pos[:2]
            )
        if distsqr_pt_edge[0] < EdgeDistanceThreshold( self )**2:
            return ( SymmetryXSectionPart.SilhouetteCurve, )
    
    ## Are we inside the symmetry curve contour?
    if shape.symmetry_xsection is not None and contour_contains( shape.symmetry_xsection, world_ray[0][:2] ):
        return ( SymmetryXSectionPart.SymmetryCurveInterior, )
    
    ## Are we inside the silhouette contour?
    if shape.silhouette is not None and contour_contains( shape.silhouette, world_ray[0][:2] ):
        return ( SymmetryXSectionPart.SilhouetteCurveInterior, )
    
    return None

def PickClosestSymmetryXSectionHandle( annotations, world_pos, world_dir = (0,0,-1) ):
    '''
    Similar to PickSymmetryXSectionHandle(), except takes a sequence of SymmetryXSection objects
    'annotations' and returns the tuple ( SymmetryXSection object, handle ) of the most specific
    picked handle over all SymmetryXSection objects in 'annotations'.
    '''
    
    handles = [ ( PickSymmetryXSectionHandle( shape, world_pos, world_dir ), shape ) for shape in annotations ]
    handles = [ handle for handle in handles if handle[0] is not None ]
    
    if len( handles ) == 0:
        return None
    else:
        ## Handles have enumerated values in increasing order of specificity, so return the
        ## most-specific (largest numbered) handle.
        closest_handle, closest_shape = max( handles )
        return closest_shape, closest_handle

def PickSymmetryXSections( annotations, world_pos, world_dir = (0,0,-1) ):
    '''
    Given a sequence of SymmetryXSection 'annotations', a 2+d world-space position 'world_pos',
    and optional world-space vector 'world_dir' (defaults to -z),
    returns a list of unsorted (t, annotation) tuples
    corresponding to intersections of the world-ray with annotations.
    Note that each annotation will appear at most once, with the smallest 't'.
    Note that t's may be negative.
    '''
    
    ray = ( zerosf( 3 ), asarrayf( world_dir ) )
    ray[0][:2] = world_pos
    
    t_anns = []
    for shape in annotations:
        try:
            ts = shape.ray_surface_intersection_ts( ray, ray_is_line = True )
        except RuntimeError:
            ts = None
        
        if ts is not None and len( ts ) > 0:
            t_anns.append( ( min(ts), shape ) )
        elif shape.silhouette is not None and contour_contains( shape.silhouette, ray[0][:2] ):
                t_anns.append( ( 0., shape ) )
    
    return t_anns

class SymmetryXSectionAnnotationCreateMode( ModeProxy ):
    class SilhouetteDrawer( StrokeDrawer ):
        def __init__( self, canvas ):
            ## Create a new annotation.
            self.shape = wx_annotate_symmetry_xsection.SymmetryXSectionShape()
            #self.annotations.anns.append( self.shape )
            
            self.shape.silhouette = []
            
            StrokeDrawer.__init__( self, canvas, self.shape.silhouette )
        
        def doPointAdded( self ):
            self.Refresh()
        
        def doStrokeFinished( self ):
            if len( self.shape.silhouette ) >= 3:
                self.shape.silhouette = smooth_closed_stroke( self.shape.silhouette )
                self.shape.update_3d()
                self.SaveStateForUndo()
                self.annotations.anns.append( self.shape )
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            self.shape = None
            self.Refresh()
        
        def doSetGLLineStyle( self ):
            glColor4ub( *( Color[ 'silhouette' ] + (200,) ) )
            self.glLineWidthWorld( self.doGetLineWidthWorld() )
    
    class SymmetryXSectionDrawer( StrokeDrawer ):
        def __init__( self, canvas, shape ):
            canvas.SaveStateForUndo()
            
            self.shape = shape
            self.shape.symmetry_xsection = []
            
            ## Don't pass self.shape.symmetry_xsection, because we don't want the drawing
            ## code to draw the unfinished loop.
            StrokeDrawer.__init__( self, canvas, [] )
        
        def doPointAdded( self ):
            self.Refresh()
        
        def doStrokeFinished( self ):
            if len( self.stroke ) >= 3:
                self.shape.symmetry_xsection = smooth_closed_stroke( self.stroke )
                self.shape.update_3d()
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            self.shape = None
            self.Refresh()
        
        def doSetGLLineStyle( self ):
            glColor4ub( *( Color[ 'symmetry_xsection' ] + (200,) ) )
            self.glLineWidthWorld( self.doGetLineWidthWorld() )
    
    class CircularCrossSectionCreator( StrokeDrawer ):
        def __init__( self, canvas, shape ):
            ## The implementation of this class assumes that both the silhouette and
            ## the symmetry_xsection exist.
            assert len( shape.silhouette ) > 2
            assert len( shape.symmetry_xsection ) > 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.SaveStateForUndo()
            self.shape = shape
            ## TODO: Until I can handle multiple circular cross sections, keep replacing the last one.
            del self.shape.circular_xsections[:]
            self.shape.circular_xsections.append( None )
        
        def doPointAdded( self ):
            if self.stroke[0] != self.stroke[-1]:
                ## 1 Intersect a perpendicular to the line between the stroke endpoints with both
                ##   the symmetry_xsection, to find the circular cross section's center,
                ##   and with the silhouette, to find the circular cross section's radius.
                stroke01 = asarray( self.stroke[-1] ) - asarray( self.stroke[0] )
                perp_ray = ( self.stroke[0], ( stroke01[1], -stroke01[0] ) )
                sym_isects = ray_lineloop_intersection_ts( perp_ray, self.shape.symmetry_xsection, ray_is_line = True )
                sil_isects = ray_lineloop_intersection_ts( perp_ray, self.shape.silhouette, ray_is_line = True )
                
                ## 2 Store the two closest +- intersection points in a tuple with the stroke[-1] endpoint.
                sym_isects.sort()
                sil_isects.sort()
                sym_isect_pos_t = [ t_etc for t_etc in sym_isects if t_etc[0] >= 0. ][0][0]
                sil_isect_pos_t = [ t_etc for t_etc in sil_isects if t_etc[0] >= 0. ][0][0]
                sym_isect_neg_t = [ t_etc for t_etc in sym_isects if t_etc[0] < 0. ][-1][0]
                sil_isect_neg_t = [ t_etc for t_etc in sil_isects if t_etc[0] < 0. ][-1][0]
                perp_ray = list(map( asarray, perp_ray ))
                ## Use the silhouette intersection points to find the circle's center and radius.
                ## TODO Q: Should I use a different combination of the intersections' information,
                ##         such as using the intersections with the symmetry curve to find the
                ##         circle's center, and the intersections with the silhouette to find
                ##         its radius?
                diam1 = perp_ray[0] + sil_isect_pos_t * perp_ray[1]
                diam2 = perp_ray[0] + sil_isect_neg_t * perp_ray[1]
                self.shape.circular_xsections[-1] = (
                    diam1,
                    diam2,
                    ## TODO Q: Should I use the last stroke point as the terminus means or
                    ##         the midpoint of the two intersection points + stroke01?
                    self.stroke[-1]
                    #.5 * ( diam1 + diam2 ) + stroke01
                    )
                self.shape.update_3d_rot_after()
            self.Refresh()
        
        def doStrokeFinished( self ):
            if self.shape.circular_xsections[-1] is None:
                del self.shape.circular_xsections[-1]
            
            ## We always delete some cross sections when we start this mode, so we always
            ## need to update_3d() and garbage collect and place in depth, etc.
            GarbageCollectMetaAnnotations( self.annotations.anns )
            self.shape.update_3d()
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            self.shape = None
            self.Refresh()
        
        def doDraw( self ):
            ## We override doDraw() so that we don't draw the stroke.  Instead, draw the circle.
            ## Actually, calling self.shape.update_3d_curves() means that the regular shape drawing
            ## code draws the cross section.
            return
    
    def doMouseDown( self, evt ):
        ## If we pick a handle (circular cross section or curve to overdraw), manipulate/overdraw it.
        ## Otherwise, create a new annotation.
        
        picked_handle = PickClosestSymmetryXSectionHandle(
            self.annotations.anns.of_type( wx_annotate_symmetry_xsection.SymmetryXSectionShape ),
            self.event2world( evt.GetPosition() )
            )
        if picked_handle is not None:
            #print "Picked shape, handle:", picked_handle
            if picked_handle[1][0] == SymmetryXSectionPart.SilhouetteCurve:
                forward_events( self, self.SymmetryXSectionDrawer( self, picked_handle[0] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] == SymmetryXSectionPart.CircularCrossSection:
                forward_events( self, self.CircularCrossSectionDragger( self, picked_handle[0], picked_handle[1][1] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] == SymmetryXSectionPart.SymmetryCurveInterior:
                forward_events( self, self.CircularCrossSectionCreator( self, picked_handle[0] ) )
                self.doMouseDown( evt )
        else:
            forward_events( self, self.SilhouetteDrawer( self ) )
            self.doMouseDown( evt )
    
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        if keycode in kDeleteKeys:
            
            ## evt.GetPosition() returns garbage.
            world_pos = self.event2world( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            ## If the mouse is not over any annotations, check if we're over a handle.
            ann = None
            
            ## TODO: Convenience function for this annotation picking.
            #ann = self.PickAnnotationOfType( screen_pos, annotations_ellipsoid.AnnotationData )
            t_anns = PickSymmetryXSections(
                self.annotations.anns.of_type( wx_annotate_symmetry_xsection.SymmetryXSectionShape ),
                world_pos
                )
            if len( t_anns ) > 0:
                ann = min( t_anns )[1]
            
            if ann is not None:
                self.SaveStateForUndo()
                print('deleting an annotation')
                RemoveAnnotation( self.annotations.anns, ann )
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
        
        else:
            evt.Skip()
    
    def doMouseMove( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doDraw( self ):
        def draw_child():
            child = event_child( self )
            if child is not None: child.doDraw()
            else:
                self.prepare2D()
                
                text = ''
                picked_handle = PickClosestSymmetryXSectionHandle(
                    self.annotations.anns.of_type( wx_annotate_symmetry_xsection.SymmetryXSectionShape ),
                    self.event2world( self.GetMousePosition() )
                    )
                if picked_handle is not None:
                    if picked_handle[1][0] == SymmetryXSectionPart.SilhouetteCurve:
                        text = 'Draw the symmetry cross section'
                    elif picked_handle[1][0] == SymmetryXSectionPart.CircularCrossSection:
                        text = 'Draw a circular cross section'
                    elif picked_handle[1][0] == SymmetryXSectionPart.SymmetryCurveInterior:
                        text = 'Draw a circular cross section'
                else:
                    text = 'Draw a silhouette'
                draw_mode_label_text( text )
        
        def draw_callback():
            #draw_handles()
            draw_child()
        
        self.DrawScene( callback = draw_callback )
#===================

#===================
class SweepShapePart(object):
    ## BoundaryInterior is the whole shape
    BoundaryInterior = 0
    CrossSectionInterior = 1
    CrossSection = 2
    Boundary = 3
    SymmetryPlane = 4
    Spine = 5
    Diameter = 6
    CrossSectionTiltHandle = 7
    RadiusHandle = 8
    CapRadiusHandle = 9
    SymmetryPlaneHandle = 10
def PickSweepShapeHandle( shape, canvas, world_pos, world_dir = (0,0,-1), parts = None ):
    '''
    Given a SweepShape object 'shape' and a world-space position 'world_pos' and optional
    world-space direction 'world_dir' (defaults to -z),
    returns the tuple of one or more elements ( picked SweepShapePart, per-part info... )
    corresponding to the picked part, or None if no part was picked.
    If optional parameter 'parts' is not None, then only SweepShapePart's in 'parts' will be
    considered.
    '''
    world_ray = ( zerosf( 3 ), asarrayf( world_dir ) )
    world_ray[0][:2] = world_pos[:2]
    
    if parts is None or SweepShapePart.SymmetryPlaneHandle in parts:
        ## See if we're close to a symmetry plane handle.
        for t, radians in shape.symmetry_radians.items():
            ray = shape.symmetry_plane_ray_at_t( t )
            ## Is the mouse near the point where the symmetry plane ray
            ## intersects the circular cross section?
            sym_handle_pos = ray[0] + shape.radius_at_t( t ) * dir( ray[1] )
            if mag2( world_ray[0][:2] - sym_handle_pos[:2] ) < PointDistanceThreshold( canvas )**2:
                return ( SweepShapePart.SymmetryPlaneHandle, t )
    
    if parts is None or SweepShapePart.CapRadiusHandle in parts:
        ## See if we're close to a cap radius handle.
        for t in (0,1):
            pos = shape.spine2d_at_t( t )
            tangent = shape.tangent2d_at_t( t )
            cap_radius = shape.cap_radius_at_t( t )
            
            ## Is the mouse near the point cap_radius away from the start or end point?
            cap_handle_pos = pos + cap_radius * ( -1. if t == 0 else 1 ) * dir( tangent )
            if mag2( world_ray[0][:2] - cap_handle_pos ) < PointDistanceThreshold( canvas )**2:
                return ( SweepShapePart.CapRadiusHandle, t )
    
    if parts is None or SweepShapePart.RadiusHandle in parts:
        ## For each cross section, are we near boundary perpendicular to the spine?
        for t in shape.radii.keys():
            boundary1, boundary2 = shape.boundary_perpendicular_to_t( t )
            ## Check if we're close to these perpendicular points.
            if mag2( world_ray[0][:2] - boundary1 ) < PointDistanceThreshold( canvas )**2:
                return ( SweepShapePart.RadiusHandle, t, True )
            if mag2( world_ray[0][:2] - boundary2 ) < PointDistanceThreshold( canvas )**2:
                return ( SweepShapePart.RadiusHandle, t, False )
    
    if parts is None or SweepShapePart.CrossSectionTiltHandle in parts:
        ## For each tilted cross section, are we near the intersection point of the surface
        ## with the tilt direction?
        for t in shape.tilt_radians.keys():
            ray = shape.tilt_ray_at_t( t )
            radius = shape.radius_at_t( t )
            ## Check if we're close to the intersection of the ray with the sphere at the same
            ## center and given radius.
            sphere_pt, intersected = closest_point_to_sphere_on_ray( ( ray[0], radius ), ray )
            assert intersected
            
            if mag2( world_ray[0][:2] - sphere_pt[:2] ) < PointDistanceThreshold( canvas )**2:
                return ( SweepShapePart.CrossSectionTiltHandle, t )
    
    if parts is None or SweepShapePart.Diameter in parts:
        ## Are we near a diameter line of the (perhaps approximating) circle?
        for t in shape.radii.keys():
            center = shape.spine2d_at_t( t )
            normal = dir( shape.normal2d_at_t( t ) )
            radius = shape.radius_at_t( t )
            
            distsqr, closest_line_pt, edge_index = closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
                ## A line strip composed of two points, the positive and negative
                ## normal-to-the-spine directions
                [ center + radius * normal, center - radius * normal ],
                world_ray[0][:2]
                )
            if mag2( world_ray[0][:2] - closest_line_pt ) < EdgeDistanceThreshold( canvas )**2:
                return ( SweepShapePart.Diameter, t )
    
    if parts is None or SweepShapePart.Spine in parts:
        ## Are we close to the spine curve?
        distsqr, closest_circle_pt, edge_index = closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
            shape.spine2d, world_pos[:2]
            )
        if distsqr < EdgeDistanceThreshold( canvas )**2:
            return ( SweepShapePart.Spine, annotations_sweep.edge_index2t( edge_index, shape.spine2d ) )
    
    ## TODO Q: Does anyone need this?  If so, take the drawing code that builds the
    ##         symmetry plane's line loop and find the distance to that.
    if parts is not None and SweepShapePart.SymmetryPlane in parts: raise NotImplementedError
    
    if parts is None or SweepShapePart.Boundary in parts:
        ## Are we close to the boundary curve?
        distsqr_pt_edge = closest_distsqr_and_point_and_edge_index_on_line_loop_to_point(
            shape.boundary2d, world_pos[:2]
            )
        if distsqr_pt_edge[0] < EdgeDistanceThreshold( canvas )**2:
            return ( SweepShapePart.Boundary, )
    
    if parts is None or SweepShapePart.CrossSection in parts:
        ## For each cross section (specified as a circle of given radius or arbitrary closed curve),
        ## are we near the elliptical projection of the (perhaps approximating) circle?
        for t in shape.tilt_radians.keys():
            center = shape.spine3d_at_t( t )
            normal = shape.tangent3d_at_t( t )
            radius = shape.radius_at_t( t )
            
            ## Near the ellipse boundary
            try:
                axis0 = radius * orthogonal( normal )
                axis1 = radius * dir( cross( normal, axis0 ) )
                closest_circle_pt = closest_point_on_ellipse_to_point(
                    ( center[:2], ( axis0[:2], axis1[:2] ) ), world_ray[0][:2]
                    )
            ## In case of a LinAlg error with a singular axes matrix.
            except linalg.LinAlgError:
                distsqr, closest_circle_pt, edge_index = closest_distsqr_and_point_and_edge_index_on_line_loop_to_point(
                    ## A flattened (2d by stripping z) line loop from the circle
                    [ pt[:2] for pt in line_loop_from_tilted_circle( center, radius, normal ) ],
                    world_ray[0][:2]
                    )
            if mag2( world_ray[0][:2] - closest_circle_pt ) < EdgeDistanceThreshold( canvas )**2:
                return ( SweepShapePart.CrossSection, t )
        
        for t in shape.cross_sections.keys():
            distsqr, closest_circle_pt, edge_index = closest_distsqr_and_point_and_edge_index_on_line_loop_to_point(
                ## A flattened (2d by stripping z) line loop from the circle
                [ pt[:2] for pt in shape.cross_section_line_loop3d_at_t( t ) ],
                world_ray[0][:2]
                )
            if mag2( world_ray[0][:2] - closest_circle_pt ) < EdgeDistanceThreshold( canvas )**2:
                return ( SweepShapePart.CrossSection, t )
    
    if parts is None or SweepShapePart.CrossSectionInterior in parts:
        ## For each cross section (specified as a circle of given radius or arbitrary closed curve),
        ## are we inside the elliptical projection of the (perhaps approximating) circle?
        for t in chain( iter(shape.tilt_radians.keys()), iter(shape.radii.keys()), iter(shape.cross_sections.keys()) ):
            center = shape.spine3d_at_t( t )
            normal = shape.tangent3d_at_t( t )
            radius = shape.radius_at_t( t )
            
            ## Inside the ellipse
            plane_isect = ray_plane_intersection( world_ray, ( center, normal ), ray_is_line = True )
            if plane_isect is not None and mag2( asarray( center ) - asarray( plane_isect ) ) <= radius**2:
                return ( SweepShapePart.CrossSectionInterior, t )
    
    if parts is None or SweepShapePart.BoundaryInterior in parts:
        ## Are we in the interior of the shape?
        if contour_contains( shape.boundary2d, world_ray[0][:2] ):
            return ( SweepShapePart.BoundaryInterior, )
    
    return None

def DrawPickedSweepShapeHandle( shape, handle ):
    '''
    Given a SweepShape object 'shape' and a picked handle, as would be returned
    from PickSweepShapeHandle(),
    draws a pick preview as in Illustrator's smart guide tool.
    '''
    
    if handle is None: return
    
    ## This is Illustrator CS3's Smart Guides color and line width.
    glLineWidth( 1 )
    glColor3ub( 135, 229, 228 )
    glPointSize( 5 )
    glDisable( GL_POINT_SMOOTH )
    
    if handle[0] == SweepShapePart.SymmetryPlaneHandle:
        rays = shape.symmetry_plane_rays_at_vertices()
        rays = [ ( asarrayf( ray[0] ), dir( asarrayf( ray[1] ) ) ) for ray in rays ]
        
        ## Draw a "square" running along the top and bottom of the shape, as well as the
        ## sides (end-caps).
        ## Top and bottom.
        draw_linestrip( [ ray[0] + radius * ray[1] for ray, radius in zip( rays, shape.dense_radii ) ] )
        draw_linestrip( [ ray[0] - radius * ray[1] for ray, radius in zip( rays, shape.dense_radii ) ] )
        ## Sides (at end-caps).
        glBegin( GL_LINES )
        for ray, radius in [ ( rays[0], shape.dense_radii[0] ), ( rays[-1], shape.dense_radii[-1] ) ]:
            glVertex3fv( ray[0] + radius * ray[1] )
            glVertex3fv( ray[0] - radius * ray[1] )
        glEnd()
        
        ## Draw the symmetry plane handle.
        for t in (handle[1],):
            sym_handle_i = shape.vertex_index_at_t( t )
            sym_handle_pos = rays[ sym_handle_i ][0] + shape.dense_radii[ sym_handle_i ] * rays[ sym_handle_i ][1]
            glBegin( GL_POINTS )
            glVertex3fv( sym_handle_pos )
            glEnd()
        
    elif handle[0] == SweepShapePart.CapRadiusHandle:
        for t in (handle[1],):
            pos = shape.spine2d_at_t( t )
            tangent = shape.tangent2d_at_t( t )
            cap_radius = shape.cap_radius_at_t( t )
            
            ## Is the mouse near the point cap_radius away from the start or end point?
            cap_handle_pos = pos + cap_radius * ( -1. if t == 0 else 1 ) * dir( tangent )
            glBegin( GL_POINTS )
            glVertex2fv( cap_handle_pos )
            glEnd()
    
    elif handle[0] == SweepShapePart.RadiusHandle:
        ## For each cross section, are we near boundary perpendicular to the spine?
        for t in (handle[1],):
            center2d = asarray( shape.spine2d_at_t( t ) )
            ## I can't use the radius_at_t() accessor in case the radius handle is being
            ## manipulated by the user, because the radius it returns isn't being updated in
            ## real-time.
            ## UPDATE: I can; these accessors now check to see if the value was specifically overridden.
            radius = shape.radius_at_t( t )
            normal = shape.tangent3d_at_t( t )
            
            perp = zerosf(2)
            ## 90-degree right turn
            perp[0] = normal[1]
            perp[1] = -normal[0]
            perp = dir( perp )
            perp *= radius
            
            glBegin( GL_POINTS )
            glVertex2fv( center2d + perp )
            glVertex2fv( center2d - perp )
            glEnd()
    
    elif handle[0] == SweepShapePart.CrossSectionTiltHandle:
        ## For each tilted cross section, are we near the intersection point of the surface
        ## with the tilt direction?
        for t in (handle[1],):
            center2d = asarray( shape.spine2d_at_t( t ) )
            radius = shape.radius_at_t( t )
            ## I can't use the normal because the tilt handle is being manipulated by the user,
            ## but the 3d spine isn't being updated in real-time.
            ## UPDATE: I can; these accessors check to see if the value was specifically overridden.
            #tangent3d = annotations_sweep.tilt_radians2tangent3d( tilt_radians, shape.tangent2d_at_t( t ) )
            tangent3d = shape.tangent3d_at_t( t )
            
            tilt_dir = dir( cross( cross( tangent3d, (0,0,1) ), tangent3d ) )
            
            glBegin( GL_POINTS )
            glVertex3fv( append( center2d, 0. ) + radius * tilt_dir )
            glEnd()
    
    elif handle[0] == SweepShapePart.Diameter:
        ## Are we near a diameter line of the (perhaps approximating) circle?
        for t in (handle[1],):
            center = shape.spine2d_at_t( t )
            normal = dir( shape.normal2d_at_t( t ) )
            radius = shape.radius_at_t( t )
            
            draw_linestrip( [ center + radius * normal, center - radius * normal ] )
    
    elif handle[0] in (SweepShapePart.CrossSection,SweepShapePart.CrossSectionInterior):
        ## For each cross section (specified as a circle of given radius or arbitrary closed curve),
        ## are we near the elliptical projection of the (perhaps approximating) circle?
        t = handle[1]
        
        lineloop = shape.cross_section_line_loop3d_at_t( t )
        draw_lineloop( lineloop )
    
    elif handle[0] == SweepShapePart.Spine:
        ## Are we close to the spine curve?
        draw_linestrip( shape.spine2d )
        
        t = handle[1]
        glBegin( GL_POINTS )
        glVertex2fv( shape.spine2d_at_t( t ) )
        glEnd()
        
        ## If we have an arbitrary cross section there, draw it.
        if t in shape.cross_sections:
            draw_lineloop( shape.cross_section_line_loop3d_at_t( t ) )
    
    ## TODO Q: Does anyone need this?  If so, take the drawing code that builds the
    ##         symmetry plane's line loop and find the distance to that.
    elif handle[0] == SweepShapePart.SymmetryPlane: raise NotImplementedError
    
    elif handle[0] in (SweepShapePart.Boundary, SweepShapePart.BoundaryInterior):
        ## Are we close to the boundary curve?
        draw_lineloop( shape.boundary2d )
    
    glEnable( GL_POINT_SMOOTH )

def PickClosestSweepShapeHandle( annotations, canvas, world_pos, world_dir = (0,0,-1), parts = None ):
    '''
    Similar to PickSweepShapeHandle(), except takes a sequence of SweepShape objects
    'annotations' and returns the tuple ( SweepShape object, handle ) of the most specific
    picked handle over all SweepShape objects in 'annotations'.
    '''
    
    handles = [ ( PickSweepShapeHandle( shape, canvas, world_pos, world_dir, parts ), shape ) for shape in annotations ]
    handles = [ handle for handle in handles if handle[0] is not None ]
    
    if len( handles ) == 0:
        return None
    else:
        ## Handles have enumerated values in increasing order of specificity, so return the
        ## most-specific (largest numbered) handle.
        closest_handle, closest_shape = max( handles )
        return closest_shape, closest_handle

def PickSweepShapes( annotations, world_pos, world_dir = (0,0,-1) ):
    '''
    Given a sequence of SweepShape 'annotations', a 2+d world-space position 'world_pos',
    and optional world-space vector 'world_dir' (defaults to -z),
    returns a list of unsorted (t, annotation) tuples
    corresponding to intersections of the world-ray with annotations.
    
    Note that each annotation will appear at most once, with the smallest 't'.
    Note that t's may be negative.
    '''
    
    return PickShapeSurfaces( annotations, world_pos, world_dir )

def GetSweepShapeCurrentRadius( canvas ):
    return .2 / canvas.GetZoom()

class SweepShapeAnnotationCreateMode( ModeProxy ):
    SAMPLES_PER_UNIT_DISTANCE = 20
    
    class SpineDrawer( StrokeDrawer ):
        def __init__( self, canvas ):
            ## Create a new annotation.
            self.shape = SweepShape()
            #self.annotations.anns.append( self.shape )
            
            self.shape.spine2d = []
            
            StrokeDrawer.__init__( self, canvas, self.shape.spine2d )
            ## Any lookups which might pass through to our proxied class must be done after
            ## calling the superclass constructor.
        
        def doPointAdded( self ):
            self.Refresh()
        
        def doStrokeFinished( self ):
            if len( self.shape.spine2d ) >= 3:
                
                resample_N = max( 3, int(ceil( SweepShapeAnnotationCreateMode.SAMPLES_PER_UNIT_DISTANCE * self.GetZoom() * line_strip_length( self.shape.spine2d ) )) )
                
                if self.GetMouseState().ShiftDown():
                    self.shape.spine2d = asarray(
                        lerps(
                            self.shape.spine2d[0],
                            self.shape.spine2d[-1],
                            linspace( 0, 1, resample_N )
                            )
                        )
                else:
                    self.shape.spine2d = smooth_open_stroke(
                        self.shape.spine2d,
                        resample_N = resample_N
                        )
                
                ## Add untilted circles at the first and last edges
                self.shape.tilt_radians[ 0. ] = 0.
                self.shape.tilt_radians[ 1. ] = 0.
                self.shape.radii[ 0. ] = GetSweepShapeCurrentRadius( self )
                self.shape.radii[ 1. ] = GetSweepShapeCurrentRadius( self )
                self.shape.symmetry_radians[ .5 ] = -.5*pi
                
                self.shape.update_3d()
                self.SaveStateForUndo()
                self.annotations.anns.append( self.shape )
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            #self.shape = None
            self.Refresh()
            
            ## Switch to arrow mode.
            #self.SetModeByID( self.GetModeIDForShortName( "Select" ) )
            self.SwitchToolbarTools( self.GetModeIDForShortName( "Select" ) )
            self.mode.selection.add( self.shape )
        
        def doSetGLLineStyle( self ):
            glColor4ub( *( Color[ 'spine' ] + (200,) ) )
            self.glLineWidthWorld( self.doGetLineWidthWorld() )
    
    class SpineRedrawer( StrokeDrawer ):
        def __init__( self, canvas, shape ):
            self._shape = shape
            
            StrokeDrawer.__init__( self, canvas, [], auto_unforward_events = False )
            ## Any lookups which might pass through to our proxied class must be done after
            ## calling the superclass constructor.
        
        def get_shape( self ):
            ## If the user performs an undo or redo, self.shape will no longer be in the list
            ## of annotations.  Find the shape in the list of annotations with the same uuid.
            return GetCurrentShapeState( self.annotations.anns, self._shape )
        
        shape = property( get_shape )
        
        def doPointAdded( self ):
            self.Refresh()
        
        def doStrokeFinished( self ):
            if len( self.stroke ) >= 2:
                
                spine_N = len( self.shape.spine2d )
                
                stroke = smooth_open_stroke(
                    self.stroke,
                    resample_N = max( 3, int(ceil( SweepShapeAnnotationCreateMode.SAMPLES_PER_UNIT_DISTANCE * line_strip_length( self.stroke ) )) )
                    )
                
                self.SaveStateForUndo()
                import curve_oversketching
                #self.shape.set_spine2d( stroke )
                old_ts = self.shape.get_all_particular_ts()
                new_ts = list( old_ts )
                ## If the oversketch curve and spine2d are not sampled uniformly, t's along the
                ## spine will be broken by any subsequent resampling;
                ## SweepShape's resample() and smooth_spine() methods are sensitive to this.
                new_spine = curve_oversketching.curve_from_base_curve_and_oversketch_curve(
                    self.shape.spine2d, stroke,
                    scale = self.GetZoom(),
                    ts = new_ts
                    )
                assert len( new_ts ) == len( old_ts )
                self.shape.set_spine2d( new_spine, dict( list(zip( old_ts, new_ts )) ) )
                self.shape.resample( max(
                        spine_N,
                        max( 3, int(ceil( SweepShapeAnnotationCreateMode.SAMPLES_PER_UNIT_DISTANCE * line_strip_length( new_spine ) )) )
                        ) )
                ## Smooth again, for good measure.
                ## UPDATE: Actually, smoothing on curves is pretty heavy; too heavy for oversketching.
                #self.shape.smooth_spine2d()
                self.shape.update_3d()
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
                #print "successfully oversketched"
            
            else:
                ## A click exits the mode.
                unforward_events( self )
            
            self.Refresh()
        
        def doSetGLLineStyle( self ):
            glColor4ub( *( Color[ 'spine' ] + (200,) ) )
            self.glLineWidthWorld( self.doGetLineWidthWorld() )
        
        def doDraw( self ):
            if self.shape is None:
                unforward_events( self )
                self.Refresh()
                return
            
            self.prepare2D()
            draw_mode_label_text( 'Oversketch' )
            StrokeDrawer.doDraw( self )
    
    class CircularCrossSectionTilter( StrokeDrawer ):
        def __init__( self, canvas, shape, t = None ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.SaveStateForUndo()
            self.shape = shape
            self.t = t
            if self.t is not None:
                assert self.t in self.shape.tilt_radians
                assert self.t >= 0.
                assert self.t <= 1.
        
        def doPointAdded( self ):
            ## Is this the first point?  If so, pick the closest spine point.
            if self.t is None:
                assert len( self.stroke ) == 1
                distsqr, pt, edge_index = closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
                    self.shape.spine2d, self.stroke[0]
                    )
                ## Make sure the closest edge doesn't already have a tilt.
                ## TODO Q: Should I assign one regardless, and leave it up to the shape object
                ##         to tessellate itself sufficiently?
                if distsqr < EdgeDistanceThreshold( self )**2 and edge_index not in [ self.shape.edge_index_at_t( t ) for t in self.shape.tilt_radians.keys() ]:
                    self.t = annotations_sweep.edge_index2t( edge_index, self.shape.spine2d )
                else:
                    unforward_events( self )
                    return
            
            ## We need to update tilt radians.
            self.shape.tilt_radians[ self.t ] = self.shape.compute_tilt_radians_from_point_at_t( self.t, self.stroke[-1] )
            
            if kRealTimeHacks:
                ## Find the shape with the same uuid and change its tilt radians.
                beauty = [
                    shape for shape in self.annotations.anns_beautified.of_type( SweepShape )
                    if shape.uuid == self.shape.uuid
                    ]
                ## There should be one shape with the same uuid.
                assert len( beauty ) == 1
                beauty = beauty[0]
                beauty.tilt_radians[ self.t ] = self.shape.tilt_radians[ self.t ]
                ## If no-pinch is enabled, temporarily turn it off and generate the mesh.
                no_pinch = annotations_sweep.kNoPinch
                if no_pinch:
                    annotations_sweep.kNoPinch = False
                    beauty.get_mesh()
                    annotations_sweep.kNoPinch = True
                self.GetTopLevelParent().Refresh()
            
            self.Refresh()
        
        def doStrokeFinished( self ):
            ## If we were updating the mesh in real-time, we temporarily turned of no-pinch.
            ## Now regenerate the mesh with no-pinch enabled.
            if kRealTimeHacks and annotations_sweep.kNoPinch:
                self.shape.gen_mesh()
            
            self.shape.update_3d()
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            self.shape = None
        
        def doDraw( self ):
            ## We override doDraw() so that we don't draw the stroke.  Instead, draw the circle.
            ## Actually, adding the tilted circle to self.shape.tilted_circles means that the
            ## regular shape drawing code draws the cross section.
            self.DrawSweepShapeHandles( [ self.shape ] )
    
    class CircularCrossSectionRadiusManip( StrokeDrawer ):
        def __init__( self, canvas, shape, t = None ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.SaveStateForUndo()
            self.shape = shape
            self.t = t
            if self.t is not None:
                assert self.t in self.shape.radii
                assert self.t >= 0.
                assert self.t <= 1.
        
        def doPointAdded( self ):
            ## Is this the first point?  If so, pick the closest spine point.
            if self.t is None:
                assert len( self.stroke ) == 1
                ## Cast a ray perpendicular to the closest point on the shape's boundary?
                offsets = annotations_sweep.boundary_offsets2d( self.shape.spine2d, self.shape.dense_radii )
                spine2d = asarray( self.shape.spine2d )
                distsqr, pt, edge_index = min(
                    closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
                        spine2d + offsets, self.stroke[0]
                        ),
                    closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(
                        spine2d - offsets, self.stroke[0]
                        )
                    )
                ## Make sure the closest edge doesn't already have a cross section.
                ## TODO Q: Should I assign one regardless, and leave it up to the shape object
                ##         to tessellate itself sufficiently?
                if distsqr < EdgeDistanceThreshold( self )**2 and edge_index not in [ self.shape.edge_index_at_t( t ) for t in self.shape.radii.keys() ]:
                    self.t = annotations_sweep.edge_index2t( edge_index, self.shape.spine2d )
                    ## Initialize the cross section to a list of one item.
                    ## (The list's item will be set properly below.)
                    self.shape.radii[ self.t ] = None
                else:
                    unforward_events( self )
                    return
            
            ## Update the circle's radius.
            center = self.shape.spine2d_at_t( self.t )
            radius = mag( asarray( self.stroke[-1] ) - center )
            self.shape.radii[ self.t ] = radius
            
            self.Refresh()
        
        def doStrokeFinished( self ):
            self.shape.update_3d()
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            self.shape = None
        
        def doDraw( self ):
            ## We override doDraw() so that we don't draw the stroke.  Instead, draw the circle.
            ## Actually, adding the tilted circle to self.shape.tilted_circles means that the
            ## regular shape drawing code draws the cross section.
            self.DrawSweepShapeHandles( [ self.shape ] )
    
    class CapRadiusManip( StrokeDrawer ):
        def __init__( self, canvas, shape, t = None ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.SaveStateForUndo()
            self.shape = shape
            self.t = int( t )
            assert self.t in (0,1)
        
        def doPointAdded( self ):
            ## Update the circle's radius.
            center = self.shape.spine2d_at_t( self.t )
            
            #radius = mag( asarray( self.stroke[-1] ) - center )
            
            tangent = dir( self.shape.tangent2d_at_t( self.t ) )
            if self.t == 0: tangent = -tangent
            radius = max( 0., dot( tangent, asarray( self.stroke[-1] ) - center ) )
            
            self.shape.cap_radii[ self.t ] = radius
            
            self.Refresh()
        
        def doStrokeFinished( self ):
            self.shape.update_3d()
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            self.shape = None
        
        def doDraw( self ):
            ## We override doDraw() so that we don't draw the stroke.  Instead, draw the circle.
            ## Actually, adding the tilted circle to self.shape.tilted_circles means that the
            ## regular shape drawing code draws the cross section.
            self.DrawSweepShapeHandles( [ self.shape ] )
    
    class SymmetryPlaneManip( StrokeDrawer ):
        def __init__( self, canvas, shape, t ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.SaveStateForUndo()
            self.shape = shape
            self.t = t
            assert self.t in self.shape.symmetry_radians
            assert self.t >= 0.
            assert self.t <= 1.
        
        def doPointAdded( self ):
            ## Ignore the first point, so clicks don't have an effect.
            if len( self.stroke ) == 1: return
            
            ## Update the symmetry plane tilt radians.
            ## We'll adapt SweepShape.compute_tilt_radians_from_point_at_t(), using a circle
            ## perpendicular to the one used there.
            pt = self.stroke[-1]
            
            center = self.shape.spine2d_at_t( self.t )
            radius = self.shape.radius_at_t( self.t )
            boundary1, boundary2 = self.shape.boundary_perpendicular_to_t( self.t )
            perp = boundary1 - center
            x = dot( asarrayf( pt ) - center, perp ) / radius**2
            #self.shape.symmetry_radians[ self.t ] = asin( max( -1., min( 1., x ) ) )
            #self.shape.symmetry_radians[ self.t ] = pi * max( -1., min( 1., x ) )
            ## Use asin_periodic(), but modularly keep it within [-pi,pi]
            if self.GetMouseState().AltDown():
                self.shape.symmetry_radians[ self.t ] = ( ( asin_periodic( x ) - pi ) % (2*pi) ) - pi
            else:
                self.shape.symmetry_radians[ self.t ] = asin( max( -1., min( 1., x ) ) )
            #print "symmetry_plane_radians:", self.shape.symmetry_plane_radians
            
            ## Don't apply geometry constraints in real-time, it's too expensive, apply on mouse up.
            self.shape.update_3d()
            #self.annotations.Beautify()
            #self.GetTopLevelParent().Refresh()
            
            self.Refresh()
        
        def doStrokeFinished( self ):
            ## We don't need to shape.update_3d() because the symmetry plane only affects
            ## geometry annotations.
            ## UPDATE: Actually, it also effects the way cross sections are drawn.
            self.shape.update_3d()
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            self.shape = None
        
        def doDraw( self ):
            ## We override doDraw() so that we don't draw the stroke.  Instead, draw the circle.
            ## Actually, adding the tilted circle to self.shape.tilted_circles means that the
            ## regular shape drawing code draws the cross section.
            self.DrawSweepShapeHandles( [ self.shape ] )
    
    class SilhouetteRedrawer( StrokeDrawer ):
        def __init__( self, canvas, shape ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.SaveStateForUndo()
            self.shape = shape
        
        def doPointAdded( self ):
            self.Refresh()
        
        def doStrokeFinished( self ):
            ## Convert the stroke to radii for all cross sections.
            ## TODO Q: Change to flat end-caps?
            silhouette = smooth_closed_stroke( self.stroke )
            #self.shape.set_radii_from_silhouette( silhouette )
            self.shape.set_spine_and_radii_from_silhouette( silhouette )
            
            self.shape.update_3d()
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            self.shape = None
        
        def doStrokeFinished( self ):
            silhouette = smooth_closed_stroke( self.stroke )
            
            def iterative_silhouette( shape, silhouette, count ):
                if self.canvas.GetShortNameForModeID( self.canvas.GetCurrentModeID() ).lower() != "select":
                    print('no more iterating silhouette fit')
                    return
                
                print('iterating silhouette fit')
                shape.set_spine_and_radii_from_silhouette( silhouette, count )
                shape.update_3d()
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
                
                from wx_util import Timers
                Timers.RequestTimer( lambda: iterative_silhouette( shape, silhouette, count+1 ), 100 )
            
            iterative_silhouette( self.shape, silhouette, 0 )
            self.shape = None
        
        def doSetGLLineStyle( self ):
            glColor4ub( *( Color[ 'silhouette' ] + (200,) ) )
            self.glLineWidthWorld( self.doGetLineWidthWorld() )
    
    class ArbitraryCrossSection( ModeProxy ):
        def __init__( self, canvas, shape, t, immediate_click = True ):
            ModeProxy.__init__( self, canvas )
            
            ## When switching into this mode, the parent event handler immediately calls
            ## doMouseDown(), which we're not interested in.
            self.ate_first_mouse_down = not immediate_click
            
            ## This class doesn't interact well with undo/redo.  Because we are a persistent
            ## mode, undo/redo change the list of shapes out from underneath us.  Our edits
            ## are pushed onto the undo stack, but if the user ever chooses undo or redo
            ## while we are in this mode, then self.shape will point somewhere bad
            ## (nowhere I think, but even worse if it points to a saved state).
            ## UPDATE: Now that I have uuid's per shape, I can track them across undo/redo,
            ##         so this is obsolete.
            #self.EditMenuUndoRedoDisallow()
            
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            self._shape = shape
            self.t = t
            ## t may not yet be in shape.cross_sections, if we're adding it here.
            #assert self.t in self.shape.cross_sections
            assert self.t >= 0.
            assert self.t <= 1.
            
            self.characteristic_xsections = []
            #self.characteristic_xsections.append( annotations_sweep.unit_circle )
            #self.characteristic_xsections.append( annotations_sweep.unit_circle )
            #self.characteristic_xsections.append( fromfile( os.path.join( '.', 'xsections', 'unit_circle' ) ) )
            #self.characteristic_xsections.append( fromfile( os.path.join( '.', 'xsections', 'ear' ) ) )
            from glob import glob
            for xsection_file in glob( os.path.join( '.', 'xsections', '*.xsection' ) ):
                #print xsection_file
                xsection = fromfile( xsection_file )
                assert ( len( xsection ) % 2 ) == 0
                xsection = xsection.reshape( ( len( xsection ) // 2, 2 ) )
                self.characteristic_xsections.append( xsection )
            ## These will get populated at draw time.
            self.characteristic_xsection_rects = []
        
        def get_shape( self ):
            ## If the user performs an undo or redo, self.shape will no longer be in the list
            ## of annotations.  Find the shape in the list of annotations with the same uuid.
            return GetCurrentShapeState( self.annotations.anns, self._shape )
        
        shape = property( get_shape )
        
        def doPickRotateHandle( self, eye_pos ):
            return mag2( asarray( eye_pos ) - self.kRotateHandlePos ) < 2.5*kPointDistanceThreshold**2
        
        def doMouseDown( self, evt ):
            ## When switching into this mode, the parent event handler immediately calls
            ## doMouseDown(), which we're not interested in.
            if not self.ate_first_mouse_down:
                self.ate_first_mouse_down = True
                return
            
            if self.shape is None:
                unforward_events( self )
                self.Refresh()
                return
            
            eye_pos = self.event2eye( evt.GetPosition() )
            
            ## If the user clicks in the canvas, switch to curve drawing mode.
            ## TODO: Add a curve deforming mode.
            if point_is_inside_rect( eye_pos, ( -self.kCanvasOuter, -self.kCanvasOuter, self.kCanvasOuter, self.kCanvasOuter ) ):
                ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
                ## contextual menu event, but the right mouse button will be down.
                if evt.RightDown():
                    print('skipping mouse down')
                    evt.Skip()
                    return
                
                forward_events( self, SweepShapeAnnotationCreateMode.ArbitraryCrossSectionDrawer( self, self.shape, self.t ) )
                self.doMouseDown( evt )
                return
            
            ## If the user clicks on the rotate handle, switch to rotating the stroke.
            if self.doPickRotateHandle( eye_pos ):
                #print "rotate handle"
                forward_events( self, SweepShapeAnnotationCreateMode.ArbitraryCrossSectionRotator( self, self.shape, self.t ) )
                self.doMouseDown( evt )
                return
            
            ## If the user clicks a pre-defined cross section, set it as the cross section.
            for xsection, rect in zip( self.characteristic_xsections, self.characteristic_xsection_rects ):
                if point_is_inside_rect( eye_pos, rect ):
                    self.SaveStateForUndo()
                    
                    self.shape.cross_sections[ self.t ] = array( xsection )
                    
                    self.shape.update_3d()
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                    return
            
            ## If the user clicks in the background, exit this mode.
            ## Save the cross section if we have one.
            if self.t in self.shape.cross_sections:
                self.shape.cross_sections[ self.t ].tofile( 'last xsection.xsection' )
            unforward_events( self )
            #self.EditMenuUndoRedoAllow()
            self.Refresh()
        
        def doMouseMove( self, evt ):
            ## Refresh so we can give hover feedback.
            self.Refresh()
        
        def doContextMenu( self, evt ):
            menu = wx.Menu()
            menu_ids = []
            
            ## NOTE: I can't use evt.GetPosition() because it returns bogus data here.
            #eye_pos = self.event2eye( self.GetMousePosition() )
            
            ## Right now, the contextual menu can only be brought up if the user right-clicks
            ## inside the current cross section stroke view.
            ## So, we don't need to pick under the mouse to determine the contextual menu items.
            
            if len( self.shape.dense_cross_sections ) != 0:
                def smooth_arbitrary_cross_section( evt ):
                    self.SaveStateForUndo()
                    current_xsection = self.shape.cross_section_at_t( self.t )
                    self.shape.cross_sections[ self.t ] = annotations_sweep.normalize_arbitrary_cross_section(
                        smooth_closed_stroke(
                            current_xsection,
                            resample_N = len( current_xsection )
                            )
                        )
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    self.shape.update_3d()
                    
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                
                item = menu.Append( wx.ID_ANY, "Smooth Cross Section" )
                self.Bind( wx.EVT_MENU, smooth_arbitrary_cross_section, id = item.GetId() )
                menu_ids.append( item.GetId() )
            
            if self.t in self.shape.cross_sections:
                ## Add delete cross section to the menu.
                def delete_arbitrary_cross_section( evt ):
                    self.SaveStateForUndo()
                    del self.shape.cross_sections[ self.t ]
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    self.shape.update_3d()
                    
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                    
                    unforward_events( self )
                    #self.EditMenuUndoRedoAllow()
                
                if len( menu_ids ) > 0: menu.AppendSeparator()
                item = menu.Append( wx.ID_ANY, "Delete Cross Section" )
                self.Bind( wx.EVT_MENU, delete_arbitrary_cross_section, id = item.GetId() )
                menu_ids.append( item.GetId() )
            
            self.PopupMenu( menu )
            menu.Destroy()
            ## I'm not sure where I got this idea of unbinding these events.
            ## I guess it's just good cleanup, so we don't leave a dangling bound callback to a menu
            ## event that can no longer be called, and also to free the id's we reserved so we
            ## don't run out of id space, since every time this is called we generate new ids.
            for id in menu_ids: self.Unbind( wx.EVT_MENU, id = id )
        
        def doDraw( self ):
            if self.shape is None:
                unforward_events( self )
                self.Refresh()
                return
            
            ## Draw the canvas, the characteristic cross section buttons, the smooth button,
            ## and the rotate slider.
            
            self.prepare2D()
            
            glMatrixMode( GL_PROJECTION )
            glPushMatrix()
            glLoadIdentity()
            
            self.drawCharacteristicCrossSections()
            self.drawCanvas()
            
            #self.drawCurrentStroke( 200./255. )
            #self.drawCurrentStroke()
            self.drawCurrentStroke( 1. )
            
            child = event_child( self )
            if child is not None: child.doDraw()
            
            glMatrixMode( GL_PROJECTION )
            glPopMatrix()
        
        kCanvasOuter = .52
        kCanvasInner = .5
        kRotateHandlePos = ( kCanvasOuter + .11, 0. )
        
        def drawCanvas( self ):
            ## Draw the white canvas.
            glColor4f( 0,0,0, .8 )
            glRectf( -self.kCanvasOuter, -self.kCanvasOuter, self.kCanvasOuter, self.kCanvasOuter )
            glColor4f( 1,1,1, .8 )
            glRectf( -self.kCanvasInner, -self.kCanvasInner, self.kCanvasInner, self.kCanvasInner )
            
            ## Draw the symmetry direction.
            glColor3ub( *Color['symmetry_sheet'] )
            self.glLineWidthWorld( WorldDefaultLineWidth( self ) )
            glBegin( GL_LINES )
            glVertex2f( 0,0 )
            glVertex2f( self.kCanvasInner, 0 )
            glEnd()
            
            eye_pos = self.event2eye( self.GetMousePosition() )
            if ( self.doPickRotateHandle( eye_pos ) and event_child( self ) is None ) or isinstance( event_child( self ), SweepShapeAnnotationCreateMode.ArbitraryCrossSectionRotator ):
                glColor3f( .5, .5, 1 )
            else:
                glColor4f( .5, .5, 1, .75 )
            ## Draw the rotate handle
            draw_widget(
                'cross_section_tilt',
                asarray( self.kRotateHandlePos ) - ( .06, 0. ),
                ( 0, 1 ),
                4*kPointDistanceThreshold
                )
        
        def drawCurrentStroke( self, alpha = 50./255. ):
            ## Draw the current stroke (if re-drawing) or a circle.
            current_xsection = self.shape.cross_section_at_t( self.t )
            
            glColor4ub( *( Color[ 'symmetry_xsection' ] + (int(round(alpha*255)),) ) )
            draw_thick_lineloop2d( self.kCanvasInner * current_xsection, 2*kWorldDefaultLineWidth )
        
        kCharacteristicCrossSectionOuter = .22
        kCharacteristicCrossSectionInner = .2
        def drawCharacteristicCrossSections( self ):
            
            width_of_one = 2*self.kCharacteristicCrossSectionOuter
            padding = ( self.kCharacteristicCrossSectionOuter - self.kCharacteristicCrossSectionInner )
            
            total_width = len( self.characteristic_xsections) * width_of_one + ( len( self.characteristic_xsections ) - 1 ) * padding
            
            origin = array( (
                -total_width / 2. + width_of_one / 2.,
                #-self.kCanvasOuter - (2*self.kCharacteristicCrossSectionOuter - self.kCharacteristicCrossSectionInner),
                self.kCanvasOuter + (2*self.kCharacteristicCrossSectionOuter - self.kCharacteristicCrossSectionInner)
                ) )
            
            self.characteristic_xsection_origins = []
            for xsection in self.characteristic_xsections:
                #self.characteristic_xsection_origins.append( origin.copy() )
                
                outer_rect = (
                    origin[0] - self.kCharacteristicCrossSectionOuter, origin[1] - self.kCharacteristicCrossSectionOuter,
                    origin[0] + self.kCharacteristicCrossSectionOuter, origin[1] + self.kCharacteristicCrossSectionOuter
                    )
                self.characteristic_xsection_rects.append( outer_rect )
                glColor4f( 0,0,0,.8 )
                glRectf( *outer_rect )
                glColor4f( 1,1,1,.8 )
                inner_rect = (
                    origin[0] - self.kCharacteristicCrossSectionInner, origin[1] - self.kCharacteristicCrossSectionInner,
                    origin[0] + self.kCharacteristicCrossSectionInner, origin[1] + self.kCharacteristicCrossSectionInner
                    )
                glRectf( *inner_rect )
                
                glMatrixMode( GL_MODELVIEW )
                glPushMatrix()
                
                glTranslatef( origin[0], origin[1], 0. )
                glColor4ub( *( Color[ 'symmetry_xsection' ] + (200,) ) )
                draw_thick_lineloop2d( self.kCharacteristicCrossSectionInner * xsection, kWorldDefaultLineWidth )
                
                glMatrixMode( GL_MODELVIEW )
                glPopMatrix()
                
                origin[0] += self.kCharacteristicCrossSectionOuter * 2 + (self.kCharacteristicCrossSectionOuter - self.kCharacteristicCrossSectionInner)
    
    class ArbitraryCrossSectionRotator( ModeProxy ):
        def __init__( self, canvas, shape, t ):
            ModeProxy.__init__( self, canvas )
            
            ## We can't assert this, because we will create it as needed
            ## in doMouseDown().
            # assert t in shape.cross_sections
            self.shape = shape
            self.t = t
        
        def doMouseDown( self, evt ):
            self.SaveStateForUndo()
            self.last_pos = self.event2eye( evt.GetPosition() )
            
            ## If the user tries to rotate the cross section at a point that
            ## hasn't had a cross section specified yet, set it.
            ## TODO Q: Should cross section rotation follow the symmetry sheet?
            if self.t not in self.shape.cross_sections:
                self.SaveStateForUndo()
                self.shape.cross_sections[ self.t ] = self.shape.cross_section_at_t( self.t )
        
        def doMouseDrag( self, evt ):
            cur_pos = self.event2eye( evt.GetPosition() )
            
            xsection = self.shape.cross_sections[ self.t ]
            rot = get_rotation_matrix( *get_axis_angle_rotation(
                append( self.last_pos, 0. ), append( cur_pos, 0. )
                ) )[:2,:2]
            self.shape.cross_sections[ self.t ] = dot( rot, xsection.T ).T
            
            self.last_pos = cur_pos
            
            self.Refresh()
        
        def doMouseUp( self, evt ):
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            unforward_events( self )
    
    class ArbitraryCrossSectionDrawer( StrokeDrawer ):
        def __init__( self, canvas, shape, t ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.shape = shape
            self.t = t
            ## t may not yet be in shape.cross_sections, if we're adding it here.
            #assert self.t in self.shape.cross_sections
            assert self.t >= 0.
            assert self.t <= 1.
        
        def event2stroke( self, event_pos ):
            return self.event2eye( event_pos )
        
        def doPointAdded( self ):
            self.Refresh()
        
        def doStrokeFinished( self ):
            ## Abort if the user didn't draw a loop.
            if len( self.stroke ) < 3:
                print("Sketching cross section cancelled: no stroke drawn.")
                self.shape = None
                self.GetTopLevelParent().Refresh()
                return
            
            self.SaveStateForUndo()
            
            xsection = asarray( smooth_closed_stroke( self.stroke, resample_N = annotations_sweep.kCrossSectionSamples ) )
            ## normalize_arbitrary_cross_section() translates to the center-of-mass and scales to
            ## a unit size, so we don't need to worry about the world-space coordinates zoom or
            ## centering.
            xsection = annotations_sweep.normalize_arbitrary_cross_section( xsection )
            
            self.shape.cross_sections[ self.t ] = xsection
            
            self.shape.update_3d()
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            self.shape = None
        
        def doSetGLLineStyle( self ):
            glColor4ub( *( Color[ 'symmetry_xsection' ] + (200,) ) )
            self.glLineWidthWorld( 2*kWorldDefaultLineWidth )
    
    class Translater( StrokeDrawer ):
        def __init__( self, canvas, shape ):
            ## The implementation of this class assumes that the spine exists.
            assert len( shape.spine2d ) >= 2
            
            StrokeDrawer.__init__( self, canvas, [] )
            self.shape = shape
        
        def doPointAdded( self ):
            if len( self.stroke ) < 2: return
            elif len( self.stroke ) == 2: self.SaveStateForUndo()
            
            delta = asarray( self.stroke[-1] ) - asarray( self.stroke[-2] )
            self.shape.translate_2d( delta )
            
            self.Refresh()
        
        def doStrokeFinished( self ):
            if len( self.stroke ) >= 2:
                ## Actually, I don't need to call this for translate_2d().
                #self.shape.update_3d()
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            else:
                self.Refresh()
            
            self.shape = None
        
        def doDraw( self ):
            ## Don't draw anything, we're translating the shape, so normal drawing shows all.
            pass
    
    def doMouseDown( self, evt ):
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return
        
        ## If we pick a handle (circular cross section or curve to overdraw), manipulate/overdraw it.
        ## Otherwise, create a new annotation.
        
        picked_handle = None
        parts = None
        ## If the command key is down, always draw a new spine.
        if evt.CmdDown():
            parts = []
        ## If the shift key is down, always move.
        elif evt.ShiftDown():
            parts = [ SweepShapePart.BoundaryInterior ]
        picked_handle = PickClosestSweepShapeHandle(
            self.annotations.anns.of_type( SweepShape ),
            self,
            self.event2world( evt.GetPosition() ),
            parts = parts
            )
        if picked_handle is not None:
            #print "Picked shape, handle:", picked_handle
            if picked_handle[1][0] == SweepShapePart.CrossSectionTiltHandle:
                forward_events( self, self.CircularCrossSectionTilter( self, picked_handle[0], picked_handle[1][1] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] == SweepShapePart.RadiusHandle:
                forward_events( self, self.CircularCrossSectionRadiusManip( self, picked_handle[0], picked_handle[1][1] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] == SweepShapePart.SymmetryPlaneHandle:
                forward_events( self, self.SymmetryPlaneManip( self, picked_handle[0], picked_handle[1][1] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] == SweepShapePart.Spine:
                forward_events( self, self.CircularCrossSectionTilter( self, picked_handle[0] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] == SweepShapePart.Boundary:
                forward_events( self, self.CircularCrossSectionRadiusManip( self, picked_handle[0] ) )
                self.doMouseDown( evt )
            elif picked_handle[1][0] in (SweepShapePart.BoundaryInterior,):
                forward_events( self, self.Translater( self, picked_handle[0] ) )
                self.doMouseDown( evt )
        else:
            forward_events( self, self.SpineDrawer( self ) )
            self.doMouseDown( evt )
    
    def doContextMenu( self, evt ):
        ## NOTE: I can't use evt.GetPosition() because it returns bogus data here.
        world_pos = self.event2world( self.GetMousePosition() )
        
        ## If we pick a SweepShape, show the contextual menu.
        picked_handle = PickClosestSweepShapeHandle(
            self.annotations.anns.of_type( SweepShape ),
            self,
            world_pos
            )
        if picked_handle is None:
            evt.Skip()
            return
        
        shape = picked_handle[0]
        
        menu = wx.Menu()
        menu_ids = []
        
        def redraw_silhouette( evt ):
            forward_events( self, self.SilhouetteRedrawer( self, shape ) )
        
        item = menu.Append( wx.ID_ANY, "Redraw Silhouette" )
        self.Bind( wx.EVT_MENU, redraw_silhouette, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        def redraw_spine( evt ):
            forward_events( self, self.SpineRedrawer( self, shape ) )
        
        item = menu.Append( wx.ID_ANY, "Redraw Spine" )
        self.Bind( wx.EVT_MENU, redraw_spine, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        if True: #wx_annotate_shell.kDebugging:
            def smooth():
                self.SaveStateForUndo()
                shape.smooth_spine2d()
                GarbageCollectMetaAnnotations( self.annotations.anns )
                shape.update_3d()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Smooth Spine" )
            self.Bind( wx.EVT_MENU, lambda evt: smooth(), id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            
            def resample( factor ):
                self.SaveStateForUndo()
                shape.resample( max( 3, int(round( factor * len( shape.spine2d ) )) ) )
                GarbageCollectMetaAnnotations( self.annotations.anns )
                shape.update_3d()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Double Spine Tesselation" )
            self.Bind( wx.EVT_MENU, lambda evt: resample( 2 ), id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            item = menu.Append( wx.ID_ANY, "Halve Spine Tesselation" )
            self.Bind( wx.EVT_MENU, lambda evt: resample( .5 ), id = item.GetId() )
            menu_ids.append( item.GetId() )
        
        ## Did we pick a tilt or radius manipulation handle?  If so, add its items to the contextual menu.
        picked = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.CrossSectionTiltHandle,) )
        if picked is not None:
            t = picked[1]
            
            menu.AppendSeparator()
            
            ## Add delete cross section to the menu.
            ## Make delete disabled if we don't want the user to delete the cross section
            ## (start or end cross section or only one left).
            
            def delete_circle( evt, t = t ):
                self.SaveStateForUndo()
                del shape.tilt_radians[ t ]
                GarbageCollectMetaAnnotations( self.annotations.anns )
                shape.update_3d()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Delete Tilt Handle" )
            self.Bind( wx.EVT_MENU, delete_circle, id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            ## Don't allow users to delete the start or end cross sections.
            ## UPDATE: Actually, I don't care.
            ## UPDATE 2: Radii falloff to zero or less with only one cross section on either end.
            ## TODO: Fix radii falling-off with only one cross section on either end.
            #if t == 0. or t == 1.: item.Enable( False )
            
            ## Don't allow users to delete the last remaining cross section.
            if len( shape.tilt_radians ) == 1: item.Enable( False )
        
        picked = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.RadiusHandle,) )
        if picked is not None:
            t = picked[1]
            
            menu.AppendSeparator()
            
            ## Add delete cross section to the menu.
            ## Make delete disabled if we don't want the user to delete the cross section
            ## (start or end cross section or only one left).
            
            def delete_circle( evt, t = t ):
                self.SaveStateForUndo()
                del shape.radii[ t ]
                GarbageCollectMetaAnnotations( self.annotations.anns )
                shape.update_3d()
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            item = menu.Append( wx.ID_ANY, "Delete Scale Handle" )
            self.Bind( wx.EVT_MENU, delete_circle, id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            ## Don't allow users to delete the start or end cross sections.
            ## UPDATE: Actually, I don't care.
            ## UPDATE 2: Radii falloff to zero or less with only one cross section on either end.
            ## TODO: Fix radii falling-off with only one cross section on either end.
            #if t == 0. or t == 1.: item.Enable( False )
            
            ## Don't allow users to delete the last remaining cross section.
            if len( shape.radii ) == 1: item.Enable( False )
        
        ## Did we pick the spine?  If so, allow the user to draw an arbitrary cross section.
        ## If there already is an arbitrary cross section draw, allow the user to delete it.
        ## TODO: Add an arbitrary cross section "pencil" handle for control-clicking on.
        ##       Or this will be handled by proper selection UI for cross sections.
        picked = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.Spine,) )
        ## If we don't pick the spine, also check if we picked an existing arbitrary cross section.
        ## TODO: Arbitrary cross sections are not drawn as circles, but they are picked as circles.
        ##       Fix this so that they are picked as arbitrary cross sections.
        if picked is None:
            picked = PickSweepShapeHandle( shape, self, world_pos, parts = (SweepShapePart.CrossSectionInterior, SweepShapePart.CrossSection) )
            ## Make sure this picked cross section is really an arbitrary one.
            if picked is not None and picked[1] not in shape.cross_sections: picked = None
        if picked is not None:
            t = picked[1]
            
            menu.AppendSeparator()
            
            ## Add draw arbitrary cross section to the menu.
            
            def draw_arbitrary_cross_section( evt, t = t ):
                ## Push the new UI mode.
                forward_events( self, self.ArbitraryCrossSection( self, shape, t, immediate_click = False ) )
                
                '''
                ## For now, just insert a dummy cross section.
                self.SaveStateForUndo()
                #shape.cross_sections[ t ] = annotations_sweep.unit_circle
                ## A cardioid.
                r = 1./3.
                shape.cross_sections[ t ] = asarray( [ (
                        2*r * ( cos( theta ) - .5 * cos( 2*theta ) ),
                        2*r * ( sin( theta ) - .5 * sin( 2*theta ) )
                    ) for theta in linspace( 0, 2*pi, kCrossSectionSamples )[:-1] ] )
                del r
                pass
                
                GarbageCollectMetaAnnotations( self.annotations.anns )
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
                '''
            
            if t in shape.cross_sections:
                txt = "Re-Draw Arbitrary Cross Section"
            else:
                txt = "Draw Arbitrary Cross Section"
            item = menu.Append( wx.ID_ANY, txt )
            self.Bind( wx.EVT_MENU, draw_arbitrary_cross_section, id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            
            ## Add smooth and delete cross-section to the menu.
            if t in shape.cross_sections:
                ## Add smooth cross section to the menu.
                def smooth_arbitrary_cross_section( evt, t = t ):
                    self.SaveStateForUndo()
                    shape.cross_sections[ t ] = annotations_sweep.normalize_arbitrary_cross_section(
                        smooth_closed_stroke(
                            shape.cross_sections[ t ],
                            resample_N = len( shape.cross_sections[ t ] )
                            )
                        )
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    shape.update_3d()
                    
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                
                item = menu.Append( wx.ID_ANY, "Smooth Arbitrary Cross Section" )
                self.Bind( wx.EVT_MENU, smooth_arbitrary_cross_section, id = item.GetId() )
                menu_ids.append( item.GetId() )
                
                ## Add delete cross section to the menu.
                def delete_arbitrary_cross_section( evt, t = t ):
                    self.SaveStateForUndo()
                    del shape.cross_sections[ t ]
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    shape.update_3d()
                    
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
                
                item = menu.Append( wx.ID_ANY, "Delete Arbitrary Cross Section" )
                self.Bind( wx.EVT_MENU, delete_arbitrary_cross_section, id = item.GetId() )
                menu_ids.append( item.GetId() )
        
        self.PopupMenu( menu )
        menu.Destroy()
        ## I'm not sure where I got this idea of unbinding these events.
        ## I guess it's just good cleanup, so we don't leave a dangling bound callback to a menu
        ## event that can no longer be called, and also to free the id's we reserved so we
        ## don't run out of id space, since every time this is called we generate new ids.
        for id in menu_ids: self.Unbind( wx.EVT_MENU, id = id )
    
    def doKeyUp( self, evt ):
        ## Some operations change depending on modifier state.
        self.Refresh()
        evt.Skip()
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        ## Some operations change depending on modifier state.
        self.Refresh()
        
        if keycode in kDeleteKeys:
            
            ## evt.GetPosition() returns garbage.
            world_pos = self.event2world( self.GetMousePosition() )
            
            ## First check if the mouse is over an annotation.
            ## If the mouse is not over any annotations, check if we're over a handle.
            ann = None
            
            ## TODO: Convenience function for this annotation picking.
            #ann = self.PickAnnotationOfType( screen_pos, annotations_ellipsoid.AnnotationData )
            t_anns = PickSweepShapes(
                self.annotations.anns.of_type( SweepShape ),
                world_pos
                )
            if len( t_anns ) > 0:
                ann = min( t_anns )[1]
            
            if ann is not None:
                self.SaveStateForUndo()
                print('deleting an annotation')
                RemoveAnnotation( self.annotations.anns, ann )
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
        
        else:
            evt.Skip()
    
    def doMouseMove( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doDraw( self ):
        def draw_child():
            child = event_child( self )
            if child is not None:
                child.doDraw()
                return True
        
        def draw_callback():
            if draw_child(): return
            
            self.prepare2D()
            
            text = ''
            
            picked_handle = None
            parts = None
            ## If the command key is down, always draw a new spine.
            if self.GetMouseState().CmdDown():
                parts = []
            ## If the shift key is down, always move.
            elif self.GetMouseState().ShiftDown():
                parts = [ SweepShapePart.BoundaryInterior ]
            picked_handle = PickClosestSweepShapeHandle(
                self.annotations.anns.of_type( SweepShape ),
                self,
                self.event2world( self.GetMousePosition() ),
                parts = parts
                )
            if picked_handle is not None:
                draw_handles = [ picked_handle[0] ]
                ## XXX TODO: This is ugly.
                #if wx_annotate_shell.kPaperFigureMode: draw_handles = None
                self.DrawSweepShapeHandles( draw_handles )
                
                if picked_handle[1][0] == SweepShapePart.CrossSectionTiltHandle:
                    text = 'Change tilt'
                elif picked_handle[1][0] == SweepShapePart.RadiusHandle:
                    ## This text is too long
                    #text = 'Change circular cross section radius'
                    text = 'Change thickness'
                elif picked_handle[1][0] == SweepShapePart.SymmetryPlaneHandle:
                    text = 'Change symmetry plane'
                elif picked_handle[1][0] == SweepShapePart.Spine:
                    text = 'Tilt'
                elif picked_handle[1][0] == SweepShapePart.Boundary:
                    text = 'Change thickness'
                elif picked_handle[1][0] in (SweepShapePart.BoundaryInterior,):
                    text = 'Move'
            else:
                text = 'Draw a spine'
            draw_mode_label_text( text )
        
        self.DrawScene( callback = draw_callback )
#===================

#===================
class ConnectionCurveMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
        self.excluded = set()
        
        self.kConnectionManipulationMode = 'drag'
        assert self.kConnectionManipulationMode in ('draw', 'drag', 'draw-then-drag')
    
    def doMouseDown( self, evt ):
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return self
        
        world_pos = self.event2world( self.GetMousePosition() )
        ## If we don't pick exactly two regions under the mouse, don't draw a connection.
        shapes = [ ann for t,ann in PickShapeSurfaces( [ ann for ann in self.annotations.anns if ann not in self.locked_shapes ], append( world_pos, 0. ) ) ]
        ## Ignore regions in self.excluded.
        shapes = [ shape for shape in shapes if shape not in self.excluded ]
        if len( shapes ) != 2: return False
        
        ## If there's already a connection between the two regions, overdraw it.
        ## TODO Q: Should this drag the connection curve instead?
        previous_connection = [ con for con in self.annotations.anns.of_type( Connection ) if set( con.regions ) == set( shapes ) ]
        assert len( previous_connection ) <= 1
        previous_connection = None if len( previous_connection ) == 0 else previous_connection[0]
        
        ## Check if we would create a cycle in the connection graph; if so, bail.
        if previous_connection is None:
            component = set()
            def visit( con, node ):
                component.add( node )
            run_dfs( self.annotations.anns.of_type( Connection ), pre_visit = visit, roots = [shapes[0]] )
            if shapes[1] in component:
                print("Refusing to add a cycle to the graph of connected shapes.")
                return
        
        ## Forward events to the connection manipulation class.
        if self.kConnectionManipulationMode == 'draw':
            ## Stroke drawing:
            forward_events( self, self.ConnectionDrawer( self, shapes, previous_connection ) )
        elif self.kConnectionManipulationMode == 'drag':
            ## Curve dragging:
            forward_events( self, self.ConnectionDragger( self, shapes, previous_connection ) )
        elif self.kConnectionManipulationMode == 'draw-then-drag':
            if previous_connection is None:
                ## Stroke drawing:
                forward_events( self, self.ConnectionDrawer( self, shapes, previous_connection ) )
            else:
                ## Curve dragging:
                forward_events( self, self.ConnectionDragger( self, shapes, previous_connection ) )
        else:
            raise RuntimeError("Unknown kConnectionManipulationMode: " + self.kConnectionManipulationMode)
        self.doMouseDown( evt )
    
    class ConnectionDrawer( StrokeDrawer ):
        def __init__( self, canvas, shapes, previous_connection = None ):
            StrokeDrawer.__init__( self, canvas, [] )
            self.SaveStateForUndo()
            
            assert len( shapes ) == 2
            
            self.shapes = shapes
            self.previous_connection = previous_connection
        
        def doPointAdded( self ):
            self.Refresh()
        
        def doStrokeFinished( self ):
            if len( self.stroke ) > 1 and tuple( self.stroke[0] ) != tuple( self.stroke[-1] ):
                
                stroke = smooth_open_stroke( self.stroke, resample_N = 10 )
                if self.previous_connection:
                    self.previous_connection.ReplaceConnectionCurve( stroke )
                else:
                    connection = Connection( stroke, self.shapes[0], True, self.shapes[1], True )
                    self.annotations.anns.append( connection )
                
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            self.Refresh()
        
        def doSetGLLineStyle( self ):
            glColor4ub( *( Color[ 'connection' ] + (200,) ) )
            self.glLineWidthWorld( self.doGetLineWidthWorld() )
    
    class ConnectionDragger( ModeProxy ):
        def __init__( self, canvas, shapes, previous_connection = None ):
            ModeProxy.__init__( self, canvas )
            self.SaveStateForUndo()
            
            assert len( shapes ) == 2
            
            self.shapes = shapes
            self.connection = previous_connection
            if self.connection is None:
                self.connection = Connection( [], self.shapes[0], True, self.shapes[1], True )
                self.annotations.anns.append( self.connection )
        
        def doMouseDown( self, evt ):
            world_pos = self.event2world( evt.GetPosition() )
            self.SetConnectionFromPoint( world_pos )
            
            self.Refresh()
        
        def doMouseDrag( self, evt ):
            world_pos = self.event2world( evt.GetPosition() )
            self.SetConnectionFromPoint( world_pos )
            
            self.Refresh()
            
            ## Real-time updates are too slow.
            #self.annotations.Beautify()
            #self.GetTopLevelParent().Refresh()
        
        def doMouseUp( self, evt ):
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
            
            unforward_events( self )
        
        def SetConnectionFromPoint( self, world_pos ):
            assert len( world_pos ) == 2
            
            ## 1 Pick both meshes at world_pos.
            ## 2 Translate the meshes in depth so that the points have the same z.
            ## 3 Set the connection curve to the mesh-mesh intersection curve (in 2D).
            
            ## 1
            ray = ( append( world_pos, 0. ), (0,0,-1) )
            ts0 = self.connection.regions[0].ray_surface_intersection_ts( ray, ray_is_line = True )
            ts1 = self.connection.regions[1].ray_surface_intersection_ts( ray, ray_is_line = True )
            
            if len( ts0 ) < 2 or len( ts1 ) < 2:
                print("Ray doesn't intersect meshes")
                return
            
            ts0.sort()
            ts1.sort()
            
            ## NOTE: These depth values are along the -z direction.
            depth0 = ts0[ 0 if self.connection.regions_attaches_front[0] else 1 ]
            depth1 = ts1[ 0 if self.connection.regions_attaches_front[1] else 1 ]
            
            ## 2
            ## The depths are along the -z direction, so if we want the regions to coincide under
            ## world_pos, translate region1 by ( depth1 - depth0 ).
            z_offset = depth1 - depth0
            ## TODO Q: Can I just call offset_in_depth()?  Or is there some hidden dependency on
            ##    the depth, such as the initial depth placement of the root result shape?
            ## A1: I think the initial depth placement of the root result shape is dependent;
            ##     I'm certainly not taking this optimization while debugging.
            ## NOTE: offset_in_depth() offsets along the -z direction.
            #self.connection.regions[1].offset_in_depth( -z_offset )
            mesh1new = TriMesh()
            mesh1new.faces = self.connection.regions[1].mesh.faces
            mesh1new.vs = array( self.connection.regions[1].mesh.vs )
            mesh1new.vs[:,2] += z_offset
            
            if kRealTimeHacks:
                ## Find the shape with the same uuid as regions0 and regions1 and translate
                ## one in depth.
                beauty = []
                for region in self.connection.regions:
                    beauty.append( [
                        shape for shape in self.annotations.anns_beautified
                        if hasattr( shape, 'uuid' ) and shape.uuid == region.uuid
                        ] )
                    ## There should be one shape with the same uuid.
                    assert len( beauty[-1] ) == 1
                    beauty[-1] = beauty[-1][0]
                
                ## Translate beauty1 in depth so that beauty0.z - beauty1.z == region0.z - region1.z.
                region0_z = average( asarray( self.connection.regions[0].mesh.vs )[:,2] )
                region1_z = average( mesh1new.vs[:,2] )
                
                beauty0_z = average( asarray( beauty[0].mesh.vs )[:,2] )
                beauty1_z = average( asarray( beauty[1].mesh.vs )[:,2] )
                
                beauty1new_z = beauty0_z + ( region1_z - region0_z )
                ## NOTE: offset_in_depth() offsets along the -z direction.
                beauty[1].offset_in_depth( beauty1_z - beauty1new_z )
                self.GetTopLevelParent().Refresh()
            
            ## 3
            ## ReplaceConnectionCurve() also clears .connection_curve_precise, so call it first.
            self.connection.ReplaceConnectionCurve( [world_pos] )
            from surface_intersections import mesh_mesh_intersection_lines_info
            intersection_info = mesh_mesh_intersection_lines_info(
                self.connection.regions[0].mesh,
                #self.connection.regions[1].mesh
                mesh1new
                )
            ## Ellipsoids generate the mesh on-the-fly, so let's access it once.
            face_normals0 = self.connection.regions[0].mesh.face_normals
            face_normals1 = self.connection.regions[1].mesh.face_normals
            ## Q: Do I need to do anything to this curve, such as cutting it to only use
            ##    edges on the right intersection panel (front/front or front/back, for
            ##    example) or simplifying it?
            ## A: Yes.  It confusing to see the front and back facing intersections at the same time.
            self.connection.connection_curve_precise = list(map( asarray, lines_to_line_strips( [
                segment[:,:2]
                for segment, faces in zip( intersection_info[0], intersection_info[1] )
                if (
                    ## Is the face front facing?
                    ## In other words, is the face normal's dot product with the look vector,
                    ## (0,0,-1), negative?
                    ## In other words, does the face normal point in the positive z direction?
                    bool( face_normals0[ faces[0] ][2] > 0. ) == bool( self.connection.regions_attaches_front[0] )
                    and
                    bool( face_normals1[ faces[1] ][2] > 0. ) == bool( self.connection.regions_attaches_front[1] )
                    )
                ] ) ))
            #print "connection_curve_precise:"
            #from pprint import pprint
            #pprint( self.connection.connection_curve_precise.round(3) )
    
    def doContextMenu( self, evt ):
        ## NOTE: I can't use evt.GetPosition() because it returns bogus data here.
        evt_pos = self.GetMousePosition()
        #print 'context menu:', self.event2screen( evt_pos )
        
        world_pos = self.event2world( evt_pos )
        shapes = [ ann for t,ann in PickShapeSurfaces( [ ann for ann in self.annotations.anns if ann not in self.locked_shapes ], append( world_pos, 0. ) ) ]
        #print 'context menu picked regions:', regions
        
        ## If we pick a single region, pop up a menu allowing the user to exclude or un-exclude it.
        if len( shapes ) == 1:
            menu = wx.Menu()
            menu_ids = []
            
            def exclude( evt ):
                self.excluded.add( shapes[0] )
            
            def unexclude( evt ):
                self.excluded.remove( shapes[0] )
            
            if shapes[0] not in self.excluded:
                item = menu.Append( wx.ID_ANY, "Exclude shape" )
                self.Bind( wx.EVT_MENU, exclude, id = item.GetId() )
            else:
                item = menu.Append( wx.ID_ANY, "Un-exclude shape" )
                self.Bind( wx.EVT_MENU, unexclude, id = item.GetId() )
            menu_ids.append( item.GetId() )
            
            self.PopupMenu( menu )
            menu.Destroy()
            ## I'm not sure where I got this idea of unbinding these events.
            ## I guess it's just good cleanup, so we don't leave a dangling bound callback to a menu
            ## event that can no longer be called, and also to free the id's we reserved so we
            ## don't run out of id space, since every time this is called we generate new ids.
            for id in menu_ids: self.Unbind( wx.EVT_MENU, id = id )
            
            return self
        
        ## Ignore regions in self.excluded.
        shapes = [ shape for shape in shapes if shape not in self.excluded ]
        
        ## If we don't pick exactly two regions, bail.
        if len( shapes ) != 2:
            evt.Skip()
            return
        
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
        
        def redraw( evt ):
            forward_events( self, self.ConnectionDrawer( self, shapes, connection ) )
        
        def drag( evt ):
            forward_events( self, self.ConnectionDragger( self, shapes, connection ) )
        
        item = menu.Append( wx.ID_ANY, "Draw Connection Curve" )
        self.Bind( wx.EVT_MENU, redraw, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        item = menu.Append( wx.ID_ANY, "Drag Connection Curve" )
        self.Bind( wx.EVT_MENU, drag, id = item.GetId() )
        menu_ids.append( item.GetId() )
        
        menu.AppendSeparator()
        
        def disconnect( evt ):
            self.SaveStateForUndo()
            RemoveAnnotation( self.annotations.anns, connection )
            
            self.annotations.Beautify()
            self.GetTopLevelParent().Refresh()
        
        item = menu.Append( wx.ID_ANY, "Disconnect" )
        self.Bind( wx.EVT_MENU, disconnect, id = item.GetId() )
        menu_ids.append( item.GetId() )
        menu.AppendSeparator()
        
        def connect_with_signature( signature ):
            self.SaveStateForUndo()
            ## TODO Q: I need to generate the precise intersection curve, if there is one.
            ##         Is there a clean way to do it?
            connection.regions_attaches_front = signature
            
            self.annotations.Beautify()
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
        ## I'm not sure where I got this idea of unbinding these events.
        ## I guess it's just good cleanup, so we don't leave a dangling bound callback to a menu
        ## event that can no longer be called, and also to free the id's we reserved so we
        ## don't run out of id space, since every time this is called we generate new ids.
        for id in menu_ids: self.Unbind( wx.EVT_MENU, id = id )
        
        return self
    
    def doMouseMove( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doDraw( self ):
        def draw_child():
            child = event_child( self )
            if child is not None: child.doDraw()
            else:
                self.prepare2D()
                
                world_pos = self.event2world( self.GetMousePosition() )
                ## If we don't pick exactly two regions under the mouse, don't draw a connection.
                shapes = [ ann for t,ann in PickShapeSurfaces( [ ann for ann in self.annotations.anns if ann not in self.locked_shapes ], append( world_pos, 0. ) ) ]
                ## Ignore regions in self.excluded.
                shapes = [ shape for shape in shapes if shape not in self.excluded ]
                if len( shapes ) == 2:
                    if self.kConnectionManipulationMode == 'draw':
                        draw_mode_label_text( 'Draw connection curve' )
                    elif self.kConnectionManipulationMode == 'drag':
                        draw_mode_label_text( 'Drag connection curve' )
                    elif self.kConnectionManipulationMode == 'draw-then-drag':
                        if ConnectionBetween( self.annotations.anns, shapes[0], shapes[1] ) is None:
                            draw_mode_label_text( 'Draw connection curve' )
                        else:
                            draw_mode_label_text( 'Drag connection curve' )
                    else:
                        raise RuntimeError("Unknown kConnectionManipulationMode: " + self.kConnectionManipulationMode)
                else:
                    draw_mode_label_text( 'Connect two shapes' )
                
                ## Draw a picking preview (smart guide).
                for shape in shapes:
                    if isinstance( shape, SweepShape ):
                        DrawPickedSweepShapeHandle( shape, (SweepShapePart.BoundaryInterior,) )
                    elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
                        ## TODO: Smart guides for ellipses.
                        pass
        
        def draw_callback():
            ## XXX TODO: This is ugly.
            if wx_annotate_shell.kPaperFigureMode:
                self.DrawAllGeometryAnnotations()
            else:
                self.DrawConnections()
            
            draw_child()
        
        self.DrawScene( callback = draw_callback )
#===================

#===================
class SamePropertyMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def GetPropertyName( self ):
        '''
        Returns the lower-case name of the property.
        '''
        
        raise NotImplementedError
    
    def PickManip( self, world_pos, world_dir = (0,0,-1) ):
        '''
        Given a world-space position 'world_pos' and optional world-space
        direction 'world_dir' (defaults to -z),
        returns a *Manip object or None for the picked manipulable object.
        '''
        
        raise NotImplementedError
    
    def DrawManipSmartGuide( self, manip ):
        '''
        Draws smart guide for the given manip.
        '''
        
        raise NotImplementedError
    
    def Manip2SameProperty( self, manip ):
        '''
        Returns the Same* object containing the given object 'manip' or None, if no
        Same* object contains 'manip'.
        '''
        
        raise NotImplementedError
    
    def NewSameProperty( self, manip ):
        '''
        Returns a new Same* object containing the given object 'manip'.
        '''
        
        raise NotImplementedError
    
    def doMouseDown( self, evt ):
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return
        
        ## From here, the initial state:
        ##   clicking on an unconstrained manipulable object starts a new group and we transition to assignment state
        ##   clicking on a constrained manipulable object selects that group and we transition to assignment state
        ##   right clicking a constrained manipulable object offers to remove it from its group
        ## From the assignment state:
        ##   clicking on a manipulable object assigns it to the assignment state's group
        ##   right clicking a constrained manipulable object bring up a contextual menu offering
        ##     to remove it from its group.  Actually, just use the delete key.
        
        ## Pick a manipulable object (circular cross section or symmetry plane).
        manip = self.PickManip( self.event2world( evt.GetPosition() ) )
        if manip is not None:
            self.annotations.MoveToFront( manip.shape )
            ## It is in a group?
            same_prop = self.Manip2SameProperty( manip )
            ## If it's not, start one.
            if same_prop is None:
                same_prop = self.NewSameProperty( manip )
                self.SaveStateForUndo()
                self.annotations.anns.append( same_prop )
            
            ## Switch to assignment state (mode)
            forward_events( self, self.SamePropertyAssigner( self, same_prop ) )
            ## We don't want to re-call self.doMouseDown( evt ), because the tilt_manip under
            ## the mouse is already in the group.
            ## We do want to refresh so that status text changes.
            self.Refresh()
    
    class SamePropertyAssigner( ModeProxy ):
        def __init__( self, canvas, same_prop ):
            ModeProxy.__init__( self, canvas )
            assert same_prop is not None
            self.same_prop = same_prop
        
        def doMouseDown( self, evt ):
            ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
            ## contextual menu event, but the right mouse button will be down.
            if evt.RightDown():
                print('skipping mouse down')
                evt.Skip()
                return
            
            ## Pick a manipulable object (circular cross section or symmetry plane).
            manip = self.PickManip( self.event2world( evt.GetPosition() ) )
            if manip is not None:
                self.annotations.MoveToFront( manip.shape )
                ## It is in a group?
                same_prop = self.Manip2SameProperty( manip )
                ## If it's already in our group, we're done.
                if same_prop is self.same_prop: return
                
                ## Otherwise, assign it to our group.
                self.SaveStateForUndo()
                self.same_prop.manips.append( manip )
                ## If it was in a different group, remove it from that group.
                ## NOTE: We can only garbage collect after adding it to our group,
                ##       because garbage collection might find our group with only one object in it
                ##       and remove it.
                if same_prop is not None:
                    same_prop.manips.remove( manip )
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                ## Regenerate results.
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            else:
                ## If we don't pick a tilt manipulable object, leave assignment mode.
                unforward_events( self )
            
            self.Refresh()
        
        def doKeyDown( self, evt ):
            ## Do the same thing as the main mode's keydown (delete tilt annotation under the mouse).
            event_parent( self )._doKeyDown( evt )
        
        def doMouseMove( self, evt ):
            ## Refresh so our status text changes.
            self.Refresh()
        
        def doDraw( self ):
            self.prepare2D()
            
            world_pos = self.event2world( self.GetMousePosition() )
            manip = self.PickManip( world_pos )
            if manip is not None:
                self.annotations.MoveToFront( manip.shape )
                self.DrawManipSmartGuide( manip )
                if self.Manip2SameProperty( manip ) is not self.same_prop:
                    draw_mode_label_text( 'Assign to same-' + self.GetPropertyName() + ' group' )
            else:
                draw_mode_label_text( 'De-select same-' + self.GetPropertyName() + ' group' )
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        if keycode in kDeleteKeys:
            
            ## evt.GetPosition() returns garbage.
            world_pos = self.event2world( self.GetMousePosition() )
            
            ## Are we over a manip?  If so, remove it from its group.
            manip = self.PickManip( world_pos )
            if manip is not None:
                self.annotations.MoveToFront( manip.shape )
                same_prop = self.Manip2SameProperty( manip )
                if same_prop is not None:
                    self.SaveStateForUndo()
                    same_prop.manips.remove( manip )
                    GarbageCollectMetaAnnotations( self.annotations.anns )
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
            
        else:
            evt.Skip()
    
    def doMouseMove( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doDraw( self ):
        def draw_child():
            child = event_child( self )
            if child is not None: child.doDraw()
            else:
                self.prepare2D()
                
                world_pos = self.event2world( self.GetMousePosition() )
                manip = self.PickManip( world_pos )
                if manip is not None:
                    self.annotations.MoveToFront( manip.shape )
                    self.DrawManipSmartGuide( manip )
                    if self.Manip2SameProperty( manip ) is None:
                        draw_mode_label_text( 'Start same-' + self.GetPropertyName() + ' group' )
                    else:
                        draw_mode_label_text( 'Choose same-' + self.GetPropertyName() + ' group' )
        
        def draw_callback():
            self.DrawAnnotations()
            draw_child()
        
        self.DrawScene( callback = draw_callback )

class SameTiltMode( SamePropertyMode ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def GetPropertyName( self ):
        return 'tilt'
    
    def PickManip( self, world_pos, world_dir = (0,0,-1) ):
        '''
        Given a world-space position 'world_pos' and optional world-space
        direction 'world_dir' (defaults to -z),
        returns a TiltManip object or None for the picked tilt manipulable object.
        '''
        
        ## Loop over every kind of shape annotation that has tilted circles, and see if we pick
        ## any of their tilted circles.
        for shape in self.annotations.anns:
            ## Skip locked shapes.
            if shape in self.locked_shapes: continue
            
            if isinstance( shape, SweepShape ):
                parts = [SweepShapePart.CrossSection] #, SweepShapePart.CrossSectionInterior]
                handle = PickSweepShapeHandle( shape, self, world_pos, world_dir, parts = parts )
                if handle is not None and handle[1] in shape.tilt_radians:
                    assert handle[0] in parts
                    return annotations_symbols.TiltManipSweepShape( shape, handle[1] )
            
            elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
                ## TODO: pick only the circular cross section
                handle = wx_annotate_sphere.PickAnnotationHandle(
                    shape,
                    world_pos,
                    lambda kind: PointDistanceThreshold( self )
                    )
                if handle is not None:
                    return annotations_symbols.TiltRadiusSizeManipEllipsoid( shape )
            
            else:
                ## Is it even a shape?  It could be any kind of annotation.
                #raise NotImplementedError, "No Manip for the shape type"
                pass
        
        return None
    
    def DrawManipSmartGuide( self, manip ):
        '''
        Draws smart guide for the given manip.
        '''
        
        if isinstance( manip.shape, SweepShape ):
            DrawPickedSweepShapeHandle( manip.shape, ( SweepShapePart.CrossSection, manip.t ) )
        elif isinstance( manip.shape, annotations_ellipsoid.AnnotationData ):
            pass
        else:
            raise RuntimeError("Unknown shape type")
    
    def Manip2SameProperty( self, tilt_manip ):
        '''
        Returns the annotations_symbols.SameTilt object containing the given object 'tilt_manip' or None, if no
        annotations_symbols.SameTilt object contains 'tilt_manip'.
        '''
        for same_tilt in self.annotations.anns.of_type( annotations_symbols.SameTilt ):
            if tilt_manip in same_tilt.manips:
                return same_tilt
        return None
    
    def NewSameProperty( self, manip ):
        '''
        Returns a new annotations_symbols.SameTilt object containing the given object 'manip'.
        '''
        return annotations_symbols.SameTilt( [ manip ] )
    
    def DrawAnnotations( self ):
        ## XXX TODO: This is ugly.
        if wx_annotate_shell.kPaperFigureMode:
            self.DrawAllGeometryAnnotations()
        else:
            self.DrawSameTiltAnnotations()

class SameSizeMode( SamePropertyMode ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def GetPropertyName( self ):
        return 'size'
    
    def PickManip( self, world_pos, world_dir = (0,0,-1) ):
        '''
        Given a world-space position 'world_pos' and optional world-space
        direction 'world_dir' (defaults to -z),
        returns a SizeManip object or None for the picked tilt manipulable object.
        '''
        
        ## Loop over every kind of shape annotation that has a curve, and see if we pick
        ## any of them.
        for shape in self.annotations.anns:
            ## Skip locked shapes.
            if shape in self.locked_shapes: continue
            
            if isinstance( shape, SweepShape ):
                handle = PickSweepShapeHandle( shape, self, world_pos, world_dir, parts = [SweepShapePart.Spine] )
                if handle is not None:
                    assert handle[0] == SweepShapePart.Spine
                    connections = [ con for con in self.annotations.anns.of_type( Connection ) if shape in con.regions ]
                    return annotations_symbols.SizeManipSweepShape( shape, connections[0] if len( connections ) > 0 else None )
            
            elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
                ## TODO: pick only the long axis
                handle = wx_annotate_sphere.PickAnnotationHandle(
                    shape,
                    world_pos,
                    lambda kind: PointDistanceThreshold( self )
                    )
                if handle is not None:
                    return annotations_symbols.TiltRadiusSizeManipEllipsoid( shape )
            
            else:
                ## Is it even a shape?  It could be any kind of annotation.
                #raise NotImplementedError, "No Manip for the shape type"
                pass
        
        return None
    
    def DrawManipSmartGuide( self, manip ):
        '''
        Draws smart guide for the given manip.
        '''
        
        if isinstance( manip.shape, SweepShape ):
            ## We don't have a specific 't' value along the spine; either re-pick the shape or use
            ## something arbitrary (0).
            DrawPickedSweepShapeHandle( manip.shape, ( SweepShapePart.Spine, 0. ) )
        elif isinstance( manip.shape, annotations_ellipsoid.AnnotationData ):
            pass
        else:
            raise RuntimeError("Unknown shape type")
    
    def Manip2SameProperty( self, manip ):
        '''
        Returns the annotations_symbols.SameSize object containing the given object 'manip' or None, if no
        annotations_symbols.SameSize object contains 'manip'.
        '''
        for same_size in self.annotations.anns.of_type( annotations_symbols.SameSize ):
            if manip in same_size.manips:
                return same_size
        return None
    
    def NewSameProperty( self, manip ):
        '''
        Returns a new annotations_symbols.SameSize object containing the given object 'manip'.
        '''
        return annotations_symbols.SameSize( [ manip ] )
    
    def DrawAnnotations( self ):
        ## XXX TODO: This is ugly.
        if wx_annotate_shell.kPaperFigureMode:
            self.DrawAllGeometryAnnotations()
        else:
            self.DrawSameSizeAnnotations()

class SameRadiusMode( SamePropertyMode ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def GetPropertyName( self ):
        return 'scale'
    
    def PickManip( self, world_pos, world_dir = (0,0,-1) ):
        '''
        Given a world-space position 'world_pos' and optional world-space
        direction 'world_dir' (defaults to -z),
        returns a RadiusManip object or None for the picked radius manipulable object.
        '''
        
        ## Loop over every kind of shape annotation that has tilted circles, and see if we pick
        ## any of their tilted circles.
        for shape in self.annotations.anns:
            ## Skip locked shapes.
            if shape in self.locked_shapes: continue
            
            if isinstance( shape, SweepShape ):
                parts = [SweepShapePart.Diameter]
                handle = PickSweepShapeHandle( shape, self, world_pos, world_dir, parts = parts )
                if handle is not None and handle[1] in shape.radii:
                    assert handle[0] in parts
                    return annotations_symbols.RadiusManipSweepShape( shape, handle[1] )
            
            elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
                ## TODO: pick only the circular cross section
                handle = wx_annotate_sphere.PickAnnotationHandle(
                    shape,
                    world_pos,
                    lambda kind: PointDistanceThreshold( self )
                    )
                if handle is not None:
                    return annotations_symbols.TiltRadiusSizeManipEllipsoid( shape )
            
            else:
                ## Is it even a shape?  It could be any kind of annotation.
                #raise NotImplementedError, "No Manip for the shape type"
                pass
        
        return None
    
    def DrawManipSmartGuide( self, manip ):
        '''
        Draws smart guide for the given manip.
        '''
        
        if isinstance( manip.shape, SweepShape ):
            ## We don't have a specific 't' value along the spine; either re-pick the shape or use
            ## something arbitrary (0).
            DrawPickedSweepShapeHandle( manip.shape, ( SweepShapePart.Diameter, manip.t ) )
        elif isinstance( manip.shape, annotations_ellipsoid.AnnotationData ):
            pass
        else:
            raise RuntimeError("Unknown shape type")
    
    def Manip2SameProperty( self, manip ):
        '''
        Returns the annotations_symbols.SameRadius object containing the given object 'manip' or None, if no
        annotations_symbols.SameRadius object contains 'manip'.
        '''
        for same_radius in self.annotations.anns.of_type( annotations_symbols.SameRadius ):
            if manip in same_radius.manips:
                return same_radius
        return None
    
    def NewSameProperty( self, manip ):
        '''
        Returns a new annotations_symbols.SameRadius object containing the given object 'manip'.
        '''
        return annotations_symbols.SameRadius( [ manip ] )
    
    def DrawAnnotations( self ):
        ## XXX TODO: This is ugly.
        if wx_annotate_shell.kPaperFigureMode:
            self.DrawAllGeometryAnnotations()
        else:
            self.DrawSameRadiusAnnotations()
#===================

#===================
class MirrorAcrossParentMode( ModeProxy ):
    def __init__( self, canvas ):
        ModeProxy.__init__( self, canvas )
    
    def PickManip( self, world_pos, world_dir = (0,0,-1) ):
        '''
        Given a world-space position 'world_pos' and optional world-space
        direction 'world_dir' (defaults to -z),
        returns a MirrorManipSweepShape object or None
        for the picked manipulable object.
        '''
        
        
        t_shapes = PickShapeSurfaces( [ ann for ann in self.annotations.anns if ann not in self.locked_shapes ], world_pos, world_dir )
        ## Choose shapes based on creation (list index) order, not 't' order.
        for t, shape in t_shapes:
            manip = self.ManipForShape( shape )
            if manip is not None:
                return manip
        
        return None
    
    def DrawManipSmartGuide( self, manip ):
        '''
        Draws smart guide for the given manip.
        '''
        
        if isinstance( manip.shape, SweepShape ):
            ## We don't have a specific 't' value along the spine; either re-pick the shape or use
            ## something arbitrary (0).
            DrawPickedSweepShapeHandle( manip.shape, ( SweepShapePart.Boundary, ) )
        elif isinstance( manip.shape, annotations_ellipsoid.AnnotationData ):
            pass
        else:
            raise RuntimeError("Unknown shape type")
    
    @staticmethod
    def ManipForShape( shape ):
        if isinstance( shape, SweepShape ):
            return annotations_symbols.MirrorManipSweepShape( shape )
        
        elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
            return annotations_symbols.MirrorManipEllipsoid( shape )
        
        else:
            ## Is it even a shape?  It could be any kind of annotation.
            #raise NotImplementedError, "No MirrorManip for the shape type"
            return None
    
    def doMouseDown( self, evt ):
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return
        
        ## From here, the initial state:
        ##   clicking on an unconstrained manipulable object starts a new group and we transition to assignment state
        ##   clicking on a constrained manipulable object selects that group and we transition to assignment state
        ##   right clicking a constrained manipulable object offers to remove it from its group
        ## From the assignment state:
        ##   clicking on a manipulable object assigns it to the assignment state's group
        ##   right clicking a constrained manipulable object bring up a contextual menu offering
        ##     to remove it from its group.  Actually, just use the delete key.
        
        ## Pick a manipulable object (circular cross section or symmetry plane).
        manip = self.PickManip( self.event2world( evt.GetPosition() ) )
        if manip is not None:
            self.annotations.MoveToFront( manip.shape )
            ## No matter what, we're going to change the annotations.
            self.SaveStateForUndo()
            
            ## If it's already the object-to-mirror in a mirror annotation, remove that annotation.
            ## Make a copy of the list of MirrorAcrossParent annotations, because we're removing
            ## as we iterate.
            for mirror_ann in list( self.annotations.anns.of_type( annotations_symbols.MirrorAcrossParent ) ):
                if mirror_ann.manip == manip:
                    mirror_ann.purge_mirrored_shape( self.annotations.anns )
                    RemoveAnnotation( self.annotations.anns, mirror_ann )
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
            
            ## Create a new mirror annotation.
            mirror_ann = annotations_symbols.MirrorAcrossParent( manip )
            
            ## Switch to parent assignment state (mode)
            forward_events( self, self.ParentAssigner( self, mirror_ann ) )
            ## We don't want to re-call self.doMouseDown( evt ), because the tilt_manip under
            ## the mouse is already in the group.
            ## We do want to refresh so that status text changes.
            self.Refresh()
    
    class ParentAssigner( ModeProxy ):
        def __init__( self, canvas, mirror_ann ):
            ModeProxy.__init__( self, canvas )
            assert mirror_ann is not None
            self.mirror_ann = mirror_ann
        
        def doMouseDown( self, evt ):
            ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
            ## contextual menu event, but the right mouse button will be down.
            if evt.RightDown():
                print('skipping mouse down')
                evt.Skip()
                return
            
            ## Pick a manipulable object (circular cross section or symmetry plane).
            manip = self.PickManip( self.event2world( evt.GetPosition() ) )
            if manip is not None: self.annotations.MoveToFront( manip.shape )
            if( manip is not None and
                manip != self.mirror_ann.manip and
                ConnectionBetween( self.annotations.anns, manip.shape, self.mirror_ann.manip.shape ) is not None
                ):
                ## Set it as the parent object of our mirror annotation.
                ## SaveStateForUndo() is already called before forwarding events to ParentAssigner.
                assert self.mirror_ann.parent_manip is None
                self.mirror_ann.parent_manip = manip
                self.annotations.anns.append( self.mirror_ann )
                ## Regenerate results.
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            ## No matter what, leave this mode.
            unforward_events( self )
            self.Refresh()
        
        def doMouseMove( self, evt ):
            ## Refresh so our status text changes.
            ## UPDATE: Actually, I changed doDraw() to always draw the same text.
            ## UPDATE 2: But I want mouse movement to move-to-front the shape under the mouse.
            self.Refresh()
            pass
        
        def doDraw( self ):
            self.prepare2D()
            
            ## Draw a small grey circle in the center of the to-be-aligned shape.
            #kSize = .05
            kSize = 2.5 * WorldDefaultLineWidth( self )
            col = array( ( .5, .5, .5 ) )
            
            self.glLineWidthWorld( 2*WorldDefaultLineWidth( self ) )
            
            glColor3fv( col )
            #draw_circle( self.mirror_ann.manip1.center()[:2], kSize )
            draw_circle_interior( self.mirror_ann.manip.center()[:2], kSize )
            glColor3fv( .5*col )
            draw_circle_interior( self.mirror_ann.manip.center()[:2], .5*kSize )
            
            ## Draw instructional text for the user.
            draw_mode_label_text( 'Choose mirror' )
            #return
            
            world_pos = self.event2world( self.GetMousePosition() )
            manip = self.PickManip( world_pos )
            if manip is not None:
                self.annotations.MoveToFront( manip.shape )
                self.DrawManipSmartGuide( manip )
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        if keycode in kDeleteKeys:
            
            ## evt.GetPosition() returns garbage.
            world_pos = self.event2world( self.GetMousePosition() )
            
            ## Are we over a manip in a mirror annotation?  If so, remove the mirror annotation.
            manip = self.PickManip( world_pos )
            if manip is not None:
                self.annotations.MoveToFront( manip.shape )
                for mirror_ann in self.annotations.anns.of_type( annotations_symbols.MirrorAcrossParent ):
                    if manip == mirror_ann.manip or manip == mirror_ann.parent_manip:
                        self.SaveStateForUndo()
                        mirror_ann.purge_mirrored_shape( self.annotations.anns )
                        RemoveAnnotation( self.annotations.anns, mirror_ann )
                        self.annotations.Beautify()
                        self.GetTopLevelParent().Refresh()
                        ## In case there are multiple, only remove the first one.
                        break
            
        else:
            evt.Skip()
    
    def doMouseMove( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doDraw( self ):
        def draw_child():
            child = event_child( self )
            if child is not None: child.doDraw()
            #else:
            elif True:
                self.prepare2D()
                
                world_pos = self.event2world( self.GetMousePosition() )
                manip = self.PickManip( world_pos )
                if manip is not None:
                    self.annotations.MoveToFront( manip.shape )
                    self.DrawManipSmartGuide( manip )
                    draw_mode_label_text( 'Mirror this shape about another' )
        
        def draw_callback():
            self.DrawAnnotations()
            draw_child()
        
        self.DrawScene( callback = draw_callback )
    
    def DrawAnnotations( self ):
        ## XXX TODO: This is ugly.
        if wx_annotate_shell.kPaperFigureMode:
            self.DrawAllGeometryAnnotations()
        else:
            self.DrawMirrorAcrossParentAnnotations()
#===================

#===================
class MirrorAlignAcrossParentMode( ModeProxy ):
    def __init__( self, canvas, pair_aligning_mode ):
        ModeProxy.__init__( self, canvas )
        
        self.pair_aligning_mode = pair_aligning_mode
    
    def PickManip( self, world_pos, world_dir = (0,0,-1) ):
        '''
        Given a world-space position 'world_pos' and optional world-space
        direction 'world_dir' (defaults to -z),
        returns a MirrorManipSweepShape object or None
        for the picked manipulable object.
        '''
        
        
        t_shapes = PickShapeSurfaces( [ ann for ann in self.annotations.anns if ann not in self.locked_shapes ], world_pos, world_dir )
        ## Choose shapes based on creation (list index) order, not 't' order.
        for t, shape in t_shapes:
            manip = self.ManipForShape( shape )
            if manip is not None:
                return manip
        
        return None
    
    def DrawManipSmartGuide( self, manip ):
        '''
        Draws smart guide for the given manip.
        '''
        
        if isinstance( manip.shape, SweepShape ):
            ## We don't have a specific 't' value along the spine; either re-pick the shape or use
            ## something arbitrary (0).
            DrawPickedSweepShapeHandle( manip.shape, ( SweepShapePart.Boundary, ) )
        elif isinstance( manip.shape, annotations_ellipsoid.AnnotationData ):
            pass
        else:
            raise RuntimeError("Unknown shape type")
    
    @staticmethod
    def ManipForShape( shape ):
        if isinstance( shape, SweepShape ):
            return annotations_symbols.MirrorManipSweepShape( shape )
        
        elif isinstance( shape, annotations_ellipsoid.AnnotationData ):
            return annotations_symbols.MirrorManipEllipsoid( shape )
        
        else:
            ## Is it even a shape?  It could be any kind of annotation.
            #raise NotImplementedError, "No MirrorManip for the shape type"
            return None
    
    def doMouseDown( self, evt ):
        ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
        ## contextual menu event, but the right mouse button will be down.
        if evt.RightDown():
            print('skipping mouse down')
            evt.Skip()
            return
        
        ## From here, the initial state:
        ##   clicking on an unconstrained manipulable object starts a new group and we transition to assignment state
        ##   clicking on a constrained manipulable object selects that group and we transition to assignment state
        ##   right clicking a constrained manipulable object offers to remove it from its group
        ## From the assignment state:
        ##   clicking on a manipulable object assigns it to the assignment state's group
        ##   right clicking a constrained manipulable object bring up a contextual menu offering
        ##     to remove it from its group.  Actually, just use the delete key.
        
        ## Pick a manipulable object (circular cross section or symmetry plane).
        manip = self.PickManip( self.event2world( evt.GetPosition() ) )
        if manip is not None:
            self.annotations.MoveToFront( manip.shape )
            ## No matter what, we're going to change the annotations.
            self.SaveStateForUndo()
            
            ## If it's already the object-to-mirror in a mirror annotation, remove that annotation.
            ## Make a copy of the list of MirrorAcrossParent annotations, because we're removing
            ## as we iterate.
            for mirror_ann in list( self.annotations.anns.of_type( annotations_symbols.MirrorAlignAcrossParent ) ):
                if manip in ( mirror_ann.manip1, mirror_ann.manip2 ):
                    RemoveAnnotation( self.annotations.anns, mirror_ann )
                    self.annotations.Beautify()
                    self.GetTopLevelParent().Refresh()
            
            ## Create a new mirror annotation.
            mirror_ann = annotations_symbols.MirrorAlignAcrossParent( manip )
            
            ## Switch to parent assignment state (mode)
            forward_events( self, self.SecondAssigner( self, mirror_ann ) )
            ## We don't want to re-call self.doMouseDown( evt ), because the tilt_manip under
            ## the mouse is already in the group.
            ## We do want to refresh so that status text changes.
            self.Refresh()
    
    class SecondAssigner( ModeProxy ):
        def __init__( self, canvas, mirror_ann ):
            ModeProxy.__init__( self, canvas )
            assert mirror_ann is not None
            assert mirror_ann.manip1 is not None
            assert mirror_ann.manip2 is None
            
            self.mirror_ann = mirror_ann
        
        def doMouseDown( self, evt ):
            ## This is crucial for contextual menus to work.  Mouse down gets first dibs on a potential
            ## contextual menu event, but the right mouse button will be down.
            if evt.RightDown():
                print('skipping mouse down')
                evt.Skip()
                return
            
            ## Whether to regenerate results.
            did_something = False
            
            ## Pick a manipulable object (one with a symmetry plane).
            manip = self.PickManip( self.event2world( evt.GetPosition() ) )
            if manip is not None: self.annotations.MoveToFront( manip.shape )
            if( manip is not None and
                manip != self.mirror_ann.manip1
                ):
                if not self.pair_aligning_mode:
                    ## If manip is directly connected to mirror_ann.manip1, then
                    ## mirror-align manip1 on manip as the symmetry parent.
                    if ConnectionBetween( self.annotations.anns, manip.shape, self.mirror_ann.manip1.shape ) is not None:
                        self.mirror_ann.manip2 = self.mirror_ann.manip1
                        self.mirror_ann.parent_manip = manip
                        did_something = True
                ## Otherwise, if there is a bridging shape between them, use the bridging shape as
                ## the symmetry parent.
                else:
                    assert self.pair_aligning_mode
                    ## We can only use 'manip' if it is one node away from mirror_ann.manip1.
                    ## First we must search for such a bridging (parent) node.
                    ## Collect all nodes connected to manip and self.mirror_ann.manip1.
                    
                    manip_neighbors = neighbors( self.annotations.anns, manip.shape )
                    manip1_neighbors = neighbors( self.annotations.anns, self.mirror_ann.manip1.shape )
                    
                    if self.mirror_ann.manip1.shape not in manip_neighbors and manip.shape not in manip1_neighbors:
                        set_intersection = manip_neighbors.intersection( manip1_neighbors )
                        assert len( set_intersection ) <= 1
                        if len( set_intersection ) == 1:
                            bridge_shape = list(set_intersection)[0]
                            
                            ## Set manip as the second shape and the bridging shape as the parent object
                            ## of our mirror annotation.
                            assert self.mirror_ann.manip2 is None
                            self.mirror_ann.manip2 = manip
                            assert self.mirror_ann.parent_manip is None
                            self.mirror_ann.parent_manip = self.ManipForShape( bridge_shape )
                            did_something = True
            
            if did_something:
                ## SaveStateForUndo() is already called before forwarding events to SecondAssigner.
                self.annotations.anns.append( self.mirror_ann )
                ## Regenerate results.
                self.annotations.Beautify()
                self.GetTopLevelParent().Refresh()
            
            ## No matter what, leave this mode.
            unforward_events( self )
            self.Refresh()
        
        def doMouseMove( self, evt ):
            ## Refresh so our status text changes.
            ## UPDATE: Actually, I changed doDraw() to always draw the same text.
            ## UPDATE 2: Actually, I changed it back because of smart guides.
            self.Refresh()
        
        def doDraw( self ):
            self.prepare2D()
            
            ## Draw a small grey circle in the center of the to-be-aligned shape.
            #kSize = .05
            kSize = 2.5 * WorldDefaultLineWidth( self )
            col = array( ( .5, .5, .5 ) )
            
            self.glLineWidthWorld( 2*WorldDefaultLineWidth( self ) )
            
            glColor3fv( col )
            #draw_circle( self.mirror_ann.manip1.center()[:2], kSize )
            draw_circle_interior( self.mirror_ann.manip1.center()[:2], kSize )
            glColor3fv( .5*col )
            draw_circle_interior( self.mirror_ann.manip1.center()[:2], .5*kSize )
            
            ## Draw instructional text for the user.
            draw_mode_label_text( 'Choose second shape' )
            
            world_pos = self.event2world( self.GetMousePosition() )
            manip = self.PickManip( world_pos )
            if manip is not None:
                self.annotations.MoveToFront( manip.shape )
                self.DrawManipSmartGuide( manip )
    
    def doKeyDown( self, evt ):
        #print evt.GetKeyCode()
        keycode = evt.GetKeyCode()
        
        if keycode in kDeleteKeys:
            
            ## evt.GetPosition() returns garbage.
            world_pos = self.event2world( self.GetMousePosition() )
            
            ## Are we over a manip in a mirror annotation?  If so, remove the mirror annotation.
            manip = self.PickManip( world_pos )
            if manip is not None:
                self.annotations.MoveToFront( manip.shape )
                for mirror_ann in self.annotations.anns.of_type( annotations_symbols.MirrorAlignAcrossParent ):
                    if manip in ( mirror_ann.manip1, mirror_ann.manip2, mirror_ann.parent_manip ):
                        self.SaveStateForUndo()
                        RemoveAnnotation( self.annotations.anns, mirror_ann )
                        self.annotations.Beautify()
                        self.GetTopLevelParent().Refresh()
                        ## In case there are multiple, only remove the first one.
                        break
            
        else:
            evt.Skip()
    
    def doMouseMove( self, evt ):
        ## Refresh so our status text changes.
        self.Refresh()
    
    def doDraw( self ):
        def draw_child():
            child = event_child( self )
            if child is not None: child.doDraw()
            #else:
            elif True:
                self.prepare2D()
                
                world_pos = self.event2world( self.GetMousePosition() )
                manip = self.PickManip( world_pos )
                if manip is not None:
                    self.annotations.MoveToFront( manip.shape )
                    self.DrawManipSmartGuide( manip )
                    if self.pair_aligning_mode:
                        draw_mode_label_text( 'Align two shapes about another' )
                    else:
                        draw_mode_label_text( 'Center one shape on another' )
        
        def draw_callback():
            self.DrawAnnotations()
            draw_child()
        
        self.DrawScene( callback = draw_callback )
    
    def DrawAnnotations( self ):
        ## XXX TODO: This is ugly.
        if wx_annotate_shell.kPaperFigureMode:
            self.DrawAllGeometryAnnotations()
        else:
            self.DrawMirrorAlignAcrossParentAnnotations()

class MirrorAlignAcrossParentModePaired( MirrorAlignAcrossParentMode ):
    def __init__( self, canvas ):
        MirrorAlignAcrossParentMode.__init__( self, canvas, True )
class MirrorAlignAcrossParentModeSingle( MirrorAlignAcrossParentMode ):
    def __init__( self, canvas ):
        MirrorAlignAcrossParentMode.__init__( self, canvas, False )
#===================

#===================
class AnnotateCanvasCollection( AnnotateCanvas ):
    def __init__( self, *args, **kwargs ):
        AnnotateCanvas.__init__( self, *args, **kwargs )
        
        import illustrator_mode
        import cuboid.ui
        #mode_id = self.RegisterMode( illustrator_mode.ArrowMode, "Select and manipulate shapes", "Select" )
        mode_id = self.RegisterMode( illustrator_mode.ArrowMode2, "Select and manipulate shapes", "Select" )
        self.RegisterMode( illustrator_mode.DrawSweepShapeMode, "Create a new cylindrical shape with a spine", "Spine" )
        self.RegisterMode( cuboid.ui.DrawCuboidMode, "Create a new cuboid shape from a curve", "Cuboid" )
        self.RegisterMode( illustrator_mode.DrawEllipsoidShapeMode, "Create a new ellipsoid", "Ellipse" )
        
        #mode_id = self.RegisterMode( SweepShapeAnnotationCreateMode, "Create a new cylindrical shape with a spine", "Spine" )
        #self.RegisterMode( EllipseAnnotationCreateMode, "Create a new ellipsoid", "Ellipse" )
        #self.RegisterMode( SymmetryXSectionAnnotationCreateMode, "Create a new rotund shape with a symmetry cross section", "SymmetryXSection" )
        
        self.RegisterMode( ConnectionCurveMode, "Draw connection curves", "Connect" )
        self.RegisterMode( SameTiltMode, "Mark same tilt angles", "Same Tilt" )
        self.RegisterMode( SameSizeMode, "Mark same length shapes", "Same Length" )
        self.RegisterMode( SameRadiusMode, "Mark same radii", "Same Radius" )
        self.RegisterMode( MirrorAcrossParentMode, "Mirror shapes across their parents", "Mirror Shape" )
        self.RegisterMode( MirrorAlignAcrossParentModePaired, "Align shapes across their parent", "Align Shape" )
        self.RegisterMode( MirrorAlignAcrossParentModeSingle, "Center a shape on its parent", "Center Shape" )
        
        self.RegisterMode( illustrator_mode.ZoomMode, "Zoom in or out", "Zoom" )
        self.RegisterMode( illustrator_mode.PanMode, "Pan the view", "Pan" )
        
        self.SetModeByID( mode_id )
        
        self.shape_fill = False
        
        self.draw_shape_cache = dict()
        
        self.locked_shapes = set()

    def PickShapeHandlesByType( self, screen_pos, shape_type ):
        print("TODO Q: What's the right convenience method?")
        for shape in self.annotations.anns.of_type( SymmetryXSectionShape ):
            pass
    
    def DrawScene( self, **kwargs ):
        # print( "AnnotateCanvasCollection.DrawScene()" )
        # import traceback
        # traceback.print_stack( limit = 5 )
        
        kwargs.setdefault( 'callback', None )
        kwargs.setdefault( 'shape_annotations', True )
        kwargs.setdefault( 'background', True )
        ## A flag for whether-or-not to draw handles for shape annotations.
        ## TODO Q: How should I organize my drawing code?  Should these two default to False?
        kwargs.setdefault( 'shape_annotation_handles', False )
        kwargs.setdefault( 'geometry_annotations', False ) ## Change this flag for the other paper figure mode
        ## TODO: Flags for all the things different modes would want drawn (glow under mouse).
        
        #glClearColor( .9,.9,.9,0 )
        glClearColor( 1,1,1,1 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        self.prepare2D()
        
        ## Draw the gltex if there is one.
        if kwargs['background'] and self.annotations.gltex( self.context ) is not None:
            glColor4f( 1,1,1, .5 )
            draw_gltexture( self.annotations.gltex( self.context ), self.annotations.tex_pos, self.annotations.tex_size )
        
        if self.show_annotations:
            ## Draw the solid shape annotations
            if kwargs['shape_annotations']:
                ## TODO: Obey kwargs['shape_annotation_handles']
                for shape in self.annotations.anns[::-1]:
                    ## Ellipsoid shapes
                    if isinstance( shape, annotations_ellipsoid.AnnotationData ):
                        if shape in self.locked_shapes:
                            self.DrawEllipseLocked( shape, **kwargs )
                        else:
                            self.DrawEllipse( shape, **kwargs )
                    
                    ## symmetry_xsection shapes
                    elif isinstance( shape, SymmetryXSectionShape ):
                        if shape in self.locked_shapes:
                            raise NotImplementedError
                        else:
                            self.DrawSymmetryXSectionShape( shape, **kwargs )
                    
                    ## sweep shapes
                    elif isinstance( shape, SweepShape ):
                        if shape in self.locked_shapes:
                            self.DrawSweepShapeLocked( shape, **kwargs )
                        else:
                            self.DrawSweepShape( shape, **kwargs )
            
            ## Draw the geometry and connection annotations
            if kwargs['geometry_annotations']:
                self.DrawAllGeometryAnnotations()
        
        if kwargs['callback'] is not None: kwargs['callback']()
        
        self.SwapBuffers()
    
    def DrawEllipse( self, ellipse_ann, **kwargs ):
        if wx_annotate_shell.kPaperFigureMode:
            self.glLineWidthWorld( WorldDefaultLineWidth( self ) )
            glColor3ub( *Color["silhouette"] )
            draw_ellipse( *ellipse_ann.ellipse )
            return
        
        if self.shape_fill:
            glColor4fv( ub2f( Color["shape_fill"] ) )
            draw_ellipse_interior( *ellipse_ann.ellipse )
        
        self.glLineWidthWorld( WorldDefaultLineWidth( self ) )
        glColor3ub( *Color["silhouette"] )
        draw_ellipse( *ellipse_ann.ellipse )
        
        self.glLineWidthWorld( .75*WorldDefaultLineWidth( self ) )
        
        lineloop = line_loop_from_tilted_circle( ellipse_ann.ellipsoid.center, mag( ellipse_ann.ellipsoid.axes[1] ), ellipse_ann.ellipsoid.axes[0] )
        
        #set_linear_alpha_for_points( [z for x,y,z in lineloop], ub2f( Color["circular_xsection"] ) )
        ##glColor3fv( (213/255., 99/255., 191/255.) )
        #draw_lineloop( lineloop )
        #unset_linear_alpha_for_points()
        
        colors = colors_from_linear_alpha_for_points( [z for x,y,z in lineloop], ub2f( Color["circular_xsection"] ) )
        draw_thick_lineloop2d( lineloop, .75*WorldDefaultLineWidth( self ), colors )
    
    def DrawEllipseLocked( self, ellipse_ann, **kwargs ):
        glLineWidth( 1 )
        glColor3ub( *Color["silhouette"] )
        draw_ellipse( *ellipse_ann.ellipse )
    
    def DrawEllipseHandles( self, which = None ):
        self.prepare2D()
        
        #self.glLineWidthWorld( .75*WorldDefaultLineWidth( self ) )
        #glColor3f( .75, .75, .75 )
        glColor3f( 1,1,1 )
        for ann in self.annotations.anns.of_type( annotations_ellipsoid.AnnotationData ):
            if which is not None and ann not in which: continue
            
            for handle in ann.Handles:
                ## Don't draw the whole ellipse handle.
                if handle.kind == annotations_ellipsoid.AnnotationData.WholeEllipseHandle:
                    pass
                ## TODO: Either don't draw the unpickable, obscured axis 0 (long axis) handle,
                ##       or make it pickable.
                elif handle.kind == annotations_ellipsoid.AnnotationData.Axis0Handle:
                    sign = 1. if ann.ellipsoid.axes[0][2] >= 0. else -1.
                    draw_widget(
                        'ellipse_rotscale',
                        ann.ellipsoid.center + sign*ann.ellipsoid.axes[0],
                        sign*ann.ellipsoid.axes[0],
                        8*PointDistanceThreshold( self )
                        )
                    if abs( ann.ellipsoid.axes[0][2] ) < kEps:
                        draw_widget( 
                            'ellipse_rotscale',
                            ann.ellipsoid.center - sign*ann.ellipsoid.axes[0],
                            -sign*ann.ellipsoid.axes[0],
                            8*PointDistanceThreshold( self )
                            )
                elif handle.kind == annotations_ellipsoid.AnnotationData.Axis1Handle:
                    draw_widget(
                        'radius_manip',
                        ann.ellipsoid.center + ann.ellipsoid.axes[1],
                        ann.ellipsoid.axes[1],
                        6*PointDistanceThreshold( self )
                        )
                    draw_widget(
                        'radius_manip',
                        ann.ellipsoid.center - ann.ellipsoid.axes[1],
                        -ann.ellipsoid.axes[1],
                        6*PointDistanceThreshold( self )
                        )
                elif handle.kind == annotations_ellipsoid.AnnotationData.CrossSectionHandle:
                    draw_widget(
                        'cross_section_tilt',
                        ann.ellipsoid.center + ann.ellipsoid.axes[2],
                        ann.ellipsoid.axes[0],
                        4*PointDistanceThreshold( self )
                        )
                else:
                    pass
                    #draw_circle( handle.get_pos(), PointDistanceThreshold( self ) )
    
    def DrawSymmetryXSectionShape( self, sym_ann, **kwargs ):
        self.glLineWidthWorld( WorldDefaultLineWidth( self ) )
        
        if sym_ann.silhouette is not None:
            glColor3ub( *Color["silhouette"] )
            draw_lineloop( sym_ann.silhouette )
        
        for circular_xsection_3d_line_loop in sym_ann.circular_xsections_3d_line_loop:
            set_linear_alpha_for_points( [z for x,y,z in circular_xsection_3d_line_loop], ub2f( Color["circular_xsection"] ) )
            #glColor3fv( (213/255., 99/255., 191/255.) )
            draw_lineloop( circular_xsection_3d_line_loop )
            unset_linear_alpha_for_points()
        
        if sym_ann.symmetry_xsection_3d is not None:
            set_linear_alpha_for_points( [z for x,y,z in sym_ann.symmetry_xsection_3d], ub2f( Color["symmetry_xsection"] ) )
            draw_lineloop( sym_ann.symmetry_xsection_3d )
            unset_linear_alpha_for_points()
        elif sym_ann.symmetry_xsection is not None:
            glColor3ub( *Color["symmetry_xsection"] )
            draw_lineloop( sym_ann.symmetry_xsection )
    
    def DrawSweepShape( self, shape, **kwargs ):
        
        ## Make a weak keyed dictionary from shapes to display lists.
        ## To check if the shape changed, we will use the hash used to generate the mesh.
        ## UPDATE: Translation doesn't invalidate any hashes, so we add the first point
        ##         of the spine to the hash.
        ## UPDATE 2: Zoom affects how we draw, so we need to hash it, too.
        import weakref
        
        ## NOTE: This ignores kwargs!
        if shape not in [ key() for key in list(self.draw_shape_cache.keys()) ]:
            ## Make the first item a new object(), so it will always fail comparison with anything
            ## else, causing us to recompute the display list in the next if statement.
            value = [ ( object() ), glGenLists( 1 ) ]
            self.draw_shape_cache[ weakref.ref( shape, lambda x: glDeleteLists( value[1], 1 ) ) ] = value
        else:
            value = self.draw_shape_cache[ weakref.ref( shape ) ]
        
        cur_hash = ( shape.gen_mesh_hash(), tuple( shape.spine2d[0] ), self.GetZoom(), self.shape_fill, wx_annotate_shell.kPaperFigureMode )
        if value[0] == cur_hash:
            glCallList( value[1] )
        else:
            glNewList( value[1], GL_COMPILE_AND_EXECUTE )
            self._DrawSweepShape( shape, **kwargs )
            glEndList()
            value[0] = cur_hash
        
        ## For debugging:
        if hasattr( shape, 'debug_draw' ):
            for data in shape.debug_draw:
                color = ( 255, 0, 0 )
                
                kind = data[0]
                pts = data[1]
                if len( data ) > 2: color = data[2]
                
                glColor3ub( *color )
                if kind == 'lineloop':
                    draw_thick_lineloop2d( pts, WorldDefaultLineWidth( self ) )
                elif kind == 'linestrip':
                    draw_thick_linestrip2d( pts, WorldDefaultLineWidth( self ) )
                else:
                    print('ERROR: unknown type')
    
    def _DrawSweepShape( self, shape, **kwargs ):
        if wx_annotate_shell.kPaperFigureMode:
            glColor3ub( *Color["silhouette"] )
            draw_thick_lineloop2d( shape.boundary2d, WorldDefaultLineWidth( self ) )
            return
        
        ## Alpha under the shape.
        if self.shape_fill: #True: #not wx_annotate_shell.kPaperFigureMode:
            from gluTesselated import gluTesselated
            glColor4fv( ub2f( Color["shape_fill"] ) )
            gluTesselated( shape.boundary2d ).draw()
        
        
        if shape.spine2d is not None:
            glColor3ub( *Color["spine"] )
            draw_thick_linestrip2d( shape.spine2d, WorldDefaultLineWidth( self ) )
        
        ## Comment out this next block if making paper figures.
        if not wx_annotate_shell.kPaperFigureMode:
            glColor3ub( *Color["diameter"] )
            for t in shape.radii.keys():
                center = shape.spine2d_at_t( t )
                normal = dir( shape.normal2d_at_t( t ) )
                radius = shape.radius_at_t( t )
                linestrip = [ center + radius * normal, center - radius * normal ]
                
                draw_thick_linestrip2d( linestrip, .25*WorldDefaultLineWidth( self ) )
        
        self.glLineWidthWorld( .75*WorldDefaultLineWidth( self ) )
        for t in chain( iter(shape.tilt_radians.keys()), iter(shape.cross_sections.keys()) ):
            lineloop = shape.cross_section_line_loop3d_at_t( t )
            
            if t in shape.cross_sections:
                color = ub2f( Color["symmetry_xsection"] )
            else:
                color = ub2f( Color["circular_xsection"] )
            
            #set_linear_alpha_for_points( [z for x,y,z in lineloop], color )
            ##glColor3fv( (213/255., 99/255., 191/255.) )
            #draw_lineloop( lineloop )
            #unset_linear_alpha_for_points()
            
            colors = colors_from_linear_alpha_for_points( [z for x,y,z in lineloop], color )
            draw_thick_lineloop2d( lineloop, .75*WorldDefaultLineWidth( self ), colors )
        
        glColor3ub( *Color["silhouette"] )
        draw_thick_lineloop2d( shape.boundary2d, WorldDefaultLineWidth( self ) )
        
        
        if wx_annotate_shell.kDebugging:
            self.prepare2D()
            
            glMatrixMode( GL_MODELVIEW )
            glPushMatrix()
            glTranslatef( shape.spine2d[0][0], shape.spine2d[0][1], 0. )
            scale = .2 * 1./72.
            glScalef( scale, scale, scale )
            
            glColor3ub( 41, 60, 127 )
            glEnable( GL_TEXTURE_2D )
            font.Render( str( self.annotations.anns.index( shape ) ) )
            glDisable( GL_TEXTURE_2D )
            
            glPopMatrix()
    
    def DrawSweepShapeLocked( self, shape, **kwargs ):
        glColor3ub( *Color["silhouette"] )
        glLineWidth( 1 )
        draw_lineloop( shape.boundary2d )
    
    def DrawSweepShapeHandles( self, which = None ):
        self.prepare2D()
        for shape in self.annotations.anns.of_type( SweepShape ):
            if which is not None and shape not in which: continue
            
            ## Draw symmetry plane and its manipulation handle
            kThickness = .75*WorldDefaultLineWidth( self )
            self.glLineWidthWorld( kThickness )
            glColor3ub( *Color['symmetry_sheet'] )
            
            rays = shape.symmetry_plane_rays_at_vertices()
            rays = [ ( asarrayf( ray[0] ), dir( asarrayf( ray[1] ) ) ) for ray in rays ]
            
            ## Draw a "square" running along the top and bottom of the shape, as well as the
            ## sides (end-caps).
            ## Top and bottom.
            #draw_linestrip( [ ray[0] + radius * ray[1] for ray, radius in izip( rays, shape.dense_radii ) ] )
            draw_thick_linestrip2d( [ ray[0] + radius * ray[1] for ray, radius in zip( rays, shape.dense_radii ) ], kThickness )
            #draw_linestrip( [ ray[0] - radius * ray[1] for ray, radius in izip( rays, shape.dense_radii ) ] )
            draw_thick_linestrip2d( [ ray[0] - radius * ray[1] for ray, radius in zip( rays, shape.dense_radii ) ], kThickness )
            ## Sides (at end-caps).
            glBegin( GL_LINES )
            for ray, radius in [ ( rays[0], shape.dense_radii[0] ), ( rays[-1], shape.dense_radii[-1] ) ]:
                glVertex3fv( ray[0] + radius * ray[1] )
                glVertex3fv( ray[0] - radius * ray[1] )
            glEnd()
            
            ## Draw the symmetry plane handles.
            for t in shape.symmetry_radians.keys():
                sym_handle_i = shape.vertex_index_at_t( t )
                sym_handle_pos = rays[ sym_handle_i ][0] + shape.dense_radii[ sym_handle_i ] * rays[ sym_handle_i ][1]
                #draw_circle( sym_handle_pos[:2], PointDistanceThreshold( self ) )
                tangent = shape.tangent3d_at_t( t )
                glColor3f( 1,1,1 )
                draw_widget( 'symmetry_plane_manip', sym_handle_pos, ( -tangent[1], tangent[0] ), 4*PointDistanceThreshold( self ) )
            
            ## Draw the cap radius handles.
            for t in (0,1):
                pos = shape.spine2d_at_t( t )
                tangent = shape.tangent2d_at_t( t )
                cap_radius = shape.cap_radius_at_t( t )
                
                ## Is the mouse near the point cap_radius away from the start or end point?
                cap_handle_pos = pos + cap_radius * ( -1. if t == 0 else 1 ) * dir( tangent )
                
                alpha = .5 if shape.cap_radii[t] is None else 1.
                glColor4f( 1,1,1, alpha )
                draw_widget( 'cap_radius_manip', cap_handle_pos, tangent, 6*PointDistanceThreshold( self ) )
            
            ## Finally, draw the radius and tilting handles on top.
            ## Tilting handles.
            for t, tilt_radians in shape.tilt_radians.items():
                center2d = asarray( shape.spine2d_at_t( t ) )
                radius = shape.radius_at_t( t )
                ## I can't use the normal because the tilt handle is being manipulated by the user,
                ## but the 3d spine isn't being updated in real-time.
                ## UPDATE: I can; these accessors check to see if the value was specifically overridden.
                #tangent3d = annotations_sweep.tilt_radians2tangent3d( tilt_radians, shape.tangent2d_at_t( t ) )
                tangent3d = shape.tangent3d_at_t( t )
                
                tilt_dir = dir( cross( cross( tangent3d, (0,0,1) ), tangent3d ) )
                
                ## Draw radius manipulation handles.
                self.glLineWidthWorld( .75*WorldDefaultLineWidth( self ) )
                #glColor3f( .75, .75, .75 )
                glColor3f( 1,1,1 )
                
                ## Draw circle tilting handle.
                draw_widget(
                    'cross_section_tilt',
                    append( center2d, 0. ) + radius * tilt_dir,
                    tangent3d,
                    4*PointDistanceThreshold( self )
                    )
            ## Radius handles.
            for t in shape.radii.keys():
                center2d = asarray( shape.spine2d_at_t( t ) )
                ## I can't use the radius_at_t() accessor in case the radius handle is being
                ## manipulated by the user, because the radius it returns isn't being updated in
                ## real-time.
                ## UPDATE: I can; these accessors now check to see if the value was specifically overridden.
                radius = shape.radius_at_t( t )
                normal = shape.tangent3d_at_t( t )
                
                ## Draw radius manipulation handles.
                self.glLineWidthWorld( .75*WorldDefaultLineWidth( self ) )
                #glColor3f( .75, .75, .75 )
                glColor3f( 1,1,1 )
                
                perp = zerosf(2)
                ## 90-degree right turn
                perp[0] = normal[1]
                perp[1] = -normal[0]
                perp = dir( perp )
                perp *= radius
                
                #draw_circle( center2d + perp, PointDistanceThreshold( self ) )
                #draw_circle( center2d - perp, PointDistanceThreshold( self ) )
                ## Multiply PointDistanceThreshold( self ) by two, because this widget is larger
                ## than the pickable circle.
                draw_widget( 'radius_manip', center2d + perp, perp, 6*PointDistanceThreshold( self ) )
                draw_widget( 'radius_manip', center2d - perp, -perp, 6*PointDistanceThreshold( self ) )
    
    def DrawAllGeometryAnnotations( self ):
        self.DrawMirrorAlignAcrossParentAnnotations()
        self.DrawMirrorAcrossParentAnnotations()
        self.DrawSameRadiusAnnotations()
        self.DrawSameSizeAnnotations()
        self.DrawSameTiltAnnotations()
        self.DrawConnections()
    
    def DrawConnections( self ):
        self.prepare2D()
        for connection in self.annotations.anns.of_type( Connection ):
            ## Don't draw connections on locked shapes.
            if len( [ r for r in connection.regions if r in self.locked_shapes] ) > 0: continue
            
            ## Darkening ranges from 0 to 1, depending on which surface we are attaching to.  (1 for front/front, 0 for back/back.)
            darkening = ( 2*connection.regions_attaches_front[0] + connection.regions_attaches_front[1] ) / 3.
            #darkening = .5**( 3*(1. - darkening) )
            #darkening = .5*darkening + .5
            
            ## It looks like we can't (any longer --- we used to with 3.0.0b3) pass a tuple
            ## containing floating point values to this PyOpenGL function; instead, we have
            ## to pass an array.
            glColor4ubv( asarray( tuple( darkening * asarray( Color[ 'connection' ] )) + (200,) ) )
            if connection.connection_curve_precise is None:
                draw_thick_linestrip2d( connection.connection_curve, 2*WorldDefaultLineWidth( self ) )
            else:
                for segment in connection.connection_curve_precise:
                    draw_thick_linestrip2d( segment, 2*WorldDefaultLineWidth( self ) )
    
    def GetSamePropertyColor( self, property ):
        max_id = 10000
        ## Pick some nice, per-same property group color.
        col = color_spiral( property.unique_value * max_id, max_id )
        ## Darken it somewhat.
        col = .75 * asarrayf( col )
        return col
    def SetSamePropertyColor( self, property ):
        glColor3fv( self.GetSamePropertyColor( property ) )
    
    def DrawSameTiltAnnotations( self ):
        self.prepare2D()
        for same_tilt in self.annotations.anns.of_type( annotations_symbols.SameTilt ):
            
            #kSize = .04
            kSize = 2. * WorldDefaultLineWidth( self )
            
            self.SetSamePropertyColor( same_tilt )
            
            for tilt_manip in same_tilt.manips:
                ## Don't draw annotations on locked shapes.
                if tilt_manip.shape in self.locked_shapes: continue
                
                
                ## Draw an arc attached to a symmetric L shape, like in the mode icon.
                
                center, normal = list(map( asarrayf, tilt_manip.get_center_normal() ))
                normal = dir( normal )
                draw_widget( 'same_angle', center, normal, 4*kSize )
                
                '''
                right = zerosf( 2 )
                right[0] = normal[1]
                right[1] = -normal[0]
                
                
                ## The L shape.
                self.glLineWidthWorld( .75*WorldDefaultLineWidth( self ) )
                glBegin( GL_LINE_STRIP )
                
                glVertex2fv( center + kSize * normal )
                glVertex2fv( center )
                glVertex2fv( center + kSize * right )
                
                glEnd()
                
                ## The arc.
                self.glLineWidthWorld( WorldDefaultLineWidth( self ) )
                draw_ellipse(
                    center, kSize*right, kSize*normal,
                    NUM_SLICES = 5, start_radians = 0, end_radians = pi/2
                    )
                '''
    
    def DrawSameSizeAnnotations( self ):
        self.prepare2D()
        for same_size in self.annotations.anns.of_type( annotations_symbols.SameSize ):
            
            ## Make sure there is enough space between the lines, even if they get very thick.
            #kSize = max( .06, 5*WorldDefaultLineWidth( self ) )
            #kSize = max( 3*WorldDefaultLineWidth( self ), 5*WorldDefaultLineWidth( self ) )
            kSize = 5*kWorldDefaultLineWidth
            
            self.SetSamePropertyColor( same_size )
            
            for manip in same_size.manips:
                ## Don't draw annotations on locked shapes.
                if manip.shape in self.locked_shapes: continue
                
                
                ## Draw two perpendicular tick marks (short line segments).
                
                length = manip.get_size()
                center, normal = list(map( asarrayf, manip.get_center_normal() ))
                normal = dir( normal )
                
                right = zerosf( 2 )
                right[0] = normal[1]
                right[1] = -normal[0]
                
                size = 3*length * kSize
                thickness = 1.5*length * kWorldDefaultLineWidth
                
                
                ## Two tick marks.
                tick1 = (
                    center - .4*size * right - .2*size * normal,
                    center + .4*size * right - .2*size * normal
                    )
                draw_thick_linestrip2d( tick1, 2*thickness )
                tick2 = (
                    center - .4*size * right + .2*size * normal,
                    center + .4*size * right + .2*size * normal
                    )
                draw_thick_linestrip2d( tick2, 2*thickness )
    
    def DrawSameRadiusAnnotations( self ):
        self.prepare2D()
        for same_radius in self.annotations.anns.of_type( annotations_symbols.SameRadius ):
            
            ## Make sure there is enough space between the lines, even if they get very thick.
            #kSize = .04
            kSize = max( .04, 2.5*kWorldDefaultLineWidth )
            
            self.SetSamePropertyColor( same_radius )
            
            for radius_manip in same_radius.manips:
                ## Don't draw annotations on locked shapes.
                if radius_manip.shape in self.locked_shapes: continue
                
                
                ## Draw a line with arrows at both ends.
                
                radius = radius_manip.get_radius()
                center, normal = list(map( asarrayf, radius_manip.get_center_normal() ))
                normal = dir( normal )
                
                right = zerosf( 2 )
                right[0] = normal[1]
                right[1] = -normal[0]
                
                size = 10.*radius * kSize
                thickness = 5.*radius * kWorldDefaultLineWidth
                
                
                ## The length of the line.
                line = (
                    center + (radius-.5*size) * right,
                    center - (radius-.5*size) * right
                    )
                draw_thick_linestrip2d( line, 1.5*thickness )
                
                ## An arrow at one end
                arrow1 = (
                    line[0] - size*right + size*normal,
                    line[0],
                    line[0] - size*right - size*normal
                    )
                draw_thick_linestrip2d( arrow1, 1.5*thickness )
                
                ## An arrow at the other end
                arrow2 = (
                    line[1] + size*right + size*normal,
                    line[1],
                    line[1] + size*right - size*normal
                    )
                draw_thick_linestrip2d( arrow2, 1.5*thickness )
                
                ## Two 'same size' tick marks halfway along the arrow.
                tick1 = (
                    center - .7*size * normal - .2*size * right,
                    center + .7*size * normal - .2*size * right
                    )
                draw_thick_linestrip2d( tick1, 1*thickness )
                
                tick2 = (
                    center - .7*size * normal + .2*size * right,
                    center + .7*size * normal + .2*size * right
                    )
                draw_thick_linestrip2d( tick2, 1*thickness )
    
    def DrawMirrorAcrossParentAnnotations( self ):
        self.prepare3D()
        
        ## For each mirror annotation, draw an equals sign with a tilde above it at the
        ## child's origin.
        for mirror_ann in self.annotations.anns.of_type( annotations_symbols.MirrorAcrossParent ):
            assert mirror_ann.manip is not None and mirror_ann.parent_manip is not None
            ## Don't draw annotations on locked shapes.
            ## Actually, only don't draw it if the child is locked.
            #if( mirror_ann.manip.shape in self.locked_shapes or
            #    mirror_ann.parent_manip.shape in self.locked_shapes ): continue
            if mirror_ann.manip.shape in self.locked_shapes: continue
            
            self.SetSamePropertyColor( mirror_ann )
            
            origin3d = mirror_ann.manip.origin( mirror_ann.parent_manip.shape, self.annotations.anns )
            
            glMatrixMode( GL_MODELVIEW )
            glPushMatrix()
            glTranslatef( *origin3d )
            scale = .4 * 1./72.
            glScalef( scale, scale, scale )
            
            glEnable( GL_TEXTURE_2D )
            font.Render( '=' )
            ## font.BBox() returns 6 scalars.
            ## The first three are the 3d coordinates of the lower-left point.
            ## The next three are the 3d coordinates of the upper-right point.
            equals_height = font.BBox( '=' )[4] - font.BBox( '=' )[1]
            glTranslatef( 0, equals_height, 0 )
            font.Render( '~' )
            glDisable( GL_TEXTURE_2D )
            
            glPopMatrix()
    
    def DrawMirrorAlignAcrossParentAnnotations( self ):
        self.prepare2D()
        
        #kSize = 0.03
        kSize = 1.5 * WorldDefaultLineWidth( self )
        
        ## For each mirrored-alignment annotation, draw two small circles at the connection
        ## points connected by a line.
        for mirror_ann in self.annotations.anns.of_type( annotations_symbols.MirrorAlignAcrossParent ):
            assert mirror_ann.manip1 is not None
            assert mirror_ann.manip2 is not None
            assert mirror_ann.parent_manip is not None
            
            ## Don't draw annotations on locked shapes.
            if( mirror_ann.manip1.shape in self.locked_shapes or
                mirror_ann.manip2.shape in self.locked_shapes ): continue
            
            col = self.GetSamePropertyColor( mirror_ann )
            
            origin1_3d = mirror_ann.manip1.origin( mirror_ann.parent_manip.shape, self.annotations.anns )
            origin2_3d = mirror_ann.manip2.origin( mirror_ann.parent_manip.shape, self.annotations.anns )
            
            self.glLineWidthWorld( 1.5*WorldDefaultLineWidth( self ) )
            #draw_circle( origin1_3d, kSize )
            #draw_circle( origin2_3d, kSize )
            glColor3fv( col )
            draw_circle_interior( origin1_3d, kSize + .5*1.5*WorldDefaultLineWidth( self ) )
            draw_circle_interior( origin2_3d, kSize + .5*1.5*WorldDefaultLineWidth( self ) )
            
            self.glLineWidthWorld( 2*WorldDefaultLineWidth( self ) )
            glBegin( GL_LINES )
            glVertex3fv( origin1_3d )
            glVertex3fv( origin2_3d )
            glEnd()
            
            glColor3fv( .5*col )
            draw_circle_interior( origin1_3d, kSize )
            draw_circle_interior( origin2_3d, kSize )

class ResultCanvasCollection( ResultCanvas ):
    def __init__( self, *args, **kwargs ):
        ResultCanvas.__init__( self, *args, **kwargs )
        
        self.npr = True
        self.shader = None
        self.silhouettes = True
    
    def EnableNPRShader( self ):
        if self.shader is None:
            try:
                self.shader = ShaderProgram()
                self.shader.LoadVertexShaderFromString(
                    '''
                    varying vec3 normal;
                    void main() {
                        normal = gl_NormalMatrix * gl_Normal;
                        gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
                    }
                    '''
                    )
                self.shader.LoadFragmentShaderFromString(
                    '''
                    varying vec3 normal;
                    void main() {
                        float intensity;
                        vec4 color;
                        vec3 n = normalize(normal);
                        // TODO Q: Why doesn't this work?  This is what was here by default.
                        // vec3 l = normalize(gl_LightSource[0].position).xyz;
                        // A headlight:
                        vec3 l = vec3( 0,0,1 );
                        // A light coming from the +y direction in world-space (tied to the shape):
                        // vec3 l = normalize( ( gl_ModelViewMatrix * vec4( 0,0,1,0 ) ).xyz );
                        
                        // quantize to 5 steps (0, .25, .5, .75 and 1)
                        intensity = (floor(dot(l, n) * 4.0) + 1.0)/4.0;
                        intensity = max( 0., min( 1., intensity ) );
                        
                        // This gives us the standard headlight look.
                        // intensity = n.z;
                        
                        // Visualize the normal.
                        // color = vec4(n[0], n[1], n[2], 1.0);
                        // Visualize the light direction.
                        // color = vec4( l, 1.0 );
                        
                        // Pass the intensity straight through.
                        // color = vec4( intensity, intensity, intensity, 1.0 );
                        
                        // This doesn't work, probably because lighting is enabled.
                        // color = intensity * gl_Color;
                        // color[3] = 1.;
                        
                        // This works (lighting is enabled).
                        // color = intensity * gl_FrontMaterial.diffuse;
                        // color[3] = 1.;
                        
                        // Pass the color through unchanged (won't work if lighting is enabled).
                        // color = gl_Color;
                        
                        // This is what was originally in the shader.
                        // color = vec4(intensity*1.0, intensity*0.5, intensity*0.5, intensity*1.0);
                        // A "bisc" color.
                        color = vec4(intensity*0.95, intensity*0.92, intensity*0.75, 1.0);
                        
                        gl_FragColor = color;
                    }
                    '''
                    )
            
            except NotImplementedError:
                ## The NoShaderProgram class has no-ops for UseProgram() and UnUseProgram() methods.
                self.shader = NoShaderProgram()
        
        self.shader.UseProgram()
    
    def DisableNPRShader( self ):
        self.shader.UnUseProgram()
    
    def doDraw( self ):
        #glClearColor( .5,.5,.5,0 )
        #glClearColor( 1,1,1,0 )
        glClearColor( 1,1,1,1 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        self.prepare3D()
        
        ## cartoon shader
        if self.npr: self.EnableNPRShader()
        
        #glColor3ub( 200, 200, 255 )
        glColor3ub( 200, 200, 200 )
        
        def draw_thing_silhouette( draw_command ):
            kSilhouetteThickness = 3
            glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
            ## Multiply the silhouette thickness by two because half of the line is
            ## occluded by the front faces.
            glLineWidth( kSilhouetteThickness * 2 )
            glDisable( GL_LIGHTING )
            glColor3ub( 0,0,0 )
            glEnable( GL_CULL_FACE )
            glCullFace( GL_FRONT )
            
            if self.npr: self.DisableNPRShader()
            draw_command()
            
            glPolygonMode( GL_FRONT_AND_BACK, GL_FILL )
            glEnable( GL_LIGHTING )
            glColor3ub( 200, 200, 200 )
            glEnable( GL_CULL_FACE )
            glCullFace( GL_BACK )
            
            if self.npr: self.EnableNPRShader()
            draw_command()
        
        def draw_thing_ordinary( draw_command ):
            draw_command()
        
        draw_thing = draw_thing_ordinary
        if self.silhouettes: draw_thing = draw_thing_silhouette
        
        ## Draw only the resulting meshes and shapes.
        if self.show_meshes:
            for shape in self.annotations.anns_beautified:
                ## Ellipsoid shapes
                if isinstance( shape, annotations_ellipsoid.AnnotationData ):
                    #pass
                    #draw_mesh_faces_smooth( shape.mesh, True )
                    
                    #shape.ellipsoid.Draw()
                    draw_thing( lambda: shape.ellipsoid.Draw() )
                
                ## rotund and symmetry_xsection shapes
                elif hasattr( shape, "mesh" ) and shape.mesh is not None:
                    #draw_mesh_faces_smooth( shape.mesh, True )
                    draw_thing( lambda: draw_mesh_faces_smooth( shape.mesh, True ) )
                    # if kHackDrawEdges:
                    #     draw_thing( lambda: [ draw_mesh_faces_smooth( shape.mesh, True ), glLineWidth( 1 ), glColor3ub( 0,0,0 ), glDisable( GL_LIGHTING ), draw_mesh_edges( shape.mesh ), glEnable( GL_LIGHTING ) ] )
        
        if self.npr: self.DisableNPRShader()
        
        if wx_annotate_shell.kDebugging:
            DrawDebugging( self.prepare2D, self.prepare3D, self.annotations.anns_beautified )
        
        
        ## Draw the gltex if there is one.
        if self.annotations.gltex( self.context ) is not None:
            self.drawTexture(
                self.annotations.gltex( self.context ),
                self.annotations.tex_pos,
                self.annotations.tex_size
                )
        
        self.SwapBuffers()

def main():
    import sys
    argv = list( sys.argv[1:] )
    parse_options( argv )
    
    import wx
    app = wx.App()
    
    ## load the PyCrust shell
    ## NOTE: This generates the following errors and is very slow switching between the app window
    ##       and the PyCrust window.
    ## Sun Jun 15 20:19:36 yotronagon.local Python[32825] <Error>: CGContextRestoreGState: invalid context
    if False: #wx_annotate_shell.kDebugging:
        from wxBasicGLCanvas import load_PyCrust
        crust = load_PyCrust( app )
    
    annotations = AnnotationsBeautify()
    if False and wx_annotate_shell.kDebugging:
        
        ## Add a symmetry xsection shape
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
        
        
        ## Add an ellipse shape
        shape = annotations_ellipsoid.AnnotationData()
        shape.SetCircle( [
            ( .5 + .4 * cos( theta ), .5 + .5*sin( theta ) )
            for theta in linspace( 0, 2*pi, 10 )[:-1]
            ] )
        annotations.anns.append( shape )
        
        
        annotations.Beautify()
        #debugger()
    
    frame = AnnotateFrame(
        None, annotations,
        annotate_canvas_factory = AnnotateCanvasCollection,
        result_canvas_factory = ResultCanvasCollection
        )
    frame.Show(True)
    
    import os
    for fname in argv:
        ext = os.path.splitext( fname )[1].lower()
        if ext in ( '.annscene', ):
            annotations.Load( fname )
            frame.save_filename = fname
            frame.Refresh()
        elif ext in ( '.png', '.jpg', '.tiff', '.tif', '.gif', '.bmp' ):
            annotations.SetImage( fname )
            frame.Refresh()
    
    app.annotateframe = frame
    app.MainLoop()
    app.Destroy()

if __name__ == '__main__': main()
