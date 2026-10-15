from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *
from helpers import *
from vecutil import *
kEps = 1e-5

'''
Same tilt -> give two tilted circles the same tilt angle
Parallel plane -> give two planes or tilted circles the same tilt angle and the same axis of rotation
    same axis of rotation means
        averaging axes of rotation (as unit vectors)
        and then,
        in the case of sweep surfaces, deforming the 2d spine curves a la laplacian editing so it changes as little as possible but has those gradients implied by axis and tilt
Same length -> finds average length of 3d curves (spines or ellipse axes), then scales shapes so that curves have that length
Same radii -> averages all radii
Symmetric -> all connected objects (that don't intersect the symmetry plane in 3d?) must have same length and same radii and same tilt and optionally (with a toggle) same axis reflected about symmetry plane

ApplyConstraints( anns ) -> Modifies anns to apply constraints.  (Callee should worry about making a copy if needed.)
                            Can I build a linear system
                            (with DOFs: tilt angle, axis of rotation, length, and radii)
                            to solve all of these constraints at once?
'''

def ApplyConstraints( anns ):
    '''
    Modifies sequence of annotations 'anns' by applying constraints defined in this module.
    '''
    
    for same_tilt in anns.of_type( SameTilt ):
        same_tilt.apply()
    
    for same_radius in anns.of_type( SameRadius ):
        same_radius.apply()
    
    ## TODO Q: Should this happen first?  I would need to change update_3d() to respect z offsets.
    ## A1: It would make sense since the user drew the connection curve on the un-beautified shape.
    from wx_annotate_symmetry_xsection import PlaceInDepth
    try:
        PlaceInDepth( anns )
    except AssertionError:
        pass
    
    for same_size in anns.of_type( SameSize ):
        same_size.apply( anns )
    
    for mirror_ann in anns.of_type( MirrorAlignAcrossParent ):
        mirror_ann.apply( anns )
    
    ## Copy the list of mirror annotations in anns, because mirror_anns' apply() modifies anns itself.
    for mirror_ann in list( anns.of_type( MirrorAcrossParent ) ):
        mirror_ann.apply( anns )


class SameProperty( object ):
    def __init__( self, manips = None, master_manips = None ):
        if manips is None: manips = []
        if master_manips is None: master_manips = []
        
        ## 'manips' stores a list of objects with methods required by subclasses.
        self.manips = manips
        
        ## 'master_manips' stores a list of manips (that must also be present in manips).
        ## These manips will not be affected by other manipulators.
        ## If there is only one, it will determine the size for all manipulators in manips.
        ## If there is more than one, their average will be propagated to manipulators in manips.
        self.master_manips = list( master_manips )
        
        from random import Random
        ## unique_value is used to determine the color of this annotation, and nothing else.
        self.unique_value = Random().uniform(0.,1.)
        
        self.assert_self()
    
    def __setstate__( self, unpickled ):
        self.__dict__.update( unpickled )
        
        ## To allow unpickling from older versions of the class
        if not hasattr( self, 'master_manips' ): self.master_manips = []
        if not hasattr( self, 'unique_value' ):
            from random import Random
            self.unique_value = Random().uniform(0.,1.)
    
    def assert_self( self ):
        assert False not in [ manip in self.manips for manip in self.master_manips ]
    
    def garbage_collect( self, anns ):
        '''
        Removes manips that have become invalid.
        '''
        remove_me = [ manip for manip in self.manips if manip.is_garbage( anns ) ]
        for manip in remove_me: self.manips.remove( manip )
        
        remove_me = [ manip for manip in self.master_manips if manip.is_garbage( anns ) ]
        for manip in remove_me: self.master_manips.remove( manip )
        
        self.assert_self()
    
    def is_garbage( self, anns ):
        '''
        Returns whether or not the entire object has become invalid.
        '''
        return len( self.manips ) < 2
    
    def apply( self ):
        raise NotImplementedError


class TiltManipSweepShape( object ):
    def __init__( self, shape, t ):
        ## Keep a reference to which tilted cross section.
        self.t = t
        assert self.t >= 0.
        assert self.t <= 1.
        ## Keep a reference to the shape so we can
        ## convert to center/normal representation and
        ## update the 3d mesh and
        ## check if the tilted circle has been removed.
        self.shape = shape
    
    def __eq__( self, other ):
        if type( self ) is not type( other ): return NotImplemented
        
        return(
                ## The same shape
                self.shape is other.shape and
                ## The circle is attached to the same location
                self.t == other.t
                )
    
    def __ne__( self, other ):
        return not self == other
    
    def __setstate__( self, unpickled ):
        ## Convert older-version objects.
        if 'tilted_circle' in unpickled:
            from annotations_sweep import edge_index2t
            unpickled['t'] = edge_index2t( unpickled['tilted_circle'][0], unpickled['shape'].spine2d )
            del unpickled['tilted_circle']
        
        assert 'tilted_circle' not in unpickled
        assert 't' in unpickled
        self.__dict__.update( unpickled )
    
    def get_tilt( self ):
        return self.shape.tilt_radians[ self.t ]
    
    def set_tilt( self, new_tilt ):
        self.shape.tilt_radians[ self.t ] = new_tilt
        ## TODO Q: Do I need to call update_3d() after every type of constraint is applied,
        ##         or can I defer calling it until later.  If I can defer, I may want
        ##         some way to call update_3d automatically, on-demand.
        self.shape.update_3d()
    
    def get_center_normal( self ):
        '''
        Returns the 2d center and normal direction for this annotation.
        This is useful for drawing the annotation.
        '''
        center = self.shape.spine2d_at_t( self.t )
        normal = self.shape.tangent2d_at_t( self.t )
        return center, dir( normal )
    
    def is_garbage( self, anns ):
        ## Is our shape or tilted circle gone?
        return self.shape not in anns or self.t not in self.shape.tilt_radians

class SameTilt( SameProperty ):
    '''
    'self.manips' objects must have methods
    get_tilt(), which returns the tilt, and
    set_tilt( new_tilt ), which sets the tilt, and
    get_center_normal(), which returns a tuple of the 2d center and 2d normal direction.
    Note that get_center_normal() is used for UI purposes only.
    '''
    
    def apply( self ):
        '''
        Average the absolute value of all tilt radians, and set them all to the average value,
        keeping the original sign.
        '''
        
        ## Sweep shapes and ellipsoids have an orientation, so using a signed tilt angle doesn't do what the user wants.
        ## As a result, we'll use absolute values and preserve signs where we can.
        
        self.assert_self()
        
        if len( self.master_manips ) > 0:
            manips_to_average = self.master_manips
        else:
            manips_to_average = self.manips
        
        average_tilt = 0.
        for tilt_manip in manips_to_average:
            average_tilt += abs( tilt_manip.get_tilt() )
        
        if len( manips_to_average ) > 0:
            average_tilt /= len( manips_to_average )
        
        for tilt_manip in self.manips:
            tilt_manip.set_tilt( average_tilt * ( -1. if tilt_manip.get_tilt() < 0 else 1. ) )


class SizeManipSweepShape( object ):
    def __init__( self, shape, parenting_connection = None ):
        self.shape = shape
        ## parenting_connection is the connection that determines the origin of scaling.
        self.parenting_connection = parenting_connection
    
    def __eq__( self, other ):
        if type( self ) is not type( other ): return NotImplemented
        
        return(
                ## The same shape
                self.shape is other.shape
                ## We don't check if parenting_connection is the same, because only one
                ## size constraint is allowed per shape.
                )
    
    def __ne__( self, other ):
        return not self == other
    
    def get_origin( self ):
        if self.parenting_connection is None:
            return self.shape.center()
        else:
            from annotations_rotund import connection_frozen_com
            parent = list(self.parenting_connection.regions)
            parent.remove( self.shape )
            ## Use the connection curve's center-of-mass as the origin.
            #return connection_frozen_com( self.parenting_connection, parent )
            
            ## Use the endpoint closest to the connection curve's center-of-mass.
            connection_com = connection_frozen_com( self.parenting_connection, parent )
            '''
            print "used candidate origin", argmin( [
                mag2( pt[:2] - connection_com[:2] )
                for pt in self.shape.candidate_origins()
                ] )
            '''
            '''
            print "dist2candidate_origins:"
            from pprint import pprint
            pprint( [
                ( mag2( pt[:2] - connection_com[:2] ), tuple( pt ) )
                for pt in self.shape.candidate_origins()
                ] )
            '''
            ## I need to convert to a tuple and back to an array because min() will compare the (incomparable) arrays if the distances are the same.
            return asarray( min( [
                ( mag2( pt[:2] - connection_com[:2] ), tuple( pt ) )
                for pt in self.shape.candidate_origins()
                ] )[1] )
    
    def get_size( self ):
        return self.shape.length3d()
    
    def set_size( self, new_size, anns ):
        assert new_size > 0.
        
        scale = new_size / self.get_size()
        assert scale > 0.
        
        origin = self.get_origin()
        
        ## First, translate all attached shapes, not crossing parent_connection.
        ## We must do this first because they will compare to our geometry to find the origin.
        def translate_shape( con, shape, offset ):
            if con is not None: con.ReplaceConnectionCurve( asarray( con.connection_curve ) + offset[:2] )
            shape.translate( offset )
        
        from annotations_rotund import Connection
        from wx_annotate_collection import ConnectionBetween, run_dfs
        
        ## Collect all direct children (and the incoming connection),
        ## and the list of edges excluding those connecting self.shape to its children.
        connection_children = []
        cons_to_remove = []
        cons = anns.of_type( Connection )
        for con in cons:
            if self.shape in con.regions:
                other = set( con.regions )
                other.remove( self.shape )
                assert len( other ) == 1
                other = list(other)[0]
                
                cons_to_remove.append( con )
                ## Don't cross the parenting connection.
                if con is not self.parenting_connection:
                    connection_children.append( ( con, other ) )
        
        for con in cons_to_remove: cons.remove( con )
        
        ## Initiate a dfs for every child node as root, with offset given by the new position
        ## of the root's origin after being transformed by the scale with given origin.
        for connection_to_child, child in connection_children:
            ## Don't offset in z; we haven't placed in depth, so perform 2d scaling and offsetting.
            child_origin = SizeManipSweepShape( child, connection_to_child ).get_origin()
            offset = scale * ( asarrayf( child_origin ) - origin ) + origin - child_origin
            
            ## Since we are seeding each child separately, manually translate the connection to child.
            connection_to_child.ReplaceConnectionCurve( asarray( connection_to_child.connection_curve ) + offset[:2] )
            
            run_dfs( cons, pre_visit = lambda c,n: translate_shape( c, n, offset ), roots = [child] )
        
        
        ## Now transform self.shape.
        self.shape.set_length3d( new_size, origin )
    
    def get_center_normal( self ):
        '''
        Returns the 2d center and normal direction for this annotation.
        This is useful for drawing the annotation.
        '''
        
        ## Place it in the middle of the spine.
        halfway = len( self.shape.spine2d ) // 2
        center = self.shape.spine2d[ halfway ]
        normal = asarrayf( self.shape.spine2d[ halfway+1 ] ) - center
        return center, dir( normal )
    
    def is_garbage( self, anns ):
        return self.shape not in anns

class TiltRadiusSizeManipEllipsoid( SizeManipSweepShape ):
    def __init__( self, *args, **kwargs ):
        SizeManipSweepShape.__init__( self, *args, **kwargs )
    
    def get_tilt( self ):
        return self.shape.GetTiltRadians()
    
    def set_tilt( self, new_tilt ):
        self.shape.SetTiltRadians( new_tilt )
    
    def get_radius( self ):
        return self.shape.GetRadius()
    
    def set_radius( self, radius ):
        self.shape.SetRadius( radius )
    
    def get_center_normal( self ):
        '''
        Returns the 2d center and normal direction for this annotation.
        This is useful for drawing the annotation.
        '''
        return self.shape.ellipse[0], dir( self.shape.ellipse[1] )

class SameSize( SameProperty ):
    '''
    'self.manips' objects must have methods
    get_size(), which returns the size, and
    set_size( new_size, origin ), which sets the size, and
    get_center_normal(), which returns a tuple of the 2d center and 2d normal direction.
    Note that get_center_normal() is used for UI purposes only.
    '''
    
    def apply( self, anns ):
        '''
        Average all sizes, and set them all to the average value.
        '''
        
        self.assert_self()
        
        if len( self.master_manips ) > 0:
            manips_to_average = self.master_manips
        else:
            manips_to_average = self.manips
        
        average_size = 0.
        for size_manip in manips_to_average:
            average_size += abs( size_manip.get_size() )
        
        if len( manips_to_average ) > 0:
            average_size /= len( manips_to_average )
        
        for size_manip in self.manips:
            size_manip.set_size( average_size, anns )

class RadiusManipSweepShape( object ):
    def __init__( self, shape, t ):
        ## Keep a reference to which radius.
        self.t = t
        assert self.t >= 0.
        assert self.t <= 1.
        ## Keep a reference to the shape so we can
        ## convert to center/normal representation and
        ## update the 3d mesh and
        ## check if the tilted circle has been removed.
        self.shape = shape
    
    def __eq__( self, other ):
        if type( self ) is not type( other ): return NotImplemented
        
        return(
                ## The same shape
                self.shape is other.shape and
                ## The circle is attached to the same location
                self.t == other.t
                )
    
    def __ne__( self, other ):
        return not self == other
    
    def __setstate__( self, unpickled ):
        ## Convert older-version objects.
        if 'tilted_circle' in unpickled:
            from annotations_sweep import edge_index2t
            unpickled['t'] = edge_index2t( unpickled['tilted_circle'][0], unpickled['shape'].spine2d )
            del unpickled['tilted_circle']
        
        assert 'tilted_circle' not in unpickled
        assert 't' in unpickled
        self.__dict__.update( unpickled )
    
    def get_radius( self ):
        return self.shape.radii[ self.t ]
    
    def set_radius( self, new_radius ):
        self.shape.radii[ self.t ] = new_radius
        ## TODO Q: Do I need to call update_3d() after every type of constraint is applied,
        ##         or can I defer calling it until later.  If I can defer, I may want
        ##         some way to call update_3d automatically, on-demand.
        self.shape.update_3d()
    
    def get_center_normal( self ):
        '''
        Returns the 2d center and normal direction for this annotation.
        This is useful for drawing the annotation.
        '''
        center = self.shape.spine2d_at_t( self.t )
        normal = self.shape.tangent2d_at_t( self.t )
        return center, dir( normal )
    
    def is_garbage( self, anns ):
        ## Is our shape or tilted circle gone?
        return self.shape not in anns or self.t not in self.shape.radii

class SameRadius( SameProperty ):
    '''
    'self.manips' objects must have methods
    get_radius(), which returns the radius, and
    set_radius( new_radius ), which sets the radius, and
    get_center_normal(), which returns a tuple of the 2d center and 2d normal direction.
    Note that get_center_normal() is used for UI purposes only.
    '''
    
    def apply( self ):
        '''
        Average the absolute value of all radii, and set them all to the average value,
        keeping the original sign.
        '''
        
        self.assert_self()
        
        if len( self.master_manips ) > 0:
            manips_to_average = self.master_manips
        else:
            manips_to_average = self.manips
        
        average_radius = 0.
        for radius_manip in manips_to_average:
            average_radius += radius_manip.get_radius()
        
        if len( manips_to_average ) > 0:
            average_radius /= len( manips_to_average )
        
        for radius_manip in self.manips:
            radius_manip.set_radius( average_radius )


class MirrorManipSweepShape( object ):
    '''
    This class implements both interfaces, manip and parent_manip, for passing to MirrorAcrossParent and MirrorAlignAcrossParent.
    '''
    
    def __init__( self, shape ):
        self.shape = shape
    
    def __eq__( self, other ):
        if type( self ) is not type( other ): return NotImplemented
        
        return(
                ## The same shape
                self.shape is other.shape
                )
    
    def __ne__( self, other ):
        return not self == other
    
    def is_garbage( self, anns ):
        '''
        Returns True if this object is no longer relevant to the annotations 'anns'.
        '''
        return self.shape not in anns
    
    def meshcopy( self ):
        '''
        Returns a copy of this object as a mesh_gui.TriMesh.
        '''
        from mesh_gui import TriMesh
        result = TriMesh()
        result.vs = asarray( self.shape.mesh.vs ).copy()
        result.faces = asarray( self.shape.mesh.faces ).copy()
        return result
    
    def meshcopy_recursive( self, parent_shape, anns ):
        '''
        Returns a copy of this object and its entire subtree
        (all connected objects without crossing 'parent_shape')
        as a single mesh_gui.TriMesh (with possibly multiple connected components).
        '''
        
        from mesh_gui import TriMesh
        result = TriMesh()
        
        def paste_mesh( con, shape ):
            result.append( shape.mesh )
        
        from annotations_rotund import Connection
        from wx_annotate_collection import ConnectionBetween, run_dfs
        
        cons = anns.of_type( Connection )
        cons.remove( ConnectionBetween( anns, self.shape, parent_shape ) )
        run_dfs( cons, pre_visit = paste_mesh, roots = [self.shape] )
        
        return result
    
    def translate_recursive( self, parent_shape, anns, offset ):
        '''
        In-place translation of this object and its entire subtree
        (all connected objects without crossing 'parent_shape')
        by 3d vector 'offset'.
        
        NOTE: This also translates connections in the subtree.
        '''
        
        def translate_shape( con, shape ):
            if con is not None: con.ReplaceConnectionCurve( asarray( con.connection_curve ) + offset[:2] )
            shape.translate( offset )
        
        from annotations_rotund import Connection
        from wx_annotate_collection import ConnectionBetween, run_dfs
        
        cons = anns.of_type( Connection )
        cons.remove( ConnectionBetween( anns, self.shape, parent_shape ) )
        run_dfs( cons, pre_visit = translate_shape, roots = [self.shape] )
    
    def subtree_shapes( self, parent_shape, anns ):
        '''
        Returns a list containing this object and its entire subtree
        (all connected objects without crossing 'parent_shape').
        '''
        
        result = []
        
        def visit( con, shape ):
            result.append( shape )
        
        from annotations_rotund import Connection
        from wx_annotate_collection import ConnectionBetween, run_dfs
        
        cons = anns.of_type( Connection )
        cons.remove( ConnectionBetween( anns, self.shape, parent_shape ) )
        run_dfs( cons, pre_visit = visit, roots = [self.shape] )
        
        return result
    
    def origin( self, parent_shape, anns ):
        '''
        Returns the 3d world-space origin of this shape, given parent shape 'parent_shape' and
        type-list of all annotations 'anns', for transformation purposes.
        '''
        from annotations_rotund import Connection
        
        connections = [
            con for con in anns.of_type( Connection )
            if self.shape in con.regions and parent_shape in con.regions
            ]
        assert len( connections ) == 1
        return SizeManipSweepShape( self.shape, connections[0] ).get_origin()
    
    def center( self ):
        '''
        Returns a 3d point to use as the center of this shape.
        NOTE: This is only used when drawing.
        '''
        
        return self.shape.center()
    
    def mirror_point( self, pt ):
        '''
        Returns the 3d point 'pt' reflected across self.shape's symmetry plane or sheet.
        '''
        return self.shape.mirror_point( pt )
    
    def mirror_plane( self, pt ):
        '''
        Returns the 3d symmetry plane closest to the given 3d point.
        '''
        return self.shape.mirror_plane( pt )

class MirrorManipEllipsoid( MirrorManipSweepShape ):
    def mirror_point( self, pt ):
        '''
        Returns the 3d point 'pt' reflected across self.shape's symmetry plane or sheet.
        
        tested
        '''
        mirror_plane = self.mirror_plane( pt )
        return reflect_point_across_plane( pt, mirror_plane )
    
    def mirror_plane( self, pt ):
        '''
        Returns the 3d symmetry plane closest to the given 3d point.
        
        tested
        '''
        return self.shape.ellipsoid.center, self.shape.ellipsoid.axes[0]

class MirrorAcrossParent( object ):
    def __init__( self, manip = None, parent_manip = None ):
        ## 'manip' stores an object with methods .origin() and .meshcopy() and .is_garbage( anns ).
        ## .origin( parent shape, anns ) returns a 3d position for the world-space origin of the
        ##      object with given parent shape and sequence of annotations 'anns'.
        ## .meshcopy() returns a copy of itself as a mesh_gui.TriMesh object.
        ## .meshcopy_recursive( parent shape, anns ) returns a copy of itself and its entire subtree as a single mesh_gui.TriMesh object.
        ## .is_garbage( anns ) returns True if it is no longer relevant to the annotations 'anns'.
        self.manip = manip
        
        ## 'parent' stores an object with methods .mirror_point( 3d point ) and .is_garbage( anns )
        ##  and property .shape.
        ## .mirror_point( 3d point ) returns the 3d point reflected across the object's symmetry.
        ## .is_garbage( anns ) returns True if it is no longer relevant to the annotations 'anns'.
        ## .shape is a shape annotation suitable, for passing to manip.origin().
        self.parent_manip = parent_manip
        
        self.mirrored_shape = None
        
        from random import Random
        ## unique_value is used to determine the color of this annotation, and nothing else.
        self.unique_value = Random().uniform(0.,1.)
        
        self.assert_self()
    
    def __setstate__( self, unpickled ):
        self.__dict__.update( unpickled )
        
        ## To allow unpickling from older versions of the class
        if not hasattr( self, 'unique_value' ):
            from random import Random
            self.unique_value = Random().uniform(0.,1.)
        
        self.assert_self()
    
    def assert_self( self ):
        assert self.manip is None or self.manip != self.parent_manip
        if self.manip is None or self.parent_manip is None: assert self.mirrored_shape is None
    
    def garbage_collect( self, anns ):
        '''
        Removes manips that have become invalid.
        '''
        
        if self.manip.is_garbage( anns ):
            self.manip = None
            self.purge_mirrored_shape( anns )
            
        if self.parent_manip.is_garbage( anns ):
            self.parent_manip = None
            self.purge_mirrored_shape( anns )
        
        from wx_annotate_collection import ConnectionBetween
        if self.manip is not None and self.parent_manip is not None and not ConnectionBetween( anns, self.manip.shape, self.parent_manip.shape ):
            self.manip = None
            self.parent_manip = None
            self.purge_mirrored_shape( anns )
        
        self.assert_self()
    
    def is_garbage( self, anns ):
        '''
        Returns whether or not the entire object has become invalid.
        '''
        return self.manip is None or self.parent_manip is None
    
    def purge_mirrored_shape( self, anns ):
        '''
        Removes the mirrored shape from the list of annotations.
        '''
        if self.mirrored_shape is not None:
            anns.remove( self.mirrored_shape )
            self.mirrored_shape = None
    
    def apply( self, anns ):
        '''
        Adds a mirrored copy of the shape in self.manip mirrored by self.parent_manip to 'anns'.
        If a mirrored copy has been previously added, it is updated.
        
        lightly tested
        '''
        assert self.manip is not None and self.parent_manip is not None
        assert self.manip is not self.parent_manip
        if self.mirrored_shape is not None: assert self.mirrored_shape in anns
        
        if self.mirrored_shape is None:
            ## We'll give it a mesh property.
            JustMeshShape = Struct
            self.mirrored_shape = JustMeshShape()
            anns.append( self.mirrored_shape )
        
        ## Get the mirroring transformation
        origin = asarrayf( self.manip.origin( self.parent_manip.shape, anns ) )
        origin_mirror = asarrayf( self.parent_manip.mirror_point( origin ) )
        reflection_origin = .5 * ( origin + origin_mirror )
        reflection_normal = dir( origin - origin_mirror )
        
        ## Mirror self.manip as a mesh.
        #self.mirrored_shape.mesh = self.manip.meshcopy()
        ## Mirror all shapes attached to self.manip's shape, too.
        self.mirrored_shape.mesh = self.manip.meshcopy_recursive( self.parent_manip.shape, anns )
        ## Reflection in this case is sending all points across the ( reflection_origin, reflection_normal ) plane.
        self.mirrored_shape.mesh.vs = asarray( self.mirrored_shape.mesh.vs )
        vs = self.mirrored_shape.mesh.vs
        o = reflection_origin
        n = reflection_normal
        self.mirrored_shape.mesh.vs = vs - 2. * dot( vs - o, n )[:,newaxis] * n
        
        ## The mesh is now inside-out (CW); flip faces to get a CCW mesh again.
        self.mirrored_shape.mesh.faces = [ (f[0], f[2], f[1]) for f in self.mirrored_shape.mesh.faces ]
        
        ## Faces changed so topology changed.
        self.mirrored_shape.mesh.topology_changed()

class MirrorAlignAcrossParent( object ):
    def __init__( self, manip1 = None, manip2 = None, parent_manip = None ):
        ## 'manip1' and 'manip2' store objects with methods .origin() and .meshcopy() and .is_garbage( anns ).
        ## .origin( parent shape, anns ) returns a 3d position for the world-space origin of the
        ##      object with given parent shape and sequence of annotations 'anns'.
        ## .meshcopy() returns a copy of itself as a mesh_gui.TriMesh object.
        ## .meshcopy_recursive( parent shape, anns ) returns a copy of itself and its entire subtree as a single mesh_gui.TriMesh object.
        ## .subtree_shapes( parent shape, anns ) returns a list of all shapes in its subtree (including itself).
        ## .is_garbage( anns ) returns True if it is no longer relevant to the annotations 'anns'.
        self.manip1 = manip1
        self.manip2 = manip2
        
        ## 'parent' stores an object with methods .mirror_point( 3d point ) and .is_garbage( anns )
        ##  and property .shape.
        ## .mirror_plane( 3d point ) returns the 3d symmetry plane closest to the given 3d point.
        ## .is_garbage( anns ) returns True if it is no longer relevant to the annotations 'anns'.
        ## .shape is a shape annotation suitable, for passing to manip.origin().
        self.parent_manip = parent_manip
        
        from random import Random
        ## unique_value is used to determine the color of this annotation, and nothing else.
        self.unique_value = Random().uniform(0.,1.)
        
        self.assert_self()
    
    def __setstate__( self, unpickled ):
        self.__dict__.update( unpickled )
        
        ## To allow unpickling from older versions of the class
        if not hasattr( self, 'master_manips' ): self.master_manips = []
        if not hasattr( self, 'unique_value' ):
            from random import Random
            self.unique_value = Random().uniform(0.,1.)
        
        self.assert_self()
    
    def assert_self( self ):
        ## Actually, if manip1 == manip2, then we place its origin on the mirror plane.
        #assert self.manip1 != self.manip2 or self.manip1 is None
        assert self.manip1 != self.parent_manip or self.manip1 is None
        assert self.manip2 != self.parent_manip or self.manip2 is None
    
    def garbage_collect( self, anns ):
        '''
        Removes manips that have become invalid.
        '''
        
        if self.manip1.is_garbage( anns ):
            self.manip1 = None
        
        if self.manip2.is_garbage( anns ):
            self.manip2 = None
        
        if self.parent_manip.is_garbage( anns ):
            self.parent_manip = None
        
        self.assert_self()
    
    def is_garbage( self, anns ):
        '''
        Returns whether or not the entire object has become invalid.
        '''
        return self.manip1 is None or self.manip2 is None or self.parent_manip is None
    
    def apply( self, anns ):
        '''
        Adds a mirrored copy of the shape in self.manip mirrored by self.parent_manip to 'anns'.
        If a mirrored copy has been previously added, it is updated.
        
        tested
        '''
        assert self.manip1 is not None
        assert self.manip2 is not None
        assert self.parent_manip is not None
        self.assert_self()
        
        ## If manip1 == manip2, then simply project its origin onto the mirror plane.
        if self.manip1 == self.manip2:
            self._align_on_mirror_plane( anns )
        else:
            self._align_across_mirror_plane( anns )
    
    def _align_across_mirror_plane( self, anns ):
        assert self.manip1 != self.manip2
        
        ## 1 Get manip1 and manip2's closest mirroring planes to their origin points.
        ## 2 Average the two planes.
        ## 3 Solve for the smallest 3d offsets for the two origin points such that they lie on a
        ##   line parallel to the average plane's normal and bisected by the plane itself.
        ## 4 Translate both subtrees by these 3d offsets.
        
        ## 1
        origin1 = asarrayf( self.manip1.origin( self.parent_manip.shape, anns ) )
        origin2 = asarrayf( self.manip2.origin( self.parent_manip.shape, anns ) )
        
        mirror_plane1 = self.parent_manip.mirror_plane( origin1 )
        mirror_plane2 = self.parent_manip.mirror_plane( origin2 )
        
        ## 2
        mirror_plane = (
            .5 * ( asarray( mirror_plane1[0] ) + asarray( mirror_plane2[0] ) ),
            dir( asarray( mirror_plane1[1] ) + asarray( mirror_plane2[1] ) ),
            )
        #print "mirror_plane:", mirror_plane
        
        ## 3
        ## The system we want to solve is very simple.
        ## Our degrees-of-freedom are origin1 and origin2.
        ## Concatenate them together into a column vector 'x'.
        ## Let 'b' be the initial values of 'x'.
        ## Our energy is E = .5 * |x - b|^2 = .5 ( x - b )^T ( x - b ).
        ## The gradient of this is I x = b, where I is the identity matrix.
        ## Our constraints are
        ##   dot( origin1 - origin2, plane normal ) = | origin1 - origin2 |
        ## <=>
        ##   dot( origin1 - origin2, all vectors orthogonal to plane normal ) = 0
        ## <=>
        ##   dot( origin1 - origin2, ortho( plane normal ) ) = 0
        ##   dot( origin1 - origin2, cross( plane normal, ortho( plane normal ) ) ) = 0
        ## <=>
        ##   dot( origin1, ortho( plane normal ) ) - dot( origin2, ortho( plane normal ) ) = 0
        ##   dot( origin1, cross( plane normal, ortho( plane normal ) ) ) - dot( origin2, cross( plane normal, ortho( plane normal ) ) ) = 0
        ##
        ## and
        ##
        ##   dot( origin1 - plane point, plane normal ) = - dot( origin2 - plane point, plane normal )
        ## <=>
        ##   dot( origin1, plane normal ) - dot( plane point, plane normal ) = - dot( origin2, plane normal ) + dot( plane point, plane normal )
        ## <=>
        ##   dot( origin1, plane normal ) + dot( origin2, plane normal ) = 2 dot( plane point, plane normal )
        
        ## If we rotate such that plane normal equals (0,0,1), these reduce to
        ##   rot_origin1.x - rot_origin2.x = 0
        ##   rot_origin1.y - rot_origin2.y = 0
        ## and
        ##   rot_origin1.z + rot_origin2.z = 2 * rot_plane point.z
        
        b = emptyf( 6 )
        b[:3] = origin1
        b[3:] = origin2
        
        constraints = []
        ortho1 = orthogonal( mirror_plane[1] )
        ortho2 = cross( mirror_plane[1], ortho1 )
        constraints.append( ( list(range(6)), list( ortho1 ) + list( -ortho1 ), 0. ) )
        constraints.append( ( list(range(6)), list( ortho2 ) + list( -ortho2 ), 0. ) )
        constraints.append( ( list(range(6)), list( mirror_plane[1] ) * 2, 2. * dot( mirror_plane[0], mirror_plane[1] ) ) )
        #from pprint import pprint
        #print "constraints:"
        #pprint( constraints )
        
        import mymath
        x = mymath.solve_system_with_constraints(
            identity(6), b,
            linear_constraints = constraints,
            linear_constraint_type = 'hard',
            #debugging = True
            )
        
        offset1 = x[:3] - origin1
        offset2 = x[3:] - origin2
        
        ## 4
        ## Mirror self.manip1 and self.manip2 and their attached shapes as meshes.
        ## Translate each by their respective offsets.
        self.manip1.translate_recursive( self.parent_manip.shape, anns, offset1 )
        self.manip2.translate_recursive( self.parent_manip.shape, anns, offset2 )
    
    def _align_on_mirror_plane( self, anns ):
        assert self.manip1 == self.manip2
        
        origin = asarrayf( self.manip1.origin( self.parent_manip.shape, anns ) )
        mirror_plane = self.parent_manip.mirror_plane( origin )
        
        ## Project the origin onto the mirror plane.
        plane_pt = project_point_to_plane( origin, mirror_plane )
        offset = plane_pt - origin
        
        ## Mirror the manip and its attached shapes as meshes.
        ## Translate by the offset.
        self.manip1.translate_recursive( self.parent_manip.shape, anns, offset )

def symmetrize( anns, symmetry_root, frozen_shapes = None ):
    '''
    Given a sequence of annotations and a shape 'symmetry_root' to be used as the root of
    symmetrization,
    returns a list of shape annotations with symmetrization-beautification applied.
    
    Optional parameter 'frozen_shapes' contains a list of shapes to be unaffected by
    symmetrization.
    '''
    
    '''
    Make a tree rooted at the symmetric object.
    Each child has an origin, which is the attachment point to its parent.
    This origin is stored in parent-space and is scaled (or otherwise transformed)
    along with the parent when propogating symmetry.
    Nodes are scaled and tilted one step (in the graph) at a time,
    and the origin for these transformations is the attachment point to the parent,
    which remains at the origin.
    The attachment points of children are transformed along with the rest of the shape,
    and children are translated the same way as their attachment points (so they stay attached).
    This child's attachment point (from the parent's point-of-view) / origin (from the child's
    point-of-view) is the center-of-mass of the connection curve (in 3d).
    Finally, if a child's intersection curve with its parent intersects the symmetry plane,
    don't symmetrize it (skip it).  It's probably a head or tail, something that shouldn't be
    duplicated.  We could align it so its symmetry plane and the parent's are the same.
    '''
    
    ## Make a copy of all annotations so we can modify them.
    from copy import deepcopy
    anns_symmetrized = deepcopy( anns )
    
    ## A dictionary from a node to a list of its children,
    ## where each entry in its list of children is a tuple
    ## ( child node, origin in world-space on the surface of the parent ).
    node2children = gen_symmetry_node2children( anns_symmetrized )
    
    def recursive_symmetrize( node2children, symmetry_root, subtree1, subtree2 ):
        '''
        Given a dictionary mapping a shape to a list of its children,
        each child stored as the tuple ( child shape, origin in world-space ), and
        the root node as a pair ( shape that is the symmetry root node, origin in world-space ) and
        two subtree pairs ( shape that is the subtree root node, subtree root node origin in world-space ),
        modifies all shapes in both subtrees to symmetrize them.
        
        Dies if the subtrees don't have identical topology
        (each node has an identical number of children).
        '''
        
        symmetrize_pair( node2children, symmetry_root, subtree1, subtree2, frozen_shapes )
        children1 = node2children[ subtree1[0] ]
        children2 = node2children[ subtree2[0] ]
        assert len( children1 ) == len( children2 )
        matched_children = match_symmetry_children( children1, children2 )
        assert len( matched_children ) == len( children1 )
        
        for child1, child2 in matched_children:
            recursive_symmetrize( node2children, symmetry_root, child1, child2 )
    
    root_children = node2children( symmetry_root )
    
    ## Find and remove children that shouldn't be duplicated.
    shouldnt_symmetrize = find_nonsymmetry_children( symmetry_root, root_children )
    for child in shouldnt_symmetrize:
        root_children.remove( child )
    
    ## The root origin is the center of mass of its mesh's vertices.
    #assert len( symmetry_root.mesh.vs ) > 0
    #symmetry_root_origin = average( symmetry_root.mesh.vs, axis = 0, weight = symmetry_root.mesh.vertex_areas )
    
    ## Find children on either side of the symmetry plane.
    ## Match them with match_symmetry_children()
    matched_children = match_symmetry_children( [], [] )
    for child1, child2 in matched_children:
        recursive_symmetrize( node2children, symmetry_root, child1, child2 )
    
    return anns_symmetrized

def find_nonsymmetry_children( symmetry_node, children ):
    print("TODO: find_nonsymmetry_children()")
    return []

def match_symmetry_children( children1, children2 ):
    print("TODO: match_symmetry_children()")
    assert len( children1 ) == len( children2 )
    return list(zip( children1, children2 ))

def symmetrize_pair( node2children, symmetry_root, subtree1, subtree2, frozen_shapes = None ):
    print("TODO: symmetrize_pair()")

def gen_symmetry_node2children( anns ):
    '''
    Given a list of annotations, returns a dictionary from a node to a list of its children,
    where each entry in its list of children is a tuple
      ( child node, origin in world-space on the surface of the parent ).
    '''
    result = {}
    
    from annotations_sweep import SweepShape
    for shape in anns.by_type( SweepShape ):
        result[ shape ] = []
    
    return result
