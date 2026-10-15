from math import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *



from helpers import Struct, solve_for_fixed_points
from vecutil import *

kEps = 1e-8
kIterations = 10
kDebugging = False

kTestNewOverlap = False

class Annotations( list ):
    def __init__( self ):
        self.Reset()
    
    def AnnotationsFromFile( filename ):
        result = Annotations()
        result.Load( filename )
        return result
    AnnotationsFromFile = staticmethod( AnnotationsFromFile )
    
    def Reset( self, keep_texture = False ):
        del self[:]
        self.scale = 10.
        if not keep_texture:
            self.image_path = None
            self.gltex = None
    
    def sortz( self ):
        '''
        Similar to sort(), except sorts according to discrete z depth.
        '''
        ## A topological sort according to the parentage graph.
        zorder = []
        def postorder( ann, parent ):
            zorder.append( ann )
        self._run_dfs( None, postorder )
        assert len( zorder ) == len( self )
        #zorder.reverse()
        self[:] = zorder[:]
        for i, ann in enumerate( self ): ann.z = i
        
        self.sort( key = lambda entry: entry.z )
        return
    
    def MaxZ( self ):
        '''
        Returns the front-most layer's discrete z depth.
        '''
        
        if len( self ) == 0: return 0
        return max( [ ann.z for ann in self ] )
    
    def MinZ( self ):
        '''
        Returns the back-most layer's discrete z depth.
        '''
        
        if len( self ) == 0: return 0
        return min( [ ann.z for ann in self ] )
    
    def PlaceAnnotationsInDepth( self ):
        '''
        Sets the z values of AnnotationData objects' centers according to their shapes and
        discrete z depth values.
        '''
        
        if len( self ) == 0: return
        
        return self._PlaceAnnotationsInDepth_parentage()
        
        ## TODO Q: Some of my SetDepthRelativeToAnnotation() routines can handle the same
        ##         discrete z's but some may not.  Do I search in the current layer?
        ##         Should everyone just be placed according to their biggest overlaps,
        ##         and marked and then continue, and hope it is a tree?
        ## A1: mark an arbitrary object, and expand outwards to "neighboring" objects,
        ##     where neighboring means same or off-by-one z-layer and overlapping.
        ##     expand breadth-first: first place all neighbors, then recurse to each one.
        marked = set()
        remaining = set(self)
        
        self.sortz()
        seed = self[0]
        seed.ellipsoid.center[2] = 0.
        marked.add( seed )
        remaining.remove( seed )
        
        def bfs( ann ):
            assert ann in marked
            assert ann not in remaining
            
            done = False
            bfs_next_level = []
            while not done:
                done = True
                
                for child in remaining:
                    if abs( ann.z - child.z ) <= 1 and ann.Overlaps( child ):
                        child.SetDepthRelativeToAnnotation( ann, (0,0,-1) )
                        marked.add( child )
                        remaining.remove( child )
                        bfs_next_level.append( child )
                        done = False
                        break
            
            for child in bfs_next_level: bfs( child )
        
        bfs( seed )
        
        while len( remaining ) > 0:
            #for rem in remaining: break
            #assert rem == list( remaining )[0]
            rem = min( remaining, key = lambda entry: entry.z )
            
            ## If I have to place the front-most or back-most annotation, offset it by some
            ## multiple of half its "thickness" (as if shooting a ray) in z.
            padding = abs( rem.ellipsoid.axes[2,2] )
            
            ## Place it at the average of others in the layer.
            zs = [ ann.ellipsoid.center[2] for ann in marked if ann.z == rem.z ]
            ## If there are no others in the layer, look for the next-closer and farther layers.
            if len( zs ) == 0:
                ## z_m are annotations farther away than rem; they have discrete z values < rem.z,
                ## so the difference z_m.z - rem.z is negative.
                z_m = [ ann for ann in marked if ann.z < rem.z ]
                ## z_p are annotations closer than rem; they have discrete z values > rem.z,
                ## so the difference z_p.z - rem.z is positive.
                z_p = [ ann for ann in marked if ann.z > rem.z ]
                
                next_z_m_zs = []
                if len( z_m ) > 0:
                    next_z_m = max( z_m, key = lambda entry: entry.z )
                    next_z_m_zs = [ ann.ellipsoid.center[2] for ann in marked if ann.z == next_z_m.z ]
                
                next_z_p_zs = []
                if len( z_p ) > 0:
                    next_z_p = min( z_p, key = lambda entry: entry.z )
                    next_z_p_zs = [ ann.ellipsoid.center[2] for ann in marked if ann.z == next_z_p.z ]
                
                ## At least one annotation must have been placed by now!
                ## If it was placed in the same layer, we would never have taken this branch.
                ## So there must be one either above or below.
                assert len( next_z_m_zs ) > 0 or len( next_z_p_zs ) > 0
                if len( next_z_m_zs ) > 0 and len( next_z_p_zs ) > 0:
                    zs = list( next_z_m_zs ) + list( next_z_p_zs )
                elif len( next_z_m_zs ) == 0:
                    zs = [ centerz + padding*(rem.z - next_z_p.z) for centerz in next_z_p_zs ]
                elif len( next_z_p_zs ) == 0:
                    zs = [ centerz + padding*(rem.z - next_z_m.z) for centerz in next_z_m_zs ]
            
            rem.ellipsoid.center[2] = mean( zs )
            marked.add( rem )
            remaining.remove( rem )
            
            bfs( rem )
        
        return
        
        ## 1 Choose an arbitrary object in the back-most depth layer to start with.
        ## 2 Place it arbitrarily in depth (center.z = 0 is good) and mark it as "done".
        ## 3 Loop over all annotations.  If an annotations discrete depth layer is off-by-one
        ##   from a "done" annotation, place it in depth relative to the done annotation.
        ## 4 If we placed any annotations in step 3, repeat step 3.
        ## 5 Of the remaining annotations, place the back-most one arbitrarily in depth.
        ## 6 Go to step 3.
        
        ## For now, do something simpler:
        ## 1a Place all annotations in the back-most layer at depth center.z = 0.
        ## 2a For each closer layer in depth, place its annotations relative to the next-farther layer.
        
        ## 1a
        z2ann = {}
        for ann in self:
            z2ann.setdefault( ann.z, set() ).add( ann )
        z2ann = list(z2ann.items())
        z2ann.sort()
        
        for ann in z2ann[0][1]:
            ann.ellipsoid.center[2] = 0.
        print('zdepth:', z2ann[0][0])
        prev_layer = z2ann[0][1]
        
        for zdepth, layer in z2ann[1:]:
            print('zdepth:', zdepth)
            prev_layer_max_z = max( [
                ann.ellipsoid.center[2] + abs(ann.ellipsoid.axes[:,2]).sum()
                for ann in prev_layer
                ] )
            cur_layer_max_thickness = max( [ 2*abs(ann.ellipsoid.axes[:,2]).sum() for ann in prev_layer ] )
            
            for ann in layer:
                for prev_ann in prev_layer:
                    #print 'Overlaps( %s, %s )' % ( self.index(ann), self.index( prev_ann ) ),
                    if ann.Overlaps( prev_ann ):
                        ann.SetDepthRelativeToAnnotation( prev_ann, (0,0,-1) )
                        break
                else:
                    for prev_ann in layer:
                        if ann is prev_ann: continue
                        if ann.Overlaps( prev_ann ):
                            ann.SetDepthRelativeToAnnotation( prev_ann, (0,0,-1) )
                            break
                    else:
                        ann.ellipsoid.center[2] = prev_layer_max_z + .5 * cur_layer_max_thickness
            prev_layer = layer
    
    def _PlaceAnnotationsInDepth_parentage( self ):
        '''
        Sets the z values of AnnotationData objects' centers according to their parentage.
        '''
        
        self._run_dfs( lambda ann, parent: ann.SetDepthRelativeToParent(), None )
    
    def _run_dfs( self, pre_func, post_func ):
        '''
        Runs a DFS on the graph of nodes as determined by parentage.
        pre_func( ann, parent ) is called before visiting the children of ann.
        post_func( ann, parent ) is called after visiting all children of ann.
        The parameter 'parent' will be None if the node has no parent.
        '''
        
        ## Make a map from a node to nodes which are set relative to it.
        node2children = dict( [ ( ann, set() ) for ann in self ] )
        for ann in self:
            if ann.parentage is None: continue
            
            parent = ann.parentage.parent
            assert parent is not None
            
            ## Add ann to the set of its parent's children.
            node2children[ parent ].add( ann )
        
        ## Roots are nodes which have no children.
        ## TODO: Detect cycles.
        roots = [ ann for ann in self if ann.parentage is None ]
        
        def dfs( ann, parent ):
            if pre_func is not None: pre_func( ann, parent )
            
            for child in node2children[ ann ]:
                dfs( child, ann )
            
            if post_func is not None: post_func( ann, parent )
        
        for root in roots: dfs( root, None )
    
    
    def GetEdgesInfo( self ):
        edge2parentage = []
        for ann in self:
            if ann.parentage is None: continue
            
            assert ann.parentage.parent is not None
            edge = ( ann, ann.parentage.parent )
            parentage = ann.parentage
            
            ## This edge shouldn't yet be in the list (edges are undirected!)
            assert frozenset(edge) not in list(map( frozenset, list(edge2parentage.keys()) ))
            edge2parentage[ edge ] = parentage
    ## A dictionary mapping edges, which are (annotation, annotation) tuples, to the corresponding
    ## parentage structure of the first element of the edge.
    ## The order of the annotations in the edge corresponds to the point-of-view of the
    ## parentage structure, so it can be modified without checking which one is the parentage.parent.
    edges2info = property( GetEdgesInfo )
    
    def RemoveEdge( self, unordered_edge ):
        '''
        Removes the connection between two Annotation objects,
        unordered_edge[0] and unordered_edge[1].
        '''
        
        for ann in self:
            if ann is unordered_edge[0] and ann.parentage and ann.parentage.parent is unordered_edge[1]:
                ann.parentage = None
                return
            elif ann is unordered_edge[1] and ann.parentage and ann.parentage.parent is unordered_edge[0]:
                ann.parentage = None
                return
    
    def SetImage( self, image_path ):
        from GLTexture import GLTexture
        self.image_path = image_path
        try:
            self.gltex = GLTexture( image_path )
        except IOError:
            print('Failed to load texture:', image_path)
            self.image_path = None
            return
        
        ## Set the rectangle's size according to the texture's aspect ratio.
        self.tex_size = arrayf( ( self.gltex.width, self.gltex.height ) )
        self.tex_size *= 2./self.tex_size.max()
        
        ## Center the image.
        self.tex_pos = zeros(2) - .5*self.tex_size
    
    FFORMAT_V1 = 'wx_annotate_sphere file format v1'
    ## Added z to annotations
    FFORMAT_V2 = 'wx_annotate_sphere file format v2'
    ## Added parentage to annotations
    FFORMAT_V3 = 'wx_annotate_sphere file format v3'
    
    def Save( self, filename, overwrite = False ):
        import os
        if not overwrite and os.path.exists( filename ):
            print('File exists, not overwriting:', filename)
            return
        
        if overwrite and os.path.isdir( filename ):
            print('Directory with same name as desired file already exists, can\'t overwrite:', filename)
            return
        
        ## If we made it here, we're fine.
        if overwrite and os.path.exists( filename ):
            print('Overwriting file:', filename)
        
        ## This is our file format:
        fformat = ( self.FFORMAT_V3, list( self ), self.image_path )
        ## We can't save handles though because of the lambdas, so kill them.
        ## They will be regenerated anyways.
        for ann in fformat[1]: ann._handles = None
        
        import pickle
        f = open( filename, 'wb' )
        pickle.dump( fformat, f )
    
    def Load( self, filename ):
        import os
        if not os.path.isfile( filename ):
            print('Can\'t load from non-existent file:', filename)
            return
        
        self.Reset()
        import pickle
        f = open( filename, 'rb' )
        ## Specify encoding = 'latin1' to load Python 2 pickles with Python 3.
        ## https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3/28218598#28218598
        fformat = pickle.load( f, encoding = 'latin1' )
        
        if fformat[0] not in (self.FFORMAT_V1, self.FFORMAT_V2, self.FFORMAT_V3):
            print('Wrong file format, cannot load.')
            return
        self[:] = fformat[1]
        
        if fformat[0] == self.FFORMAT_V1:
            ## Fixup for scenes saved containing an older version of the annotation class.
            for ann in self:
                if not hasattr( ann, 'z' ): ann.z = 0
                if not hasattr( ann, '_handles' ): ann._handles = None
        
        if fformat[0] in (self.FFORMAT_V1, self.FFORMAT_V2):
            ## HACK to at least load scenes lacking parentage.
            for ann in self:
                assert not hasattr( ann, '_parentage' )
                ann._parentage = None
        
        self.SetImage( fformat[2] )


class AnnotationData( object ):
    def __init__( self, stroke = None, ellipse = None, ellipsoid = None, z = None ):
        self.stroke = stroke
        self.ellipse = ellipse
        self.ellipsoid = ellipsoid
        ## 2.5D depth placement
        self.z = z
        ## If not None, parentage stores a struct:
        ## ( .parent: parent annotation,
        ##   .ann_parent_front_back_signature: whether ann's front or back is attached to its parent's front or back,
        ##   .center_offset: xy offset from this annotation's center to a point on the intersection curve
        ## )
        ## 'ann_parent_front_back_signature' is a 2-element tuple of integers \in {0,1}.
        ## The first element determines whether the attachment point
        ## is ann's front or back (0 means front, 1 means back).
        ## The second element determines whether the attachment point
        ## is ann.parent's front or back (0 means front, 1 means back).
        ## By default, make ann's front attach to its parent's front.
        self._parentage = None
        self._handles = None
        
        import uuid
        self.uuid = uuid.uuid4()
    
    def __getstate__( self ):
        '''
        This method is called by the pickler, and we use it to not pickle the unpicklable handles,
        which are regenerated as needed.
        '''
        state = dict( self.__dict__ )
        state['_handles'] = None
        return state
    
    def __setstate__( self, unpickled ):
        '''
        This method is called by the unpickler, and we use it to update old-versions of this class.
        '''
        
        if 'uuid' not in unpickled:
            import uuid
            unpickled['uuid'] = uuid.uuid4()
        
        self.__dict__.update( unpickled )
    
    ## The .mesh property is needed by annotations_symbols.MirrorManip*
    def GenMesh( self ):
        '''
        Returns a TriMesh object sampling the ellipsoid.
        '''
        
        return self.ellipsoid.GenMesh()
    mesh = property( GenMesh )
    
    def GetParentage( self ): return self._parentage
    def SetParentage( self, parentage ):
        self._parentage = parentage
        self._GenHandles()
    
    parentage = property( GetParentage, SetParentage )
    
    def SetCircle( self, stroke ):
        #pt0 = asarray( stroke[0] )
        #pt1 = asarray( stroke[-1] )
        #r = .5 * mag( pt1 - pt0 )
        #self.ellipse = ( .5 * (pt0 + pt1), r * array((1,0)), r * array((0,1)) )
        self.ellipse = lsq_ellipse( stroke )
        if self.ellipse is None:
            self.ellipsoid = None
            self._handles = None
            return
        self.ellipse = list(map( asarrayf, self.ellipse ))
        self.ellipse = ( self.ellipse[0], self.ellipse[1][0], self.ellipse[1][1] )
        
        from ellipsoid import Ellipsoid
        #self.ellipsoid = Ellipsoid()
        #self.ellipsoid.center[:2] = self.circle[0]
        #self.ellipsoid.axes *= self.circle[1]
        
        self.ellipsoid = Ellipsoid()
        self.ellipsoid.center[:2] = self.ellipse[0]
        self.ellipsoid.axes[0][:2] = self.ellipse[1]
        self.ellipsoid.axes[1][:2] = self.ellipse[2]
        self.ellipsoid.axes[2] = dir( cross( self.ellipsoid.axes[0], self.ellipsoid.axes[1] ) )
        self.ellipsoid.axes[2] *= min( mag( self.ellipse[1] ), mag( self.ellipse[2] ) )
    
    def Overlaps( self, other ):
        '''
        Returns True if this annotation overlaps (in 2d)
        another AnnotationData 'other', False otherwise.
        '''
        
        assert type( other ) == type( self )
        
        result = ellipses_intersect(
            ( self.ellipse[0], self.ellipse[1:] ),
            ( other.ellipse[0], other.ellipse[1:] )
            )
        '''
        print 'Overlaps:', result
        if result:
            print self.ellipse
            print other.ellipse
        '''
        return result
    
    def SetDepthRelativeToParent( self ):
        '''
        Sets depth relative to the parent stored in self.parentage.
        '''
        if self.parentage is None:
            self.ellipsoid.center[2] = 0.
            return
        
        ## self.parentage stores the parent node, an integer indicating the relative depth offset,
        ## and an offset from self.ellipse.center to a point on the intersection curve.
        parent = self.parentage.parent
        assert parent is not None
        
        ann_ts_index = self.parentage.ann_parent_front_back_signature[0]
        parent_ts_index = self.parentage.ann_parent_front_back_signature[1]
        
        point_on_intersection = self.parentage.center_offset + self.ellipse[0]
        
        ## 1 Intersect a ray passing through point_on_intersection with both ellipsoids.
        ## 2 Set self.ellipsoid.center.z to the value that puts the 'parent_ts_index's intersection
        ## with the parent equal to the 'ann_ts_index's intersection with self.
        
        
        ## 1
        ray_pt = zeros(3)
        ray_pt[:2] = point_on_intersection[:2]
        ray_pt[2] = 1 + max( self.ellipsoid.BoundingBox()[0][2], parent.ellipsoid.BoundingBox()[0][2] )
        r = ( ray_pt, asarray( (0,0,-1) ) )
        
        if kDebugging:
            self._decs = []
            self._decs.append( ( 'line', (1,0,0), append( point_on_intersection, 0 ) + array((0,0,10)), append( point_on_intersection, 0 ) - array((0,0,10)) ) )
        
        self_ts = ray_ellipse_intersection_ts( ( self.ellipsoid.center, self.ellipsoid.axes ), r )
        self_ts.sort()
        if len( self_ts ) < 2:
            print('ERROR: connection handle is outside of the ellipse')
            ## Leave it at the last OK value.
            #self.ellipsoid.center[2] = parent.ellipsoid.center[2]
            return
        ## ts must be sorted.
        assert False not in [ t0 < t1 for t0, t1 in zip( self_ts[:-1], self_ts[1:] ) ]
        ## ts should all be positive.
        assert False not in [ t > 0 for t in self_ts ]
        
        parent_ts = ray_ellipse_intersection_ts( ( parent.ellipsoid.center, parent.ellipsoid.axes ), r )
        parent_ts.sort()
        if len( parent_ts ) < 2:
            print('ERROR: connection handle is outside of the parent ellipse')
            ## Leave it at the last OK value.
            #self.ellipsoid.center[2] = parent.ellipsoid.center[2]
            return
        ## ts must be sorted.
        assert False not in [ t0 < t1 for t0, t1 in zip( parent_ts[:-1], parent_ts[1:] ) ]
        ## ts should all be positive.
        assert False not in [ t > 0 for t in parent_ts ]
        
        if kDebugging:
            self_center = .5 * ( r[0] + self_ts[0] * r[1] + r[0] + self_ts[1] * r[1] )
            parent_center = .5 * ( r[0] + parent_ts[0] * r[1] + r[0] + parent_ts[1] * r[1] )
            
            print('---')
            print('self:', self)
            print('parent:', parent)
            print('parent_center:', parent_center)
            print('old self_center:', self_center)
            
            self._decs.append( ( 'line', (.8,.8,.4),
                self_center + (0,10,0),
                self_center - (0,10,0)
                ) )
            self._decs.append( ( 'line', (.4,.8,.8),
                parent_center + (0,10,0),
                parent_center - (0,10,0)
                ) )
        
        ## 2
        ## We want to set self.ellipsoid.center.z to the value that makes
        ## self_ts[ ann_ts_index ] equal to parent_ts[ parent_ts_index ].
        ## To get this, we want to offset self_ts[ ann_ts_index ]
        ## by ( parent_ts[ parent_ts_index ] - self_ts[ ann_ts_index ] ).
        ## The 'ts' values are along the look vector, so we need to multiply the difference
        ## by that vector.
        self.ellipsoid.center += ( parent_ts[ parent_ts_index ] - self_ts[ ann_ts_index ] ) * r[1]
    
    def SetDepthRelativeToAnnotation( self, other, look ):
        '''
        Given another AnnotationData 'other' and a 3d look vector 'look',
        sets the center of this annotation's ellipsoid so that this annotation and 'other' are
        nicely touching.
        '''
        
        assert type( other ) == type( self )
        
        look = asarray( look ).ravel()
        assert len( look ) == 3
        assert mag2( look ) > kEps
        look = dir( look )
        ## Assume the look vector is the negative z-axis.
        assert mag2( look - asarray( (0,0,-1) ) ) < kEps
        
        ## TODO Q: Should I scrap the 'in_front' parameter and place them relative to each other
        ##         according to their discrete z's, perhaps only if they Overlaps()?
        ## A: Yes.  Place according to discrete z's.  If they don't overlap,
        ##    place them so that they don't intersect.
        #raise NotImplementedError, 'TODO'
        #self.ellipsoid.center[2] = other.ellipsoid.center[2] + .5 * abs(other.ellipsoid.axes[:,2]).sum()
        
        #self.__SetDepthRelativeToAnnotation_fancy( other, look )
        #self.__SetDepthRelativeToAnnotation_simple( other, look )
        self.__SetDepthRelativeToAnnotation_medium( other, look )
    
    def __SetDepthRelativeToAnnotation_simple( self, other, look ):
        '''
        Overlap the two ellipses in z via a simple calculation.
        '''
        
        ## 1 Find the two axes of 'self' and 'other' with the largest absolute value z component.
        ## 2 Place self.ellipsoid.center.z such that
        ##  self.center.z + dz * magic * self.largest_axis_z_offset = other.center.z + dz * magic * other.largest_axis_z_offset
        ## where self.largest_axis_z_offset is the z-offset from center.z to the tip of the axis
        ## with the largest z component in the +look direction,
        ## other.largest_axis_z_offset is the same in the -look direction,
        ## the sign dz = self.discrete_z - other.discrete_z (ok an integer not really a sign),
        ## and magic = .85.
        ## TODO: clamp the largest_axis_z_offset values so that the tip of the ellipsoid never crosses center.z.
        ## <=>
        ##  self.center.z = other.center.z + dz * magic * ( so * other.largest_axis_z_offset - self.largest_axis_z_offset )
        
        self_max_abs_z_index = argmax( abs( self.ellipsoid.axes[:,2] ) )
        other_max_abs_z_index = argmax( abs( other.ellipsoid.axes[:,2] ) )
        
        self_largest_axis_z_offset = abs( self.ellipsoid.axes[ self_max_abs_z_index, 2 ] )
        other_largest_axis_z_offset = -abs( self.ellipsoid.axes[ other_max_abs_z_index, 2 ] )
        
        dz = self.z - other.z
        magic = .5
        ## TODO: This isn't very satisfying.
        
        print('old center.z:', self.ellipsoid.center[2])
        self.ellipsoid.center[2] = other.ellipsoid.center[2] + dz * magic * ( self_largest_axis_z_offset - other_largest_axis_z_offset )
        print('new center.z:', self.ellipsoid.center[2])
    
    def __SetDepthRelativeToAnnotation_medium( self, other, look ):
        '''
        Overlap the ellipses in z via a binary search that makes a ray shot through both be
        a fraction of the thickness of 'self'.
        '''
        
        ## 1 Find com, the center-of-mass of the 2d overlap of self.ellipse and other.ellipse.
        ## 2 Construct a ray r through com along the look vector.
        ## 3 Find self_thickness, the thickness of self.ellipsoid along the ray r,
        ##   and other other_thickness, the thickness of other.ellipsoid along the ray r.
        ## 4 Set self.ellipsoid.center.z to such that they overlap along the ray r by the amount
        ##   .5 * min( self_thickness, other_thickness )
        
        
        ## 1
        com = center_of_mass_of_ellipses_intersection(
            ( self.ellipse[0], self.ellipse[1:] ),
            ( other.ellipse[0], other.ellipse[1:] )
            )
        com = asarray( com )
        
        ## 2
        com3 = zeros(3)
        com3[:2] = com[:2]
        com3[2] = 1 + max( self.ellipsoid.BoundingBox()[0][2], other.ellipsoid.BoundingBox()[0][2] )
        r = ( com3, look )
        
        if kDebugging:
            self._decs = []
            self._decs.append( ( 'line', (1,0,0), append( com, 0 ) + array((0,0,10)), append( com, 0 ) - array((0,0,10)) ) )
        
        ## 3
        self_ts = ray_ellipse_intersection_ts( ( self.ellipsoid.center, self.ellipsoid.axes ), r )
        self_ts.sort()
        assert len( self_ts ) == 2
        ## ts must be sorted.
        assert False not in [ t0 < t1 for t0, t1 in zip( self_ts[:-1], self_ts[1:] ) ]
        ## ts should all be positive.
        assert False not in [ t > 0 for t in self_ts ]
        self_thickness = self_ts[1] - self_ts[0]
        
        other_ts = ray_ellipse_intersection_ts( ( other.ellipsoid.center, other.ellipsoid.axes ), r )
        other_ts.sort()
        assert len( other_ts ) == 2
        ## ts must be sorted.
        assert False not in [ t0 < t1 for t0, t1 in zip( other_ts[:-1], other_ts[1:] ) ]
        ## ts should all be positive.
        assert False not in [ t > 0 for t in other_ts ]
        other_thickness = other_ts[1] - other_ts[0]
        
        if kDebugging:
            self_center = .5 * ( r[0] + self_ts[0] * r[1] + r[0] + self_ts[1] * r[1] )
            other_center = .5 * ( r[0] + other_ts[0] * r[1] + r[0] + other_ts[1] * r[1] )
            
            print('---')
            print('self:', self)
            print('other:', other)
            print('other_center:', other_center)
            print('old self_center:', self_center)
            
            self._decs.append( ( 'line', (.8,.8,.4),
                self_center + (0,10,0),
                self_center - (0,10,0)
                ) )
            self._decs.append( ( 'line', (.4,.8,.8),
                other_center + (0,10,0),
                other_center - (0,10,0)
                ) )
        
        ## 4
        if not kTestNewOverlap:
            magic = .5
            overlap = magic * min( self_thickness, other_thickness )
        else:
            ## Steps 3/4 alt
            
            ''''
            ## Find which (normalized) axis of each ellipsoid is most orthogonal to the look vector,
            ## and overlap by a fraction of that axis's thickness.
            self_max_orthogonal_axis_index = argmax( [ abs( dot( look, dir( axis ) ) ) for axis in self.ellipsoid.axes ] )
            other_max_orthogonal_axis_index = argmax( [ abs( dot( look, dir( axis ) ) ) for axis in other.ellipsoid.axes ] )
            overlap = min(
                abs( self.ellipsoid.axes[ self_max_orthogonal_axis_index ][2] ),
                abs( other.ellipsoid.axes[ other_max_orthogonal_axis_index ][2] )
                )
            
            axis_mags = list( map( mag, self.ellipsoid.axes ) ) + list( map( mag, other.ellipsoid.axes ) )
            assert overlap >= min( axis_mags )
            overlap = min( axis_mags )
            overlap = min( [ overlap, .5 * self_thickness, .5 * other_thickness ] )
            '''
            
            self_mags_axes = [ ( mag( axis ), tuple(axis) ) for axis in self.ellipsoid.axes ]
            other_mags_axes = [ ( mag( axis ), tuple(axis) ) for axis in other.ellipsoid.axes ]
            ## Multiply the axis by two to get the long length of the ellipsoid.
            shorter_long_axis = min( max( self_mags_axes ), max( other_mags_axes ) )
            
            ## Overlap in z by 15 percent of the shorter long axis.
            magic = .15
            overlap = magic * 2. * abs( dot( asarrayf( shorter_long_axis[1] ), look ) )
            
            ## But don't overlap more than half of the shortest short axis.
            ## NOTE: This doesn't prevent the ellipsoids from sticking through each other.
            ##       To do that, I would also need to clip the overlap to min( self/other_thickness ).
            shorter_short_axis = min( min( self_mags_axes ), min( other_mags_axes ) )
            if shorter_short_axis[0] < magic * shorter_long_axis[0]:
                overlap = abs( dot( asarrayf( shorter_short_axis[1] ), look ) )
        
        ## With this setting, the two ellipsoid end up just barely touching.
        #overlap = 0.
        
        ## We're going to fold the cases where they aren't on neighboring layers into the off-by-one
        ## case by changing overlap.
        if abs( self.z - other.z ) > 1:
            padding = .1 * min( self_thickness, other_thickness )
            overlap = -( abs( self.z - other.z ) - 1 ) * padding
        
        ## If self.z == other.z + 1,
        ## we want the ray to leave self.ellipsoid 'overlap' after it enters other.ellipsoid.
        ## In other words, we want self_ts[1] = other_ts[0] + overlap.
        ## To get this, we want to offset self_ts[1] by ( other_ts[0] + overlap - self_ts[1] ).
        ## To offset self_ts[1], we offset the self.ellipsoid.center.z.
        if self.z > other.z:
            self.ellipsoid.center += ( other_ts[0] + overlap - self_ts[1] ) * look
        
        ## If self.z == other.z - 1,
        ## we want the ray to leave other.ellisoid 'overlap' after it enters self.ellipsoid.
        ## <=>
        ##  other_ts[1] = self_ts[0] + overlap
        ## <=>
        ##  self_ts[0] = other_ts[1] - overlap
        ## <=> offset self.ellipsoid.center.z by
        ##  ( other_ts[1] - overlap - self_ts[0] ) * look
        elif self.z < other.z:
            self.ellipsoid.center += ( other_ts[1] - overlap - self_ts[0] ) * look
        
        ## If self.z == other.z,
        ## we want to align the midpoints of the ray's interval in self.ellipsoid and other.ellipsoid.
        ## <=>
        ##  .5 * ( self_ts[0] + self_ts[1] ) == .5 * ( other_ts[0] + other_ts[1] )
        ## where the left-hand-side can be written as an offset from self.ellipsoid.center.z
        ## <=> offset self.ellipsoid.center.z by
        ##  .5 * ( ( other_ts[0] + other_ts[1] ) - ( self_ts[0] + self_ts[1] ) )
        elif self.z == other.z:
            self.ellipsoid.center += .5 * ( ( other_ts[0] + other_ts[1] ) - ( self_ts[0] + self_ts[1] ) ) * look
        
        if kDebugging:
            self_ts = ray_ellipse_intersection_ts( ( self.ellipsoid.center, self.ellipsoid.axes ), r )
            self_ts.sort()
            self_center = .5 * ( r[0] + self_ts[0] * r[1] + r[0] + self_ts[1] * r[1] )
            print('new self_center:', self_center)
    
    def __SetDepthRelativeToAnnotation_fancy( self, other, look ):
        ## Sample uniformly the intersection of their bounding boxes.  If the point is inside both
        ## ellipsoids, contribute its normal (project to the surface).  Take the average normal for
        ## each.
        ## Do this for all the relative placements of 'self.ellipsoid.center.z'.  Do a binary search
        ## (or two binary searches, since I think there are in general two placements --- or maybe
        ## cull the behind one if we want front).
        
        other_bbox = other.ellipsoid.BoundingBox()
        
        def evaluate_z( z ):
            old_center_z = self.ellipsoid.center[2]
            self.ellipsoid.center[2] = z
            self_bbox = self.ellipsoid.BoundingBox()
            
            both_bbox = (
                asarray( [ self_bbox[0], other_bbox[0] ] ).max(0),
                asarray( [ self_bbox[1], other_bbox[1] ] ).min(0)
                )
            
            N = 500
            from random import uniform
            self_normal = zeros(3)
            other_normal = zeros(3)
            hits = 0
            for i in range(500):
                pt = (
                    uniform( both_bbox[1][0], both_bbox[0][0] ),
                    uniform( both_bbox[1][1], both_bbox[0][1] ),
                    uniform( both_bbox[1][2], both_bbox[0][2] )
                    )
                if(
                    ellipse_contains_point( ( self.ellipsoid.center, self.ellipsoid.axes ), pt )
                    and
                    ellipse_contains_point( ( other.ellipsoid.center, other.ellipsoid.axes ), pt )
                    ):
                    sn = ellipse_normal_at_point( ( self.ellipsoid.center, self.ellipsoid.axes ), pt )
                    on = ellipse_normal_at_point( ( other.ellipsoid.center, other.ellipsoid.axes ), pt )
                    hits += 1
                    self_normal += sn
                    other_normal += on
            print('hits:', hits)
            if hits == 0:
                pt = .5 * ( self.ellipsoid.center + other.ellipsoid.center )
                hits = 1
                self_normal = ellipse_normal_at_point( ( self.ellipsoid.center, self.ellipsoid.axes ), pt )
                other_normal = ellipse_normal_at_point( ( other.ellipsoid.center, other.ellipsoid.axes ), pt )
            
            self_normal = dir( self_normal )
            other_normal = dir( other_normal )
            ## This isn't an energy, it is a function which is 0 when the normals point
            ## opposite each other, which is what we want.  Use find_zero.
            #result = dot( self_normal, other_normal ) - -1
            ## This is an energy.  Use find_minimum.
            result = ( dot( self_normal, other_normal ) - -1 )**2
            
            self.ellipsoid.center[2] = old_center_z
            return result
        
        #from optutil import find_zero as optim
        from optutil import find_minimum as optim
        z, converged = optim( evaluate_z, other.ellipsoid.center[2] + .01*(self.z - other.z), .01, verbose = True )
        self.ellipsoid.center[2] = z
    
    def SetAxis0FromRay3D( self, ray ):
        '''
        Axis 0 defaults to the X axis.
        Given a ray (pt and direction tuple), there is only one configuration that preserves
        the on-screen orientation of Axis 1, which defaults to the y axis,
        and the edge lengths of axes 0 and 2.
        In fact, there is no solution if the ray doesn't intersect the sphere with radius axis 0.
        '''
        
        ## 1 The new axis 0 is the intersection of the ray with the sphere with radius axis 0.
        ## 2 The new axis 2 is the cross product of axes 0 and 1.
        ## 3 The new axis 1 is the cross product of axes 2 and 0.
        
        ## 1
        mags_before = [ mag( self.ellipsoid.axes[i] ) for i in range(3) ]
        #print 'before axis0:', self.ellipsoid.axes[0]
        self.ellipsoid.axes[0] = mags_before[0] * dir( closest_point_to_sphere_on_ray(
            ( self.ellipsoid.center, mag( self.ellipsoid.axes[0] ) ),
            ray
            )[0] - self.ellipsoid.center )
        #print 'after axis0:', self.ellipsoid.axes[0]
        ## 2
        #self.ellipsoid.axes[1] = mags_before[1] * dir( cross( self.ellipsoid.axes[2], self.ellipsoid.axes[0] ) )
        ## Note that because we are asking for a cos(angle) of 0, closest_normal_with_cosangle_to_vector_on_plane()
        ## will always succeed.
        self.ellipsoid.axes[1] = mags_before[1] * closest_normal_with_cosangle_to_vector_on_plane(
            (1,0,0), ## I want a vector in the yz plane
            self.ellipsoid.axes[0], 0, ## perpendicular to axes[0]
            self.ellipsoid.axes[1] ## and closest to our current axis in the yz plane
            )
        ## 3
        self.ellipsoid.axes[2] = mags_before[2] * dir( cross( self.ellipsoid.axes[0], self.ellipsoid.axes[1] ) )
        # assert False not in [ fabs( mags_before[i]**2 - mag2( self.ellipsoid.axes[i] ) ) < 1e-5 for i in range(3) ]
        #print [ ( mags_before[i]**2, mag2( self.ellipsoid.axes[i] ) ) for i in range(3) ]
    
    def _StretchAxisInDepth( self, which_axis, new_length, look ):
        '''
        Scales the axis specified by 'which_axis' so that the axis has
        length 'new_length' when viewed along the look vector 'look'.
        This does not consider the scale of the scene, which this object doesn't know about.
        '''
        
        assert which_axis >= 0 and which_axis < 3
        assert new_length > kEps
        assert( len( look ) ) == 3
        ## normalize look
        look = dir( asarrayf( look ) )
        
        ## 1 Separate the given axis into components parallel and perpendicular to the look
        ##   direction.
        ## 2 Rescale the parallel component so that the combination has the given length.
        
        ## 1
        axis = self.ellipsoid.axes[ which_axis ]
        par = dot( axis, look ) * look
        perp = axis - par
        assert dot( par, perp ) < kEps
        
        ## 2
        new_par_len = new_length**2 - mag2( perp )
        if new_par_len < kEps:
            print('Desired length is shorter than just the image plane length!  Using that length.')
            new_par_len = 0.
        else:
            new_par_len = sqrt( new_par_len )
        
        ## In case there is no perpendicular component, use the negative look vector (coming out of
        ## the screen).
        if mag2( par ) < kEps: par = -look
        else: par = dir( par )
        self.ellipsoid.axes[ which_axis ] = new_par_len * par + perp
    
    def SetAxisLengthInDepth( self, which_axis, new_length, look ):
        self._StretchAxisInDepth( which_axis, new_length, look )
        ## Find new other axes so that their projections don't change.
        
        print('Setting Axis %s length' % (which_axis,))
        look = asarrayf( look )
        mags_before = [ mag( self.ellipsoid.axes[i] ) for i in range(3) ]
        
        if which_axis == 1:
            
            ## Set axes[0] to something which projects to the same point but is perpendicular to
            ## axes[which_axis].
            ## NOTE: I am favoring axis 0 over axis 2.
            self.ellipsoid.axes[0] = closest_point_along_ray_with_cosangle_to_vector(
                (self.ellipsoid.axes[0] - 2*max( mags_before ) * look, look), ## A ray passing through the existing point along the look vector
                0, self.ellipsoid.axes[which_axis], ## perpendicular to axes[which_axis]
                True ## The closest point along the ray in either direction
                )
            ## Axes[2] is now completely determined. (See SetAxis0FromRay3D().)
            self.ellipsoid.axes[2] = mags_before[2] * dir( cross( self.ellipsoid.axes[0], self.ellipsoid.axes[1] ) )
        else:
            
            ## Adapted from SetAxis0FromRay3D().
            
            self.ellipsoid.axes[1] = mags_before[1] * closest_normal_with_cosangle_to_vector_on_plane(
                (1,0,0), ## I want a vector in the yz plane
                self.ellipsoid.axes[which_axis], 0, ## perpendicular to axes[0]
                self.ellipsoid.axes[1] ## and closest to our current axis in the yz plane
                )
            last_axis = [0,2]
            last_axis.remove( which_axis )
            last_axis = last_axis[0]
            self.ellipsoid.axes[last_axis] = (-1)**(1+last_axis//2) * mags_before[last_axis] * dir( cross( self.ellipsoid.axes[which_axis], self.ellipsoid.axes[1] ) )
    
    def GetTiltRadians( self ):
        '''
        Returns the tilt angle of the circular cross section in radians.  An angle of 0 means that the circular cross section spanned
        by axes 1 and 2 is perpendicular to the image plane.
        
        used
        '''
        ## Following the comments in SweepShape.__init__() and
        ## "discussion/SweepShape-tilt_radians.pdf":
        normal = dir( self.ellipsoid.axes[0] )
        #tilt_radians = asin( max( -1., min( 1., dir(normal)[2] ) ) )
        tilt_radians = -atan( normal[2] / mag(normal[:2]) )
        return tilt_radians
    
    def SetTiltRadians( self, radians ):
        '''
        Sets the tilt angle in radians of the circular cross section.
        
        used
        '''
        x = sin( radians ) * mag( self.ellipsoid.axes[2] )
        ## Ellipsoids' orientations flip so that axes[2] always points in positive z, so we have to multiply by the -sign of axes[0].z.
        neg = 1. if ( self.ellipsoid.axes[0][2] < 0. ) else -1.
        self.SetCrossSectionFromStroke( 1, True, [ self.ellipsoid.center[:2] + x*neg*dir( self.ellipsoid.axes[2][:2] ) ], (0,0,-1) )
        #print 'after SetTiltRadians( %s ), GetTiltRadians() is' % (radians,), self.GetTiltRadians()
    
    def GetRadius( self ):
        '''
        Returns the radius of the circular cross section.
        
        used
        '''
        return mag( self.ellipsoid.axes[1] )
    
    def SetRadius( self, radius ):
        '''
        Sets the radius of the circular cross section.
        
        used
        '''
        scale = radius / mag( self.ellipsoid.axes[1] )
        self.ellipsoid.axes[1:,:] *= scale
    
    def length3d( self ): return self.GetLength3D()
    def GetLength3D( self ):
        '''
        Returns the length of the ellipsoid's 3D long axis.
        
        used
        '''
        return 2*mag( self.ellipsoid.axes[0] )
    
    def set_length3d( self, length, origin = (0,0,0) ): return self.SetLength3D( length, origin )
    def SetLength3D( self, length, origin = (0,0,0) ):
        '''
        Sets the length of the ellipsoid's 3D long axis.
        
        used
        '''
        assert length > 0.
        
        assert len( origin ) == 3
        origin = asarrayf( origin ).ravel()
        
        scale = length / ( 2*mag( self.ellipsoid.axes[0] ) )
        self.ellipsoid.axes[0,:] *= scale
        self.ellipse = ( self.ellipse[0], scale * asarrayf( self.ellipse[1] ), self.ellipse[2] )
        
        self.ellipsoid.center = scale * ( self.ellipsoid.center - origin ) + origin
        self.ellipse = ( scale * ( asarrayf( self.ellipse[0] ) - origin[:2] ) + origin[:2], self.ellipse[1], self.ellipse[2] )
    
    def GetClosestAxisMidpointInDepth( self, which_axis, depth_along ):
        return self.ellipsoid.center + .5 * self.ellipsoid.axes[which_axis] * [ 1, -1 ][ dot( self.ellipsoid.axes[which_axis], depth_along ) > 0 ]
    
    def SetCrossSectionFromStroke( self, which_axis, positive_tip, stroke, look ):
        '''
        Updates the specified axis
        given an axis index 'which_axis', a boolean indicating whether we start at the axis's
        positive or negative tip, and a list of 2d world-space points 'stroke'
        drawn from the direction of the 3d vector 'look'.
        '''
        
        #debugger()
        
        assert which_axis in (0,1,2)
        assert positive_tip == bool(positive_tip)
        
        stroke = asarray( stroke )
        ## We only need one point in the stroke, since we only care about its maximum height
        ## relative to axis 'which_axis'.
        assert stroke.shape[0] > 0
        ## 2 dimensional points.
        assert stroke.shape[1] == 2
        
        look = asarray( look ).ravel()
        assert look.shape[0] == 3
        assert len( look.shape ) == 1
        
        assert mag2( dir( look ) - asarray( (0,0,-1) ) ) < kEps
        
        ## We want 'which_axis' to keep the same xy length, but lie in the
        ## z = 0 plane.
        ## We want axis 2 (which defaults to the z axis) to have the
        ## same xyz length as 'which_axis' but be tilted so that its tip
        ## passes through the extremum point of the stroke (or near to it
        ## to keep axis 2 perpendicular to 'which_axis').
        ## We set the positive or negative axis tip of axis 2 depending
        ## on whether the extremum is above or below the baseline.
        ## We want the axis perpendicular to 'which_axis' in the xy-plane
        ## to be perpendicular to the other two axes, and have length such that
        ## its xy-length remains unchanged across this function call.
        
        ## The recipe:
        ## 1 Remove z from axis 'which_axis'
        ## 2 Set axis 2 to the z axis, but rotate such that its (positive or
        ##   negative) tip passes through the extremum point (or close to it,
        ##   to keep axis 2 perpendicular to axis 'which_axis'.
        ## 3 Set the remaining axis to be perpendicular to axis 2 and axis
        ##   'which_axis', with length such that its xy length remains unchanged.
        
        
        if which_axis == 2:
            print("I can't set the z axis")
            return
        
        old_axes = self.ellipsoid.axes.copy()
        
        ## 1
        self.ellipsoid.axes[which_axis][2] = 0.
        mag_which_axis = mag( self.ellipsoid.axes[which_axis] )
        assert mag_which_axis > kEps
        
        ## 2
        ## The line between the first and last point of endpoints
        ## determine the baseline 0 for magic_zvals().
        #endpoints = stroke
        endpoints = [
            ( self.ellipsoid.center + self.ellipsoid.axes[which_axis] )[:2],
            ( self.ellipsoid.center - self.ellipsoid.axes[which_axis] )[:2]
            ]
        from annotations_rotund import magic_zval
        zvals = [ magic_zval( pt, endpoints ) for pt in stroke ]
        ## Get the index of the absolute value biggest zvalue.
        extremumi = argmax( list(map( abs, zvals )) )
        ## Find the closest point to the plane perpendicular to axis 'which_axis'.
        perp_extremum = project_point_to_plane(
            list( stroke[ extremumi ] ) + [0.],
            ( self.ellipsoid.center, self.ellipsoid.axes[which_axis] )
            )
        ## axis 'which_axis' and the unprojected point have z = 0, so the result should as well.
        assert perp_extremum[2] < kEps
        
        ## Now set the tip of the default z-axis of the ellipsoid to
        ## pass through the extremum point in xy.
        perp_extremum[2] = 2*mag_which_axis + self.ellipsoid.center[2]
        closest_sphere_pt = closest_point_on_sphere_to_ray(
            ( self.ellipsoid.center, mag_which_axis ),
            ( perp_extremum, (0,0,-1) )
            )
        self.ellipsoid.axes[2] = closest_sphere_pt[0] - self.ellipsoid.center
        assert dot( self.ellipsoid.axes[which_axis], self.ellipsoid.axes[2] ) < kEps
        kMinRiseCosAngle = cos( pi/10 )
        if mag( self.ellipsoid.axes[2][:2] ) / mag( self.ellipsoid.axes[2] ) > kMinRiseCosAngle:
            print("Refusing to create an infinitely long (in depth) shape.  Clipping to the maximum allowed tilt in depth.")
            
            ## This computes the same thing as the next block.
            #self.ellipsoid.axes[2][2] = 0.
            #self.ellipsoid.axes[2] = mag_which_axis * dir( self.ellipsoid.axes[2] )
            #self.ellipsoid.axes[2] = rotate( self.ellipsoid.axes[2], self.ellipsoid.axes[which_axis], acos( kMinRiseCosAngle ) )
            #self.ellipsoid.axes[2][2] = abs( self.ellipsoid.axes[2][2] )
            
            ## Solve for z so that cosine (adjacent / hypoteneuse) equals kMinRiseCosAngle.
            xymag2 = mag2( self.ellipsoid.axes[2][:2] )
            self.ellipsoid.axes[2][2] = sqrt( ( xymag2 - kMinRiseCosAngle**2 * xymag2 ) / kMinRiseCosAngle**2 )
            self.ellipsoid.axes[2] = mag_which_axis * dir( self.ellipsoid.axes[2] )
        
        ## Determine the tilt of the ellipse.
        ## The picked point on axis[2] should always be closer to the camera (positive z)
        assert self.ellipsoid.axes[2][2] >= 0.
        
        ## 3
        other_axis = 0 if which_axis == 1 else 1
        
        def xymag_silhouette( axes, which_axis ):
            axis0 = array( ( mag( axes[ which_axis ][:2] ), axes[ which_axis ][2] ) )
            axis1 = array( ( mag( axes[ 2 ][:2] ), axes[ 2 ][2] ) )
            tan_theta = axis1[0] / axis0[0]
            theta = atan( tan_theta )
            #print 'should be 0:', -sin(theta) * axis0 + cos(theta) * axis1
            return cos(theta) * axis0[0] + sin(theta) * axis1[0]
        
        xymag_other_axis_silhouette = xymag_silhouette( old_axes, other_axis )
        xymag_other_axis = mag( self.ellipsoid.axes[other_axis][:2] )
        self.ellipsoid.axes[other_axis] = cross(
            self.ellipsoid.axes[which_axis],
            self.ellipsoid.axes[2]
            )
        
        assert dot( self.ellipsoid.axes[other_axis], self.ellipsoid.axes[2] ) < kEps
        assert dot( self.ellipsoid.axes[other_axis], self.ellipsoid.axes[which_axis] ) < kEps
        
        expect( abs( xymag_other_axis - mag( self.ellipsoid.axes[other_axis][:2] ) ) > kEps )
        
        '''
        ## Keep axis 'other_axis's xy magnitude unchanged.
        other_axis_scale = xymag_other_axis / mag( self.ellipsoid.axes[other_axis][:2] )
        assert other_axis_scale > kEps
        
        self.ellipsoid.axes[other_axis] *= other_axis_scale
        
        assert abs( xymag_other_axis - mag( self.ellipsoid.axes[other_axis][:2] ) ) < kEps
        '''
        
        ## Keep axis 'other_axis's silhouette xy magnitude unchanged.
        def xymag_for_silhouette_xymag( axes, which_axis, silhouette_xymag ):
            ## Uses Newton's Method to solve for 'xymag' such that the silhouette function
            ## evaluates to 'silhouette_xymag'.
            
            STEP_TOL = 1e-5
            MAX_ITER = 10
            
            l = silhouette_xymag
            b = mag( axes[2][:2] )
            x = mag( axes[ which_axis ][:2] )
            
            dstep = 2*STEP_TOL
            iter = 0
            while dstep > STEP_TOL and iter < MAX_ITER:
                iter += 1
                
                theta = atan( b/x )
                cos_theta = cos( theta )
                sin_theta = sin( theta )
                
                fx = cos_theta * x + sin_theta * b - l
                
                chain_term = ( 1 / ( 1 + b**2 / x**2 ) ) * ( -b / x**2 )
                
                dfx = x * ( -sin_theta ) * chain_term + cos_theta + b * cos_theta * chain_term
                
                dstep = abs( fx/dfx )
                #print 'fx:', fx
                #print 'dstep:', dstep
                x = x - fx/dfx
            
            if dstep > STEP_TOL:
                print('diverged!')
            else:
                #print 'converged in %d iterations' % (iter,)
                pass
            
            return x
        
        self.ellipsoid.axes[other_axis] *= (
            xymag_for_silhouette_xymag( self.ellipsoid.axes, other_axis, xymag_other_axis_silhouette )
            /
            mag( self.ellipsoid.axes[other_axis][:2] )
            )
        
        assert abs( xymag_other_axis_silhouette - xymag_silhouette( self.ellipsoid.axes, other_axis ) ) < kEps
        #print 'sil diff:', abs( xymag_other_axis_silhouette - xymag_silhouette( self.ellipsoid.axes, other_axis ) )
        
        
        assert dot( self.ellipsoid.axes[other_axis], self.ellipsoid.axes[2] ) < kEps
        assert dot( self.ellipsoid.axes[other_axis], self.ellipsoid.axes[which_axis] ) < kEps
        
        
        ## Update the 2d ellipse:
        self.ellipse = list( self.ellipse )
        self.ellipse[0] = self.ellipsoid.center[:2]
        self.ellipse[1] = self.ellipsoid.axes[which_axis][:2]
        self.ellipse[2] = xymag_other_axis_silhouette * dir( self.ellipsoid.axes[other_axis][:2] )
        self.ellipse = tuple( self.ellipse )
    
    
    Axis0Handle = 0
    Axis1Handle = 1
    CrossSectionHandle = 2
    DiscreteZHandle = 3
    WholeEllipseHandle = 4
    ParentIntersectionHandle = 5
    
    def _GenHandles( self ):
        '''
        Returns a list of Handle objects, each with the following properties:
            annotation: the AnnotationData object controlled by the Handle
            kind: one of { Axis0Handle, Axis1Handle, CrossSectionHandle }
            get_pos: a method returning the world-space position of the handle.
            translate: a method which takes a world-space offset vector for moving the handle.
        '''
        
        assert self.ellipse is not None and self.ellipsoid is not None
        
        self._handles = []
        
        from pprint import pprint
        
        ## Translate the whole ellipse
        def drag_whole_ellipse( offset ):
            if kDebugging: print('drag whole annotation:', self)
            offset = asarray( offset )
            self.ellipsoid.center[:2] += offset
            
            ## I don't think I need to call SetCrossSectionFromStroke(), I think all my
            ## structures are translation invariant.
            #self.ellipse = list( self.ellipse )
            #self.ellipse[0] = asarray( self.ellipse[0] ) + offset
            #self.ellipse = tuple( self.ellipse )
            ## UPDATE: I do need to call SetCrossSectionFromStroke(); the ellipsoid appears a little
            ##         off when I just translate.
            self.SetCrossSectionFromStroke( 1, True, [ self.ellipsoid.center[:2] + self.ellipsoid.axes[2][:2] ], (0,0,-1) )
        
        handle = Struct()
        handle.annotation = self
        handle.kind = self.WholeEllipseHandle
        handle.get_pos = lambda: self.ellipsoid.center
        handle.translate = drag_whole_ellipse
        self._handles.append( handle )
        
        
        ## Axis 0 handles
        def drag_axis0( offset, positive_axis ):
            if mag2( offset ) < kEps: return
            
            offset = asarray( list( offset ) + [0.] )
            xformed = solve_for_fixed_points(
                [
                    self.ellipsoid.center + self.ellipsoid.axes[0],
                    self.ellipsoid.center + self.ellipsoid.axes[1],
                    self.ellipsoid.center + self.ellipsoid.axes[2],
                    self.ellipse[0] + ( self.parentage.center_offset if self.parentage is not None else 0. ),
                    self.ellipsoid.center
                ],
                (self.ellipsoid.center, self.ellipsoid.center + ( 1 if positive_axis else -1 ) * self.ellipsoid.axes[0]),
                (self.ellipsoid.center, self.ellipsoid.center + ( 1 if positive_axis else -1 ) * self.ellipsoid.axes[0] + offset)
                )
            self.ellipsoid.center[:2] = xformed[-1]
            for i in range(3): self.ellipsoid.axes[i][:2] = asarray( xformed[i] ) - self.ellipsoid.center[:2]
            if self.parentage is not None: self.parentage.center_offset[:] = xformed[-2] - self.ellipse[0]
            
            ## Rescale axes.z appropriately.
            ## TODO: This moves the whole ellipsoid, I'm not sure why.
            ## UPDATE: This works perfectly except when axes[2] has no xy.
            self.SetCrossSectionFromStroke( 1, True, [ self.ellipsoid.center[:2] + self.ellipsoid.axes[2][:2] ], (0,0,-1) )
        
        handle = Struct()
        handle.annotation = self
        handle.kind = self.Axis0Handle
        handle.get_pos = lambda: self.ellipsoid.center + self.ellipsoid.axes[0]
        #handle.translate = lambda offset: pprint( offset )
        handle.translate = lambda offset: drag_axis0( offset, True )
        self._handles.append( handle )
        
        handle = Struct()
        handle.annotation = self
        handle.kind = self.Axis0Handle
        handle.get_pos = lambda: self.ellipsoid.center - self.ellipsoid.axes[0]
        #handle.translate = lambda offset: pprint( offset )
        handle.translate = lambda offset: drag_axis0( offset, False )
        self._handles.append( handle )
        
        
        ## Axis 1 handles
        def drag_axis1( offset, positive_axis ):
            if mag2( offset ) < kEps: return
            
            ## Dragging axis1 only allows stretching that axis, so
            ## project onto the current axis1.
            
            old_axis1 = ( 1 if positive_axis else -1 ) * self.ellipsoid.axes[1]
            
            offset_axis1 = old_axis1.copy()
            offset_axis1[:2] += asarray( offset )
            
            #scale = dot( offset_axis1, dir( old_axis1 ) ) / mag( old_axis1 )
            mag2_old_axis1 = mag2( old_axis1 )
            scale = dot( offset_axis1, old_axis1 ) / mag2_old_axis1
            ## Don't allow negative or zero scale
            if scale < kEps: return
            
            self.ellipsoid.axes[1] *= scale
            
            ## Don't allow scaling to a length longer than axis0
            if mag2( self.ellipsoid.axes[1] ) > mag2( self.ellipsoid.axes[0] ):
                print("Refusing to stretch Axis1.  Axis0 must remain the long axis.")
                self.ellipsoid.axes[1] *= sqrt( mag2( self.ellipsoid.axes[0] ) / mag2( self.ellipsoid.axes[1] ) )
            
            ## Rescale axes.z appropriately.
            self.SetCrossSectionFromStroke( 1, True, [ self.ellipsoid.center[:2] + self.ellipsoid.axes[2][:2] ], (0,0,-1) )
        
        handle = Struct()
        handle.annotation = self
        handle.kind = self.Axis1Handle
        handle.get_pos = lambda: self.ellipsoid.center + self.ellipsoid.axes[1]
        #handle.translate = lambda offset: pprint( offset )
        handle.translate = lambda offset: drag_axis1( offset, True )
        self._handles.append( handle )
        
        handle = Struct()
        handle.annotation = self
        handle.kind = self.Axis1Handle
        handle.get_pos = lambda: self.ellipsoid.center - self.ellipsoid.axes[1]
        #handle.translate = lambda offset: pprint( offset )
        handle.translate = lambda offset: drag_axis1( offset, False )
        self._handles.append( handle )
        
        
        ## Cross section handle
        def drag_crosssection( offset ):
            if mag2( offset ) < kEps: return
            
            ## Dragging the cross section point means dragging the maximum displacement of the
            ## cross section curve across around axis 1.  This point is also the 2d projection of
            ## axis 2.
            self.SetCrossSectionFromStroke( 1, True, [ self.ellipsoid.center[:2] + self.ellipsoid.axes[2][:2] + asarray( offset ) ], (0,0,-1) )
        
        handle = Struct()
        handle.annotation = self
        handle.kind = self.CrossSectionHandle
        handle.get_pos = lambda: self.ellipsoid.center + self.ellipsoid.axes[2]
        #handle.translate = lambda offset: pprint( offset )
        handle.translate = drag_crosssection
        self._handles.append( handle )
        
        
        ## Parentage intersection handle
        if self.parentage is not None:
            def drag_parent_intersection( offset ):
                assert self.parentage is not None
                
                if mag2( offset ) < kEps: return
                
                ## Dragging the parentage intersection point is easy, since it's an offset anyways.
                self.parentage.center_offset += offset
                
                ## TODO Q: Are we fast enough to recalculate depth placement here?
                ## A1: We don't have the list of all annotations, so the best we can do is
                ##     recalculating our own depth relative to our parent.
                ## A2: All annotations' depth placements are recalculated on mouse-up.
                self.SetDepthRelativeToParent()
            
            def parentage_get_pos():
                guess_pos = self.ellipsoid.center.copy()
                guess_pos[:2] += self.parentage.center_offset
                
                r = ( guess_pos, asarray( (0,0,-1) ) )
                ts = ray_ellipse_intersection_ts( ( self.ellipsoid.center, self.ellipsoid.axes ), r, ray_is_line = True )
                ## If we don't intersect, return our guess
                if len( ts ) == 0: return guess_pos
                ## Also look at the parent's intersection. This matters for back/back intersections,
                ## because we always want to show this handle.
                tsparent = ray_ellipse_intersection_ts( ( self.parentage.parent.ellipsoid.center, self.parentage.parent.ellipsoid.axes ), r, ray_is_line = True )
                
                return r[0] + min( list(ts) + list(tsparent) ) * r[1]
            
            handle = Struct()
            handle.annotation = self
            handle.kind = self.ParentIntersectionHandle
            handle.get_pos = parentage_get_pos
            #handle.translate = lambda offset: pprint( offset )
            handle.translate = drag_parent_intersection
            self._handles.append( handle )
        
        ## These are obsolete
        '''
        ## Discrete Z handle
        class HandleDiscreteZ( object ):
            def __init__( self ):
                self.accum_offset = zerosf(2)
            def translate( self, offset ):
                self.accum_offset += offset
                
                kDiscreteOne = .1
                val = self.accum_offset[1]
                dz = int( sign( val ) * floor( abs( val ) / kDiscreteOne ) )
                self.annotation.z += dz
                self.accum_offset[1] -= dz * kDiscreteOne
        
        handle = HandleDiscreteZ()
        handle.annotation = self
        handle.kind = self.DiscreteZHandle
        #rad = 3*pi/4
        rad = pi/4
        scale = 1.2
        handle.get_pos = lambda: append( ( self.ellipse[0] + cos(rad) * scale * self.ellipse[1] + sin(rad) * scale * self.ellipse[2] ), [0.] )
        #handle.get_pos = lambda: append( ( self.ellipse[0] + self.ellipse[1] + self.ellipse[2] ), [0.] )
        #handle.translate = lambda offset: drag_discretez( handle, offset )
        self._handles.append( handle )
        '''
    
    def GetHandles( self ):
        if self.ellipse is None or self.ellipsoid is None: return []
        if self._handles is None: self._GenHandles()
        return self._handles
    
    Handles = property( GetHandles )
    
    def Iterate( self, N = 1 ):
        return False
    
    def ray_surface_intersection_ts( self, world_ray, ray_is_line = False ):
        '''
        Given a world-space ray (3d point, 3d vector direction), returns the t's along the ray
        corresponding to intersections with the shape's 3d surface.
        If 'ray_is_line' is False, only intersections in the positive direction along the ray
        will be returned.
        
        untested
        '''
        
        return ray_ellipse_intersection_ts(
            ( self.ellipsoid.center, self.ellipsoid.axes ),
            world_ray,
            ray_is_line
            )
            
    
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
    
    def candidate_origins( self ):
        '''
        Returns a sequence of possible origin points, namely the ellipsoid's center.
        '''
        
        return [ self.ellipsoid.center ]
    
    def center( self ):
        '''
        Returns the center of the shape.
        '''
        
        return self.ellipsoid.center
    
    def translate( self, offset ):
        '''
        Given a 2d or 3d offset vector 'offset', translate this shape.
        '''
        
        assert len( offset ) in (2,3)
        
        offset3d = zerosf(3)
        offset3d[:len(offset)] = offset
        
        self.ellipsoid.center += offset3d
        self.ellipse = list( self.ellipse )
        self.ellipse[0] = asarray( self.ellipse[0] ) + offset3d[:2]
        self.ellipse = tuple( self.ellipse )
    
    def offset_in_depth( self, z_offset, oblique = (0.,0.) ):
        '''
        Translate the object along the look vector (oblique[0], oblique[1], -1) by 'z_offset'.
        Inverts the oblique projection, so the object maintains the same (projected) appearance.
        '''
        
        self.ellipsoid.center += z_offset * array( ( oblique[0], oblique[1], -1 ) )

def main():
    pass

if __name__ == '__main__': main()
