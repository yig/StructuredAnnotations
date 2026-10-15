#!/usr/bin/env python

import sys
from GLUTWindow import GLUTWindow
from myarray import *
from vecutil import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.

try:
    from OpenGL.GLUT import *
    from OpenGL.GL import *
    from OpenGL.GLU import *
except:
    print('''
ERROR: PyOpenGL not installed properly. 
        ''')
    sys.exit()

class TriMesh( object ):
    def __init__( self ):
        self.vs = []
        self.faces = []
        
        self.__face_normals = None
        self.__face_areas = None
        self.__vertex_normals = None
        self.__vertex_areas = None
        self.__edges = None
        
        self.__halfedges = None
        self.__vertex_halfedges = None
        self.__face_halfedges = None
        self.__edge_halfedges = None
        self.__directed_edge2he_index = None
        
        self.lifetime_counter = 0
    
    def __deepcopy__( self, memodict ):
        result = TriMesh()
        
        ## Make a copy of vs and faces using array().
        ## But in case they weren't stored as arrays, return them as the type they were.
        ## This is important if they were lists, and someone expected to be able to call
        ## .append() or .extend() on them.
        result.vs = array( self.vs )
        if not isinstance( self.vs, ndarray ):
            result.vs = type( self.vs )( result.vs )
        
        result.faces = array( self.faces )
        if not isinstance( self.faces, ndarray ):
            result.faces = type( self.faces )( result.faces )
        
        if hasattr( self, 'uvs' ):
            result.uvs = array( self.uvs )
            if not isinstance( self.uvs, ndarray ):
                result.uvs = type( self.uvs )( result.uvs )
        
        ## I could skip copying these cached values, but they are usually needed for rendering
        ## and copy quickly.
        if self.__face_normals is not None:
            result.__face_normals = self.__face_normals.copy()
        if self.__face_areas is not None:
            result.__face_areas = self.__face_areas.copy()
        if self.__vertex_normals is not None:
            result.__vertex_normals = self.__vertex_normals.copy()
        if self.__vertex_areas is not None:
            result.__vertex_areas = self.__vertex_areas.copy()
        if self.__edges is not None:
            result.__edges = list( self.__edges )
        
        ## I will skip copying these cached values, because they copy slowly and are
        ## not as commonly needed.  They'll still be regenerated as needed.
        '''
        if self.__halfedges is not None:
            from copy import copy
            result.__halfedges = [ copy( he ) for he in self.__halfedges ]
        if self.__vertex_halfedges is not None:
            result.__vertex_halfedges = list( self.__vertex_halfedges )
        if self.__face_halfedges is not None:
            result.__face_halfedges = list( self.__face_halfedges )
        if self.__edge_halfedges is not None:
            result.__edge_halfedges = list( self.__edge_halfedges )
        if self.__directed_edge2he_index is not None:
            result.__directed_edge2he_index = dict( self.__directed_edge2he_index )
        '''
        
        result.lifetime_counter = self.lifetime_counter
        
        return result
    
    def update_face_normals_and_areas( self ):
        if self.__face_normals is None: self.__face_normals = zerosf( ( len( self.faces ), 3 ) )
        if self.__face_areas is None: self.__face_areas = zerosf( len( self.faces ) )
        
        if len( self.vs ) > 0 and len( self.faces ) > 0:
            ## If we have this module, this is 3-4 times faster.
            try:
                import fibermesh_cpp
                om = fibermesh_cpp.new_openmesh_from_data( self.vs, self.faces )
                self.__face_normals = fibermesh_cpp.extract_openmesh_face_normals( om )
                self.__face_areas = fibermesh_cpp.extract_openmesh_face_areas( om )
                fibermesh_cpp.delete_openmesh( om )
                return
            except ImportError:
                pass
        
        ## We need subtraction between vertices.
        ## Convert vertices to arrays once here, or else we'd have to call asarray()
        ## ~6 times for each vertex.
        ## NOTE: If self.vs is already an array, then this code does nothing.
        ## TODO Q: Should I set self.vs = asarray( self.vs )?  It might violate someone's
        ##         assumption that self.vs is whatever indexable type they left it.
        ##         In particular, this violates the ability of someone to .append() or .extend()
        ##         self.vs.
        vs = asarray( self.vs )
        
        for f in range( len( self.faces ) ):
            face = self.faces[f]
            n = cross(
                vs[ face[1] ] - vs[ face[0] ],
                vs[ face[2] ] - vs[ face[1] ]
                )
            nmag = mag( n )
            self.__face_normals[f] = (1./nmag) * n
            self.__face_areas[f] = .5 * nmag
        
        assert len( self.faces ) == len( self.__face_normals )
        assert len( self.faces ) == len( self.__face_areas )
    
    
    def get_face_normals( self ):
        if self.__face_normals is None: self.update_face_normals_and_areas()
        return self.__face_normals
    
    face_normals = property( get_face_normals )
    
    
    def get_face_areas( self ):
        if self.__face_areas is None: self.update_face_normals_and_areas()
        return self.__face_areas
    
    face_areas = property( get_face_areas )
    
    
    def update_vertex_normals( self ):
        if self.__vertex_normals is None: self.__vertex_normals = zerosf( ( len(self.vs), 3 ) )
        
        ## If we have this module, this is so much faster than generating
        ## the whole halfedge data structure in python.
        try:
            import fibermesh_cpp
            om = fibermesh_cpp.new_openmesh_from_data( self.vs, self.faces )
            self.__vertex_normals = fibermesh_cpp.extract_openmesh_vertex_normals( om )
            fibermesh_cpp.delete_openmesh( om )
            return
        except ImportError:
            pass
        
        for vi in range( len( self.vs ) ):
            self.__vertex_normals[vi] = 0.
            
            for fi in self.vertex_face_neighbors( vi ):
                ## This matches the OpenMesh's FAST vertex normals.
                #self.__vertex_normals[vi] += self.face_normals[ fi ]
                ## Area weighted
                self.__vertex_normals[vi] += self.face_normals[ fi ] * self.face_areas[ fi ]
            
        ## Now normalize the normals
        #self.__vertex_normals[vi] *= 1./mag( self.__vertex_normals[vi] )
        self.__vertex_normals *= 1./sqrt( ( self.__vertex_normals**2 ).sum(1) ).reshape( (len(self.vs), 1) )
        
        assert len( self.vs ) == len( self.__vertex_normals )
    
    def get_vertex_normals( self ):
        if self.__vertex_normals is None: self.update_vertex_normals()
        return self.__vertex_normals
    
    vertex_normals = property( get_vertex_normals )
    
    
    def update_vertex_areas( self ):
        if self.__vertex_areas is None: self.__vertex_areas = zerosf( len(self.vs) )
        
        for vi in range( len( self.vs ) ):
            ## Try to compute proper area (if we have laplacian editing around)
            try:
                #raise ImportError
                import laplacian_editing
                cot_alpha, cot_beta, area = laplacian_editing.cotangentWeights(
                    self.vs[ vi ],
                    [ self.vs[ vni ] for vni in self.vertex_vertex_neighbors( vi ) ],
                    self.vertex_is_boundary( vi )
                    )
                self.__vertex_areas[vi] = area
            
            ## Otherwise use 1/3 of the incident faces' areas
            except ImportError:
                self.__vertex_areas[vi] = 0.
                for fi in self.vertex_face_neighbors( vi ):
                    self.__vertex_areas[vi] += self.face_areas[ fi ]
                
                self.__vertex_areas[vi] *= 1./3.
        
        assert len( self.vs ) == len( self.__vertex_areas )
    
    def get_vertex_areas( self ):
        if self.__vertex_areas is None: self.update_vertex_areas()
        return self.__vertex_areas
    
    vertex_areas = property( get_vertex_areas )
    
    
    def update_edge_list( self ):
        ## We need a set of set-pairs of vertices, because edges are bidirectional.
        edges = set()
        for face in self.faces:
            edges.add( frozenset( ( face[0], face[1] ) ) )
            edges.add( frozenset( ( face[1], face[2] ) ) )
            edges.add( frozenset( ( face[2], face[0] ) ) )
        
        self.__edges = [ tuple( edge ) for edge in edges ]
    
    def get_edges( self ):
        if self.__edges is None: self.update_edge_list()
        return self.__edges
    
    edges = property( get_edges )
    
    
    class HalfEdge( object ):
        def __init__( self ):
            self.to_vertex = -1
            self.face = -1
            self.edge = -1
            self.opposite_he = -1
            self.next_he = -1
    
    def update_halfedges( self ):
        '''
        Generates all half edge data structures for the mesh given by its vertices 'self.vs'
        and faces 'self.faces'.
        
        untested
        '''
        
        self.__halfedges = []
        self.__vertex_halfedges = None
        self.__face_halfedges = None
        self.__edge_halfedges = None
        self.__directed_edge2he_index = {}
        
        __directed_edge2face_index = {}
        for fi, face in enumerate( self.faces ):
            __directed_edge2face_index[ (face[0], face[1]) ] = fi
            __directed_edge2face_index[ (face[1], face[2]) ] = fi
            __directed_edge2face_index[ (face[2], face[0]) ] = fi
        
        def directed_edge2face_index( edge ):
            result = __directed_edge2face_index.get( edge, -1 )
            
            ## If result is -1, then there's no such face in the mesh.
            ## The edge must be a boundary edge.
            ## In this case, the reverse orientation edge must be have a face.
            if -1 == result:
                assert edge[::-1] in __directed_edge2face_index
            
            return result
        
        self.__vertex_halfedges = [None] * len( self.vs )
        self.__face_halfedges = [None] * len( self.faces )
        self.__edge_halfedges = [None] * len( self.edges )
        
        for ei, edge in enumerate( self.edges ):
            he0 = self.HalfEdge()
            ## The face will be -1 if it is a boundary half-edge.
            he0.face = directed_edge2face_index( edge )
            he0.to_vertex = edge[1]
            he0.edge = ei
            
            he1 = self.HalfEdge()
            ## The face will be -1 if it is a boundary half-edge.
            he1.face = directed_edge2face_index( edge[::-1] )
            he1.to_vertex = edge[0]
            he1.edge = ei
            
            ## Add the HalfEdge structures to the list.
            he0index = len( self.__halfedges )
            self.__halfedges.append( he0 )
            he1index = len( self.__halfedges )
            self.__halfedges.append( he1 )
            
            ## Now we can store the opposite half-edge index.
            he0.opposite_he = he1index
            he1.opposite_he = he0index
            
            ## Also store the index in our __directed_edge2he_index map.
            assert edge not in self.__directed_edge2he_index
            assert edge[::-1] not in self.__directed_edge2he_index
            self.__directed_edge2he_index[ edge ] = he0index
            self.__directed_edge2he_index[ edge[::-1] ] = he1index
            
            ## If the vertex pointed to by a half-edge doesn't yet have an out-going
            ## halfedge, store the opposite halfedge.
            ## Also, if the vertex is a boundary vertex, make sure its
            ## out-going halfedge a boundary halfedge.
            ## NOTE: Halfedge data structure can't properly handle butterfly vertices.
            ##       If the mesh has butterfly vertices, there will be multiple outgoing
            ##       boundary halfedges.  Because we have to pick one as the vertex's outgoing
            ##       halfedge, we can't iterate over all neighbors, only a single wing of the
            ##       butterfly.
            if self.__vertex_halfedges[ he0.to_vertex ] is None or -1 == he1.face:
                self.__vertex_halfedges[ he0.to_vertex ] = he0.opposite_he
            if self.__vertex_halfedges[ he1.to_vertex ] is None or -1 == he0.face:
                self.__vertex_halfedges[ he1.to_vertex ] = he1.opposite_he
            
            ## If the face pointed to by a half-edge doesn't yet have a
            ## halfedge pointing to it, store the halfedge.
            if -1 != he0.face and self.__face_halfedges[ he0.face ] is None:
                self.__face_halfedges[ he0.face ] = he0index
            if -1 != he1.face and self.__face_halfedges[ he1.face ] is None:
                self.__face_halfedges[ he1.face ] = he1index
            
            ## Store one of the half-edges for the edge.
            assert self.__edge_halfedges[ ei ] is None
            self.__edge_halfedges[ ei ] = he0index
        
        ## Now that all the half-edges are created, set the remaining next_he field.
        ## We can't yet handle boundary halfedges, so store them for later.
        boundary_heis = []
        for hei, he in enumerate( self.__halfedges ):
            ## Store boundary halfedges for later.
            if -1 == he.face:
                boundary_heis.append( hei )
                continue
            
            face = self.faces[ he.face ]
            i = he.to_vertex
            j = face[ ( list(face).index( i ) + 1 ) % 3 ]
            
            he.next_he = self.__directed_edge2he_index[ (i,j) ]
        
        ## Make a map from vertices to boundary halfedges (indices) originating from them.
        ## NOTE: There will only be multiple originating boundary halfedges at butterfly vertices.
        vertex2outgoing_boundary_hei = {}
        for hei in boundary_heis:
            originating_vertex = self.__halfedges[ self.__halfedges[ hei ].opposite_he ].to_vertex
            vertex2outgoing_boundary_hei.setdefault(
                originating_vertex, set()
                ).add( hei )
            if len( vertex2outgoing_boundary_hei[ originating_vertex ] ) > 1:
                print('Butterfly vertex encountered')
        
        ## For each boundary halfedge, make its next_he one of the boundary halfedges
        ## originating at its to_vertex.
        for hei in boundary_heis:
            he = self.__halfedges[ hei ]
            for outgoing_hei in vertex2outgoing_boundary_hei[ he.to_vertex ]:
                he.next_he = outgoing_hei
                vertex2outgoing_boundary_hei[ he.to_vertex ].remove( outgoing_hei )
                break
        
        assert False not in [ 0 == len( out_heis ) for out_heis in vertex2outgoing_boundary_hei.values() ]
    
    def he_index2directed_edge( self, he_index ):
        '''
        Given the index of a HalfEdge, returns the corresponding directed edge (i,j).
        
        untested
        '''
        
        he = self.halfedges[ he_index ]
        return ( self.halfedges[ he.opposite_he ].to_vertex, he.to_vertex )
    
    def directed_edge2he_index( self, edge ):
        '''
        Given a directed edge (i,j), returns the index of the HalfEdge class in
        halfedges().
        
        untested
        '''
        
        if self.__directed_edge2he_index is None: self.update_halfedges()
        
        edge = tuple( edge )
        return self.__directed_edge2he_index[ edge ]
    
    def get_halfedges( self ):
        '''
        Returns a list of all HalfEdge classes.
        
        untested
        '''
        
        if self.__halfedges is None: self.update_halfedges()
        return self.__halfedges
    
    halfedges = property( get_halfedges )
    
    def vertex_vertex_neighbors( self, vertex_index ):
        '''
        Returns the vertex neighbors (as indices) of the vertex 'vertex_index'.
        
        untested
        '''
        
        ## It's important to access self.halfedges first (which calls get_halfedges()),
        ## so that we're sure all halfedge info is generated.
        halfedges = self.halfedges
        result = []
        start_he = halfedges[ self.__vertex_halfedges[ vertex_index ] ]
        he = start_he
        while True:
            result.append( he.to_vertex )
            
            he = halfedges[ halfedges[ he.opposite_he ].next_he ]
            if he is start_he: break
        
        return result
    
    def vertex_valence( self, vertex_index ):
        '''
        Returns the valence (number of vertex neighbors) of vertex with index 'vertex_index'.
        
        untested
        '''
        
        return len( self.vertex_vertex_neighbors( vertex_index ) )
    
    def vertex_face_neighbors( self, vertex_index ):
        '''
        Returns the face neighbors (as indices) of the vertex 'vertex_index'.
        
        untested
        '''
        
        ## It's important to access self.halfedges first (which calls get_halfedges()),
        ## so that we're sure all halfedge info is generated.
        halfedges = self.halfedges
        result = []
        start_he = halfedges[ self.__vertex_halfedges[ vertex_index ] ]
        he = start_he
        while True:
            if -1 != he.face: result.append( he.face )
            
            he = halfedges[ halfedges[ he.opposite_he ].next_he ]
            if he is start_he: break
        
        return result
    
    def vertex_is_boundary( self, vertex_index ):
        '''
        Returns whether the vertex with given index is on the boundary.
        
        untested
        '''
        
        ## It's important to access self.halfedges first (which calls get_halfedges()),
        ## so that we're sure all halfedge info is generated.
        halfedges = self.halfedges
        return -1 == halfedges[ self.__vertex_halfedges[ vertex_index ] ].face
    
    def boundary_vertices( self ):
        '''
        Returns a list of the vertex indices on the boundary.
        
        untested
        '''
        
        result = []
        for hei, he in enumerate( self.halfedges ):
            if -1 == he.face:
                # result.extend( self.he_index2directed_edge( hei ) )
                result.append( he.to_vertex )
                result.append( self.halfedges[ he.opposite_he ].to_vertex )
        
        return list(frozenset( result ))
    
    def boundary_edges( self ):
        '''
        Returns a list of boundary edges (i,j).  If (i,j) is in the result, (j,i) will not be.
        
        untested
        '''
        
        result = []
        for hei, he in enumerate( self.halfedges ):
            if -1 == he.face:
                result.append( self.he_index2directed_edge( hei ) )
        return result
    
    def positions_changed( self ):
        '''
        Notify the object that vertex positions changed.
        All position-related structures (normals, areas) will be marked for re-calculation.
        '''
        
        self.__face_normals = None
        self.__face_areas = None
        self.__vertex_normals = None
        self.__vertex_areas = None
        
        self.lifetime_counter += 1
    
    
    def topology_changed( self ):
        '''
        Notify the object that topology (faces or #vertices) changed.
        All topology-related structures (halfedges, edge lists) as well as position-related
        structures (normals, areas) will be marked for re-calculation.
        '''
        
        ## Set mesh.vs to an array so that subsequent calls to asarray() on it are no-ops.
        self.vs = asarray( self.vs )
        
        self.__edges = None
        self.__halfedges = None
        self.__vertex_halfedges = None
        self.__face_halfedges = None
        self.__edge_halfedges = None
        self.__directed_edge2he_index = None
        
        self.positions_changed()
    
    def get_dangling_vertices( self ):
        '''
        Returns vertex indices in TriMesh 'mesh' that belong to no faces.
        '''
        brute_vertex_face_valence = [ 0 ] * len( self.vs )
        for i,j,k in self.faces:
            brute_vertex_face_valence[ i ] += 1
            brute_vertex_face_valence[ j ] += 1
            brute_vertex_face_valence[ k ] += 1
        return [ i for i in range( len( self.vs ) ) if 0 == brute_vertex_face_valence[i] ]
    
    def remove_vertex_indices( self, vertex_indices_to_remove ):
        '''
        Removes vertices in the list of indices 'vertex_indices_to_remove'.
        Also removes faces containing the vertices.
        
        used
        '''
        
        ## I can't assert this here because I call this function recursively to remove dangling
        ## vertices.
        ## Also, someone manipulating the mesh might want to do the same thing (call this
        ## function on dangling vertices).
        #assert 0 == len( self.get_dangling_vertices() )
        
        
        if 0 == len( vertex_indices_to_remove ): return list(range( len( self.vs )))
        
        
        ## Make a map from old to new vertices.  This is the return value.
        old2new = [ -1 ] * len( self.vs )
        last_index = 0
        for i in range( len( self.vs ) ):
            if i not in vertex_indices_to_remove:
                old2new[ i ] = last_index
                last_index += 1
        
        ## Remove vertices from vs, faces, edges.
        self.vs = [ pt for i, pt in enumerate( self.vs ) if old2new[i] != -1 ]
        #self.__vertex_normals = asarrayf( [ vn for i, vn in enumerate( self.__vertex_normals ) if old2new[i] != -1 ]
        #self.__edges = [ ( old2new[i], old2new[j] ) for i,j in self.__edges ]
        #self.__edges = [ edge for edge in self.__edges if -1 not in edge ]
        self.faces = [ ( old2new[i], old2new[j], old2new[k] ) for i,j,k in self.faces ]
        #self.__face_normals = [ n for i,n in enumerate( self.__face_normals ) if -1 not in self.faces[i] ]
        #self.__face_areas = [ n for i,n in enumerate( self.__face_areas ) if -1 not in self.faces[i] ]
        self.faces = [ tri for tri in self.faces if -1 not in tri ]
        
        ## Now that we have halfedge info, just call topology changed and everything but
        ## 'vs' and 'faces' will be regenerated.
        self.topology_changed()
        
        ## Remove dangling vertices created by removing faces incident to vertices in 'vertex_indices_to_remove'.
        ## We only need to call this once, because a dangling vertex has no faces, so its removal
        ## won't remove any faces, so no new dangling vertices can be created.
        old2new_recurse = self.remove_vertex_indices( self.get_dangling_vertices() )
        assert 0 == len( self.get_dangling_vertices() )
        
        for i in range( len( old2new ) ):
            if -1 != old2new[i]: old2new[i] = old2new_recurse[ old2new[ i ] ]
        
        return old2new
    
    
    def append( self, mesh ):
        '''
        Given a mesh, with two properties,
            .vs, containing a list of 3d vertices
            .faces, containing a list of triangles as triplets of indices into .vs
        appends 'mesh's vertices and faces to self.vs and self.faces.
        '''
        
        ## mesh's vertices are going to be copied to the end of self.vs;
        ## All vertex indices in mesh.faces will need to be offset by the current
        ## number of vertices in self.vs.
        vertex_offset = len( self.vs )
        
        self.vs = list( self.vs ) + list( mesh.vs )
        self.faces = list( self.faces ) + list( asarray( mesh.faces, dtype = int ) + vertex_offset )
        
        
        ## If there are uvs, concatenate them.
        
        ## First, if self is an empty mesh (without uv's), and the mesh to append-to has uv's,
        ## create an empty .uvs property in self.
        if not hasattr( self, 'uvs' ) and hasattr( mesh, 'uvs' ) and len( self.vs ) == 0:
            self.uvs = []
        
        if hasattr( self, 'uvs' ) and hasattr( mesh, 'uvs' ):
            self.uvs = list( self.uvs ) + list( mesh.uvs )
        elif hasattr( self, 'uvs' ):
            del self.uvs
        
        
        ## We're almost done, we only need to call topology_changed().
        ## However, let's see if we can keep some properties that are slow to regenerate.
        self__face_normals = self.__face_normals
        self__face_areas = self.__face_areas
        self__vertex_normals = self.__vertex_normals
        self__vertex_areas = self.__vertex_areas
        
        self.topology_changed()
        
        if self__face_normals is not None and mesh.__face_normals is not None:
            self.__face_normals = append( self__face_normals, mesh.__face_normals, axis = 0 )
        if self__face_areas is not None and mesh.__face_areas is not None:
            self.__face_areas = append( self__face_areas, mesh.__face_areas, axis = 0 )
        if self__vertex_normals is not None and mesh.__vertex_normals is not None:
            self.__vertex_normals = append( self__vertex_normals, mesh.__vertex_normals, axis = 0 )
        if self__vertex_areas is not None and mesh.__vertex_areas is not None:
            self.__vertex_areas = append( self__vertex_areas, mesh.__vertex_areas, axis = 0 )
    
    
    def FromTriMeshes( meshes ):
        '''
        Given a sequence of meshes, each with two properties,
            .vs, containing a list of 3d vertices
            .faces, containing a list of triangles as triplets of indices into .vs
        returns a single TriMesh object containing all meshes concatenated together.
        '''
        
        result = TriMesh()
        for mesh in meshes:
            result.append( mesh )
        
        ## Reset the lifetime counter
        result.lifetime_counter = 0
        return result
    
    FromTriMeshes = staticmethod( FromTriMeshes )
    
    
    def FromOBJ_FileName( obj_fname ):
        return TriMesh.FromOBJ_Lines( open( obj_fname ) )
    
    FromOBJ_FileName = staticmethod( FromOBJ_FileName )
    
    
    def FromOBJ_Lines( obj_lines ):
        '''
        Given lines from an OBJ file, return a new TriMesh object.
        
        tested
        '''
        
        result = TriMesh()
        
        ## NOTE: We only handle faces and vertex positions.
        for line in obj_lines:
            line = line.strip()
            
            sline = line.split()
            ## Skip blank lines
            if not sline: continue
            
            elif sline[0] == 'v':
                result.vs.append( asarrayf( list(map( float, sline[1:] )) ) )
                ## Vertices must have three coordinates.
                assert len( result.vs[-1] ) == 3
            
            elif sline[0] == 'f':
                ## The split('/')[0] means we record only the vertex coordinate indices
                ## for each face.
                face_vertex_ids = [ int( c.split('/')[0] ) for c in sline[1:] ]
                ## Faces must be triangles.
                assert len( face_vertex_ids ) == 3
                
                ## Face vertex indices cannot be zero.
                assert True not in [ ind == 0 for ind in face_vertex_ids ]
                
                ## Subtract one from positive indices, and use relative addressing for negative
                ## indices.
                face_vertex_ids = [
                    ## This awkward "make a tuple with both results and select it based on
                    ## truth" technique is to make up for Python's lack of a functional 'if'.
                    ( ind-1, len(result.vs) + ind )[ ind < 0 ]
                    for ind in face_vertex_ids
                    ]
                
                assert False not in [ ind < len( result.vs ) for ind in face_vertex_ids ]
                result.faces.append( face_vertex_ids )
        
        return result
    
    FromOBJ_Lines = staticmethod( FromOBJ_Lines )
    
    
    def write_OBJ( self, fname ):
        '''
        Writes the data out to an OBJ file named 'fname'.
        
        tested
        '''
        
        import sys
        
        out = open( fname, 'w' )
        out.write( '## Written by %s\n\n' % (sys.argv[0],) )
        
        for v in self.vs:
            out.write( 'v %s %s %s\n' % tuple(v) )
        out.write( '\n' )
        
        if hasattr( self, 'uvs' ):
            for uv in self.uvs:
                out.write( 'vt %s %s\n' % tuple(uv) )
            out.write( '\n' )
            
            for f in self.faces:
                out.write( 'f %s/%s %s/%s %s/%s\n' % tuple( ( asarray(f) + 1 ).repeat(2) ) )
        else:
            for f in self.faces:
                out.write( 'f %s %s %s\n' % tuple(asarray(f) + 1) )
        
        print('[OBJ written to "%s"]' % (fname,))

## We can't pickle anything that doesn't have a name visible at module scope.
## In order to allow pickling of class TriMesh, we'll make a reference to the inner HalfEdge class
## here at the module level.
HalfEdge = TriMesh.HalfEdge

class Camera( object ):
    def __init__( self ):
        self.center = arrayf( ( 0,0,0 ) )
        self.eye = arrayf( ( 0,0,2 ) )
        self.up = arrayf( ( 0,1,0 ) )
        self.near_far = asarrayf( ( .1, 4.1 ) )
        
        self.proj = self.Perspective( fovy = 120., aspect = 1. )
        #self.proj = self.Ortho( -1., 1., -1., 1. )
    
    class Perspective( object ):
        def __init__( self, **kwargs ):
            kwargs.setdefault( 'fovy', 120. )
            kwargs.setdefault( 'aspect', 1. )
            
            self.fovy = kwargs[ 'fovy' ]
            self.aspect = kwargs[ 'aspect' ]
        
        def apply( self, view_width, view_height, near, far ):
            gluPerspective(
                self.fovy, self.aspect * float( view_width ) / float( view_height ),
                near, far
                )
    
    class Ortho( object ):
        def __init__( self, left = -1., right = 1., bottom = -1., top = 1. ):
            self.left = left
            self.right = right
            self.bottom = bottom
            self.top = top
        
        def apply( self, view_width, view_height, near, far ):
            glOrtho( self.left, self.right, self.bottom, self.top, near, far )
    
    class Oblique( object ):
        def __init__( self, left = -1., right = 1., bottom = -1., top = 1., a = 0., b = 0., eye_distance_to_zero_plane = 0. ):
            self.left = left
            self.right = right
            self.bottom = bottom
            self.top = top
            self.a = a
            self.b = b
            self.z0 = eye_distance_to_zero_plane
        
        def apply( self, view_width, view_height, near, far ):
            glOrtho( self.left, self.right, self.bottom, self.top, near, far )
            shear_matrix = identity( 4 )
            shear_matrix[0,2] = self.a
            shear_matrix[1,2] = self.b
            #print 'shear_matrix:', shear_matrix
            glTranslated( 0, 0, -self.z0 )
            glMultMatrixd( shear_matrix.T.flatten() )
            glTranslated( 0, 0, self.z0 )
    
    def Orbit( self, axis, radians ):
        '''
        Orbits the camera about the view center using the rotation defined by 'axis' and 'radians'.
        '''
        
        self.eye = self.center + rotate( self.eye - self.center, axis, radians )
        self.up = rotate( self.up, axis, radians )
    
    def Rotate2( self, left_right, up_down ):
        '''
        Rotates the camera given image plane motion 'left_right' and 'up_down',
        scalars in the range -1..1.
        '''
        
        #print 'left_right:', left_right
        #print 'up_down:', up_down
        coef = pi
        
        self.Orbit( self.up, coef * left_right )
        self.Orbit( cross( self.up, self.eye - self.center ), coef * -up_down )
    
    def Apply( self ):
        self.__ApplyProjection()
        self.__ApplyModelview()
    
    def __ApplyProjection( self ):
        glMatrixMode( GL_PROJECTION )
        glLoadIdentity()
        
        view = glGetIntegerv( GL_VIEWPORT )
        view_width = view[2]
        view_height = view[3]
        del view
        
        self.proj.apply( view_width, view_height, self.near_far[0], self.near_far[1] )
    
    def __ApplyModelview( self ):
        glMatrixMode( GL_MODELVIEW )
        glLoadIdentity()
        
        gluLookAt(
            self.eye[0], self.eye[1], self.eye[2],
            self.center[0], self.center[1], self.center[2],
            self.up[0], self.up[1], self.up[2]
            )

def set_headlight( camera ):
    light_position = [ camera.eye[0], camera.eye[1], camera.eye[2], 1.0 ]
    glLightfv( GL_LIGHT0, GL_POSITION, light_position )
    glEnable( GL_LIGHT0 )

def draw_mesh_faces_flat( mesh, will_draw_edges = True ):
    '''
    Takes a TriMesh parameter mesh and draws it, setting as little OpenGL state as possible.
    Lighting and materials and colors be specified before this function is entered.
    If the parameters indicates that edges will be drawn, then glPolygonOffset will be used.
    '''
    
    if will_draw_edges:
        ## This offset guarantees that all polygon faces have a z-buffer value at least 1 greater.
        glPolygonOffset( 1.0, 1.0 )
        glEnable( GL_POLYGON_OFFSET_FILL )
    
    
    glBegin( GL_TRIANGLES )
    for face, normal in zip( mesh.faces, mesh.face_normals ):
        glNormal3f( *normal )
        for vertex_index in face:
            glVertex3f( *mesh.vs[ vertex_index ] )
    glEnd()
    
    glDisable( GL_POLYGON_OFFSET_FILL )

dtype_type2glenum = {
    float32: GL_FLOAT,
    float64: GL_DOUBLE,
    int32: GL_INT,
    uint32: GL_UNSIGNED_INT
    }
def draw_mesh_faces_smooth( mesh, will_draw_edges = True ):
    '''
    Takes a TriMesh parameter mesh and draws it, setting as little OpenGL state as possible.
    Lighting and materials and colors be specified before this function is entered.
    If the parameters indicates that edges will be drawn, then glPolygonOffset will be used.
    '''
    
    ## I assume that mesh.vs[0] exists, so I need this.
    if len( mesh.vs ) == 0 or len( mesh.faces ) == 0: return
    assert len( mesh.vs ) == len( mesh.vertex_normals )
    
    if will_draw_edges:
        ## This offset guarantees that all polygon faces have a z-buffer value at least 1 greater.
        glPolygonOffset( 1.0, 1.0 )
        glEnable( GL_POLYGON_OFFSET_FILL )
    
    glEnableClientState( GL_VERTEX_ARRAY )
    ## Convert mesh.vs to an array() so that subsequent calls to this function are fast.
    mesh.vs = asarray( mesh.vs )
    mesh_vs = mesh.vs # asarray( mesh.vs )
    glVertexPointer( mesh_vs.shape[1], dtype_type2glenum[ mesh_vs.dtype.type ], 0, mesh_vs )
    ## NOTE: PyOpenGL's convenience glVertexPointerd(), glNormalPointerd(), and glDrawElementsui()
    ##       are painfully slow!  Why is that?  Why doesn't it just call asarray() and check the
    ##       resulting shape and dtype?
    #glVertexPointerd( mesh.vs )
    
    glEnableClientState( GL_NORMAL_ARRAY )
    ## mesh.vertex_normals is always an array, and we can't set the property to asarray() anyways.
    #mesh.vertex_normals = asarray( mesh.vertex_normals )
    mesh_normals = mesh.vertex_normals # asarray( mesh.vertex_normals )
    glNormalPointer( dtype_type2glenum[ mesh_normals.dtype.type ], 0, mesh_normals )
    #glNormalPointerd( mesh.vertex_normals )
    
    ## We can't do the dtype_type2glenum[] thing here because it must be an unsigned int.
    glDrawElements( GL_TRIANGLES, 3 * len( mesh.faces ), GL_UNSIGNED_INT, asarray( mesh.faces, dtype = uint32 ) )
    #glDrawElementsui( GL_TRIANGLES, mesh.faces )
    
    glDisableClientState( GL_NORMAL_ARRAY )
    glDisableClientState( GL_VERTEX_ARRAY )
    
    glDisable( GL_POLYGON_OFFSET_FILL )

## The default is flat shading.
draw_mesh_faces = draw_mesh_faces_flat

def draw_mesh_edges( mesh ):
    '''
    Takes a TriMesh parameter mesh and draws it, setting as little OpenGL state as possible.
    '''
    
    glBegin( GL_LINES )
    for edge in mesh.edges:
        edge = tuple( edge )
        glVertex3f( *mesh.vs[ edge[0] ] )
        glVertex3f( *mesh.vs[ edge[1] ] )
    glEnd()

class TriMeshWindow( GLUTWindow ):
    def __init__( self, **kwargs ):
        GLUTWindow.__init__( self )
        
        kwargs.setdefault( 'background_color', ( .3, .3, .3 ) )
        kwargs.setdefault( 'mesh', TriMesh() )
        kwargs.setdefault( 'camera', Camera() )
        kwargs.setdefault( 'draw_faces', True )
        kwargs.setdefault( 'draw_edges', True )
        
        self.background_color = kwargs[ 'background_color' ]
        self.mesh = kwargs[ 'mesh' ]
        self.camera = kwargs[ 'camera' ]
        self.draw_faces = kwargs[ 'draw_faces' ]
        self.draw_edges = kwargs[ 'draw_edges' ]
        
        ## I don't want dynamic binding here because
        ## subclasses' constructors haven't run yet
        TriMeshWindow.reset( self )
    
    def reset( self ):
        print('Resetting')
        
        self.postRedisplay()
    
    def displayFunc( self ):
        col = self.background_color
        glClearColor( col[0], col[1], col[2], 1.0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        glEnable( GL_DEPTH_TEST )
        
        self.camera.Apply()
        ## TODO: Something with the light.  Right now, we're getting OpenGL's default light.
        
        set_headlight( self.camera )
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE )
        glShadeModel( GL_SMOOTH )
        
        if self.draw_faces:
            glEnable( GL_LIGHTING )
            
            glEnable( GL_COLOR_MATERIAL )
            glColor3ub( 200, 200, 255 )
            draw_mesh_faces( self.mesh, self.draw_edges )
        
        if self.draw_edges:
            glDisable( GL_LIGHTING )
            
            glLineWidth( 1 )
            glColor3ub( 0,0,0 )
            draw_mesh_edges( self.mesh )
    
    def motionFunc( self, x, y ):
        ## Very simply map the camera onto +z half of the sphere x^2 + y^2 + z^2 = 1
        
        view = glGetIntegerv( GL_VIEWPORT )
        assert view[0] == 0
        assert view[1] == 0
        
        wx = view[2]/2.
        wy = view[3]/2.
        
        ## Flip y so that 0,0 is in the bottom-left.
        y = view[3] - y - 1
        
        radius = min( view[2], view[3] )/2.
        sphere_x = ( x - wx ) / radius
        sphere_y = ( y - wy ) / radius
        
        ## Project onto the closest point on the sphere
        if sphere_x*sphere_x + sphere_y*sphere_y > 1:
            denom = sqrt( sphere_x*sphere_x + sphere_y*sphere_y )
            sphere_x /= denom
            sphere_y /= denom
        
        sphere_z = sqrt( max( 0, 1 - sphere_x*sphere_x - sphere_y*sphere_y ) )
        
        eye_mag = mag( self.camera.eye )
        self.camera.eye = eye_mag * dir( arrayf( ( sphere_x, sphere_y, sphere_z ) ) )
        # print 'New eye position: ', self.camera.eye
        
        glutPostRedisplay()
    
    def mouseFunc( self, button, state, x,y ):
        if state == GLUT_DOWN: pass
        pass
    
    def passiveMotionFunc( self, x,y ):
        pass
    
    def keyboardFunc( self, key, x, y ):
        if 'r' == key: self.reset()
        elif 'w' == key: self.draw_edges = not self.draw_edges
        else: GLUTWindow.keyboardFunc( self, key, x, y )
        '''
        if key == 'q' or key == 'Q': sys.exit(0)
        if key == 'c': self.captureScreen()
        if key == '\\': self.toggleFullScreen()
        '''
        
        glutPostRedisplay()


def main():
    glutInit( sys.argv )
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH )
    
    w1 = TriMeshWindow()
    w1.setWindowTitle( "Empty" )
    
    if len( sys.argv ) > 1: w1.mesh = TriMesh.FromOBJ_FileName( sys.argv[1] )
    else: print('Usage:', sys.argv[0], 'mesh.obj')
    
    #w2 = TriMeshWindow()
    #w2pos = w2.getWindowPosition()
    #w2pos = ( w2pos[0] + w1.getWindowShape()[0], w2pos[1] )
    #w2.positionWindow( w2pos )
    
    glutMainLoop()

# call main() if the script is being run
if __name__ == '__main__': main()
