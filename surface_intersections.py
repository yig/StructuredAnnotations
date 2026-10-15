from math import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *

try:
    import ctypes
except:
    print('''
ERROR: ctypes not installed properly. (ctypes needed for laplacian_editing_1d())
        ''')
    import sys
    sys.exit()

def platform_shared_library_suffix():
    import sys
    result = 'so'
    if 'win' in sys.platform.lower(): result = 'dll'
    ## No else if, because we want darwin to override win (which is a substring of darwin)
    if 'darwin' in sys.platform.lower(): result = 'dylib'
    return result

libtritri = ctypes.cdll.LoadLibrary( 'libtritri.' + platform_shared_library_suffix() )
real_t = ctypes.c_float
#real_t = ctypes.c_double

def to_ctypes_array( sequence, ctypes_type ):
    '''
    Returns a contiguous ctypes array suitable for extracting
    a pointer (with .ctypes.data_as() or ctypes_array_ref()) to pass to a C function.
    '''
    return ascontiguousarray( asarray( sequence, dtype = ctypes_type ) )
def ctypes_array_ref( ctypes_array, ctypes_type ):
    '''
    Returns a ctypes.POINTER to the internal data of 'ctypes_array' of type 'ctypes_type'.
    
    WARNING: Do not chain the calls to to_ctypes_array() and ctypes_array_ref() as in
               ctypes_array_ref( to_ctypes_array( sequence, ctypes_type ), ctypes_type ).
             ctypes_array_ref returns a pointer to the internal numpy storage of the array and
             cannot keep a reference to the array, so it will be garbage collected if these
             functions are chained together.
             Instead, keep a temporary reference to the result of to_ctypes_array():
               contig = to_ctypes_array( sequence, ctypes_type )
               ctypes_array_ref( contig, ctypes_type )
    '''
    return ctypes_array.ctypes.data_as( ctypes.POINTER( ctypes_type ) )

__tri_tri_intersect = libtritri.NoDivTriTriIsect
__tri_tri_intersect.argtypes = \
    [
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t )
    ]
__tri_tri_intersect.restype = ctypes.c_int
def _tri_tri_intersect( v0, v1, v2, u0, u1, u2 ):
    '''
    parameters: vertices of triangle 1: V0,V1,V2
                vertices of triangle 2: U0,U1,U2
    result    : returns 1 if the triangles intersect, otherwise 0
    '''
    
    ## I'm assuming that parameters are already stored in a ctypes-friendly format, as in:
    ##   v0 = ascontiguousarray( asarray( v0, dtype = real_t ) )
    return __tri_tri_intersect(
        v0.ctypes.data_as( ctypes.POINTER( real_t ) ),
        v1.ctypes.data_as( ctypes.POINTER( real_t ) ),
        v2.ctypes.data_as( ctypes.POINTER( real_t ) ),
        u0.ctypes.data_as( ctypes.POINTER( real_t ) ),
        u1.ctypes.data_as( ctypes.POINTER( real_t ) ),
        u2.ctypes.data_as( ctypes.POINTER( real_t ) )
        ).value


__tri_tri_intersect_with_isectline = libtritri.tri_tri_intersect_with_isectline
__tri_tri_intersect_with_isectline.argtypes = \
    [
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( ctypes.c_int ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t )
    ]
__tri_tri_intersect_with_isectline.restype = ctypes.c_int
def _tri_tri_intersect_with_isectline( v0, v1, v2, u0, u1, u2 ):
    '''
    parameters: vertices of triangle 1: V0,V1,V2
                vertices of triangle 2: U0,U1,U2
    result    : a 3-element tuple:
                    1 if the triangles intersect, otherwise 0
                    whether the tris are coplanar
                    a 2-element tuple containing the two endpoints of the line of intersection
    '''
    
    ## Convert parameters into a ctypes-friendly format.
    v0 = ascontiguousarray( asarray( v0, dtype = real_t ) )
    v1 = ascontiguousarray( asarray( v1, dtype = real_t ) )
    v2 = ascontiguousarray( asarray( v2, dtype = real_t ) )
    u0 = ascontiguousarray( asarray( u0, dtype = real_t ) )
    u1 = ascontiguousarray( asarray( u1, dtype = real_t ) )
    u2 = ascontiguousarray( asarray( u2, dtype = real_t ) )
    
    ## Make the output parameters in a ctypes-friendly way.
    isectpt1 = ascontiguousarray( zeros( 3, dtype = real_t ) )
    isectpt2 = ascontiguousarray( zeros( 3, dtype = real_t ) )
    coplanar = ctypes.c_int()
    
    intersect = __tri_tri_intersect_with_isectline(
        v0.ctypes.data_as( ctypes.POINTER( real_t ) ),
        v1.ctypes.data_as( ctypes.POINTER( real_t ) ),
        v2.ctypes.data_as( ctypes.POINTER( real_t ) ),
        u0.ctypes.data_as( ctypes.POINTER( real_t ) ),
        u1.ctypes.data_as( ctypes.POINTER( real_t ) ),
        u2.ctypes.data_as( ctypes.POINTER( real_t ) ),
        ctypes.byref( coplanar ),
        isectpt1.ctypes.data_as( ctypes.POINTER( real_t ) ),
        isectpt2.ctypes.data_as( ctypes.POINTER( real_t ) )
        )
    
    return intersect.value, coplanar.value, ( isectpt1, isectpt2 )

__tritri_wrapper_intersect_ellipsoids = libtritri.tritri_wrapper_intersect_ellipsoids
__tritri_wrapper_intersect_ellipsoids.argtypes = \
    [
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( ctypes.c_int )
    ]
__tritri_wrapper_intersect_ellipsoids.restype = ctypes.POINTER( real_t )
__tritri_wrapper_free = libtritri.tritri_wrapper_free
__tritri_wrapper_free.argtypes = [ ctypes.POINTER( real_t ) ]
__tritri_wrapper_free.restype = None
def _tritri_wrapper_intersect_ellipsoids( center1, xform1, center2, xform2 ):
    '''
    Returns an array containing a flattened sequence of line segments (pairs of points) on the intersection
    curve between two ellipsoids.
    The result has the layout
    [ [ [ line0.start.x, line0.start.y, line0.start.z, ], [ line0.end.x, line0.end.y, line0.end.z, ] ], [ [ line1.start.x, ... ].
    Ellipsoids 1 and 2 are each specified by an array of 3 floats (xyz) representing the center
    'center1' and 'center2', and an array of 9 floats representing a row-major transform matrix
    (scaling and rotation, actually the three axes as columns) 'xform1' and 'xform2'.
    '''
    
    center1 = ascontiguousarray( asarray( center1, dtype = real_t ) )
    center2 = ascontiguousarray( asarray( center2, dtype = real_t ) )
    xform1 = ascontiguousarray( asarray( xform1, dtype = real_t ) )
    xform2 = ascontiguousarray( asarray( xform2, dtype = real_t ) )
    
    num_lines_out = ctypes.c_int()
    flattened_lines = __tritri_wrapper_intersect_ellipsoids(
        center1.ctypes.data_as( ctypes.POINTER( real_t ) ),
        xform1.ctypes.data_as( ctypes.POINTER( real_t ) ),
        center2.ctypes.data_as( ctypes.POINTER( real_t ) ),
        xform2.ctypes.data_as( ctypes.POINTER( real_t ) ),
        ctypes.byref( num_lines_out )
        )
    num_lines_out = num_lines_out.value
    result = emptyf( 6*num_lines_out )
    result[:] = flattened_lines[:6*num_lines_out]
    __tritri_wrapper_free( flattened_lines )
    return result.reshape( (num_lines_out,2,3) )

def ellipsoid_ellipsoid_intersection_lines( ellipsoid1, ellipsoid2 ):
    # return ellipsoid_ellipsoid_intersection_lines_slow( ellipsoid1, ellipsoid2 )
    
    xform1 = asarray( ellipsoid1.axes.T )
    xform2 = asarray( ellipsoid2.axes.T )
    center1 = asarray( ellipsoid1.center ).ravel()
    center2 = asarray( ellipsoid2.center ).ravel()
    
    return _tritri_wrapper_intersect_ellipsoids( center1, xform1, center2, xform2 )


_sphere_faces = None
_sphere_vertices = None
def _get_sphere():
    global _sphere_faces, _sphere_vertices
    if _sphere_faces is None:
        assert _sphere_vertices is None
        
        ## Load OBJ of the sphere or tesselate it or something.
        from mesh_gui import TriMesh
        import os
        #sphere_fname = 'sphere.obj'
        sphere_fname = 'sphere2.out.obj'
        mesh = TriMesh.FromOBJ_FileName( os.path.join( '.', 'ext', sphere_fname ) )
        _sphere_faces = mesh.faces
        _sphere_vertices = mesh.vs
    
    assert _sphere_faces is not None
    assert _sphere_vertices is not None
    return _sphere_faces, _sphere_vertices

def ellipsoid_ellipsoid_intersection_lines_slow( ellipsoid1, ellipsoid2 ):
    '''
    Given two Ellipsoid's, returns a list of pairs of 3d points.  Each pair of points represents
    a line segment.
    TODO: Collapse line segments that share a point into a line loop.
    
    tested
    '''
    
    xform1 = asarray( ellipsoid1.axes )
    xform2 = asarray( ellipsoid2.axes )
    center1 = asarray( ellipsoid1.center ).ravel()
    center2 = asarray( ellipsoid2.center ).ravel()
    
    fs, vs = _get_sphere()
    
    result = []
    
    '''
    ## This printed ~1e-11
    #diff_minus_1_squared = 0.
    #for v in vs: diff_minus_1_squared += ( mag(v) - 1. )**2
    #print 'diff_minus_1_squared:', diff_minus_1_squared
    
    for i in xrange( len(fs) ):
        #result.append( ( vs[ fs[i][0] ], vs[ fs[i][1] ] ) )
        #result.append( ( vs[ fs[i][1] ], vs[ fs[i][2] ] ) )
        #result.append( ( vs[ fs[i][0] ], vs[ fs[i][2] ] ) )
        
        #result.append( ( dot( xform1, vs[ fs[i][0] ] ) + center1, dot( xform1, vs[ fs[i][1] ] ) + center1 ) )
        #result.append( ( dot( xform1, vs[ fs[i][1] ] ) + center1, dot( xform1, vs[ fs[i][2] ] ) + center1 ) )
        #result.append( ( dot( xform1, vs[ fs[i][0] ] ) + center1, dot( xform1, vs[ fs[i][2] ] ) + center1 ) )
        
        result.append( ( dot( xform2, vs[ fs[i][0] ] ) + center2, dot( xform2, vs[ fs[i][1] ] ) + center2 ) )
        result.append( ( dot( xform2, vs[ fs[i][1] ] ) + center2, dot( xform2, vs[ fs[i][2] ] ) + center2 ) )
        result.append( ( dot( xform2, vs[ fs[i][0] ] ) + center2, dot( xform2, vs[ fs[i][2] ] ) + center2 ) )
    return result
    '''
    
    for i in range( len(fs) ):
        for j in range( len(fs) ):
            intersect, coplanar, line_segment = _tri_tri_intersect_with_isectline(
                dot( xform1, vs[ fs[i][0] ] ) + center1,
                dot( xform1, vs[ fs[i][1] ] ) + center1,
                dot( xform1, vs[ fs[i][2] ] ) + center1,
                dot( xform2, vs[ fs[j][0] ] ) + center2,
                dot( xform2, vs[ fs[j][1] ] ) + center2,
                dot( xform2, vs[ fs[j][2] ] ) + center2
                )
            if intersect and not coplanar:
                result.append( line_segment )
    
    return result

__tritri_wrapper_intersect_meshes = libtritri.tritri_wrapper_intersect_meshes
__tritri_wrapper_intersect_meshes.argtypes = \
    [
        ## The first mesh
        ctypes.POINTER( real_t ),
        ctypes.c_int,
        ctypes.POINTER( ctypes.c_int ),
        ctypes.c_int,
        
        ## The second mesh
        ctypes.POINTER( real_t ),
        ctypes.c_int,
        ctypes.POINTER( ctypes.c_int ),
        ctypes.c_int,
        
        ## An output parameter.
        ctypes.POINTER( ctypes.c_int )
    ]
__tritri_wrapper_intersect_meshes.restype = ctypes.POINTER( real_t )
def _tritri_wrapper_intersect_meshes( vs0, faces0, vs1, faces1 ):
    '''
    Returns an array containing a flattened sequence of line segments (pairs of points) on the intersection
    curve between two triangle meshes, each represented by a sequence of 3d vertices vs0,1
    and a sequence of faces, tuples of three ints indexing the vertices in vs0,1, faces0,1.
    The result has the layout
    [ [ [ line0.start.x, line0.start.y, line0.start.z, ], [ line0.end.x, line0.end.y, line0.end.z, ] ], [ [ line1.start.x, ... ].
    '''
    
    vs0 = ascontiguousarray( asarray( vs0, dtype = real_t ) )
    vs1 = ascontiguousarray( asarray( vs1, dtype = real_t ) )
    assert vs0.shape[1] == 3
    assert vs1.shape[1] == 3
    
    faces0 = ascontiguousarray( asarray( faces0, dtype = ctypes.c_int ) )
    faces1 = ascontiguousarray( asarray( faces1, dtype = ctypes.c_int ) )
    assert faces0.shape[1] == 3
    assert faces1.shape[1] == 3
    
    num_lines_out = ctypes.c_int()
    flattened_lines = __tritri_wrapper_intersect_meshes(
        vs0.ctypes.data_as( ctypes.POINTER( real_t ) ),
        len( vs0 ),
        faces0.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) ),
        len( faces0 ),
        
        vs1.ctypes.data_as( ctypes.POINTER( real_t ) ),
        len( vs1 ),
        faces1.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) ),
        len( faces1 ),
        
        ctypes.byref( num_lines_out )
        )
    num_lines_out = num_lines_out.value
    result = emptyf( 6*num_lines_out )
    result[:] = flattened_lines[:6*num_lines_out]
    __tritri_wrapper_free( flattened_lines )
    return result.reshape( (num_lines_out,2,3) )

def mesh_mesh_intersection_lines( mesh0, mesh1 ):
    '''
    Given two triangle meshes, mesh0 and mesh1, each with the fields .vs and .faces,
    where .vs is a sequence of 3D vertices and .faces is a sequence of tuples of three indices
    into .vs,
    returns the same result as _tritri_wrapper_intersect_meshes().
    '''
    
    return _tritri_wrapper_intersect_meshes( mesh0.vs, mesh0.faces, mesh1.vs, mesh1.faces )

class face_face_intersection_t( ctypes.Structure ):
    _fields_ = [
        ("x0", real_t),
        ("y0", real_t),
        ("z0", real_t),
        ("x1", real_t),
        ("y1", real_t),
        ("z1", real_t),
        
        ("face0", ctypes.c_int),
        ("face1", ctypes.c_int)
        ]
__tritri_wrapper_intersect_meshes_info = libtritri.tritri_wrapper_intersect_meshes_info
__tritri_wrapper_intersect_meshes_info.argtypes = __tritri_wrapper_intersect_meshes.argtypes
__tritri_wrapper_intersect_meshes_info.restype = ctypes.POINTER( face_face_intersection_t )
__tritri_wrapper_free_info = libtritri.tritri_wrapper_free_info
__tritri_wrapper_free_info.argtypes = [ ctypes.POINTER( face_face_intersection_t ) ]
__tritri_wrapper_free_info.restype = None
def _tritri_wrapper_intersect_meshes_info( vs0, faces0, vs1, faces1 ):
    '''
    Given two triangle meshes, each represented by a sequence of 3d vertices vs0,1
    and a sequence of faces, tuples of three ints indexing the vertices in vs0,1, faces0,1,
    returns two arrays, the first an array of line segments (pairs of 3d
    points) on the intersection curve between two triangle meshes with the layout
    [ [ [ line0.start.x, line0.start.y, line0.start.z, ], [ line0.end.x, line0.end.y, line0.end.z, ] ], [ [ line1.start.x, ... ]
    and the second an array of pairs of indices for which two faces (in faces0 and faces1)
    correspond to the intersection line.
    '''
    
    vs0 = ascontiguousarray( asarray( vs0, dtype = real_t ) )
    vs1 = ascontiguousarray( asarray( vs1, dtype = real_t ) )
    assert vs0.shape[1] == 3
    assert vs1.shape[1] == 3
    
    faces0 = ascontiguousarray( asarray( faces0, dtype = ctypes.c_int ) )
    faces1 = ascontiguousarray( asarray( faces1, dtype = ctypes.c_int ) )
    assert faces0.shape[1] == 3
    assert faces1.shape[1] == 3
    
    num_lines_out = ctypes.c_int()
    intersection_info = __tritri_wrapper_intersect_meshes_info(
        vs0.ctypes.data_as( ctypes.POINTER( real_t ) ),
        len( vs0 ),
        faces0.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) ),
        len( faces0 ),
        
        vs1.ctypes.data_as( ctypes.POINTER( real_t ) ),
        len( vs1 ),
        faces1.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) ),
        len( faces1 ),
        
        ctypes.byref( num_lines_out )
        )
    num_lines_out = num_lines_out.value
    
    result_lines = emptyf( ( num_lines_out, 2, 3 ) )
    result_faces = empty( ( num_lines_out, 2 ), dtype = int )
    
    for i in range( num_lines_out ):
        face_face_intersection = intersection_info[i]
        result_lines[i][0] = face_face_intersection.x0, face_face_intersection.y0, face_face_intersection.z0
        result_lines[i][1] = face_face_intersection.x1, face_face_intersection.y1, face_face_intersection.z1
        result_faces[i] = face_face_intersection.face0, face_face_intersection.face1
    
    __tritri_wrapper_free_info( intersection_info )
    return result_lines, result_faces

def mesh_mesh_intersection_lines_info( mesh0, mesh1 ):
    '''
    Given two triangle meshes, mesh0 and mesh1, each with the fields .vs and .faces,
    where .vs is a sequence of 3D vertices and .faces is a sequence of tuples of three indices
    into .vs,
    returns the same result as _tritri_wrapper_intersect_meshes().
    '''
    
    return _tritri_wrapper_intersect_meshes_info( mesh0.vs, mesh0.faces, mesh1.vs, mesh1.faces )

__intersect_triangle1 = libtritri.intersect_triangle1
__intersect_triangle1.argtypes = \
    [
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double )
    ]
__intersect_triangle1.restype = ctypes.c_int
def ray_triangle_intersection( ray, tri ):
    '''
    Given a ray (3d point, 3d direction) and 'triangle' (a list of three 3d points),
    returns a vector (t u v), where t is the distance to the plane in which the triangle lies
    and (u, v) represents the coordinates inside the triangle.
    If the ray does not intersect the triangle, returns None.
    NOTE: t may be negative, indicating that the intersection is in the negative ray direction.
    
    used
    '''
    
    ## Convert parameters into a ctypes-friendly format.
    rorig = ascontiguousarray( asarray( ray[0], dtype = ctypes.c_double ) )
    assert len( rorig ) == 3
    rdir = ascontiguousarray( asarray( ray[1], dtype = ctypes.c_double ) )
    assert len( rdir ) == 3
    v0 = ascontiguousarray( asarray( tri[0], dtype = ctypes.c_double ) )
    v1 = ascontiguousarray( asarray( tri[1], dtype = ctypes.c_double ) )
    v2 = ascontiguousarray( asarray( tri[2], dtype = ctypes.c_double ) )
    assert len( v0 ) == 3
    assert len( v1 ) == 3
    assert len( v2 ) == 3
    
    ## Make the output parameters in a ctypes-friendly way.
    t = ctypes.c_double(0)
    u = ctypes.c_double(0)
    v = ctypes.c_double(0)
    
    intersect = __intersect_triangle1(
        rorig.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        rdir.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        v0.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        v1.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        v2.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        ctypes.byref( t ),
        ctypes.byref( u ),
        ctypes.byref( v )
        )
    
    ## TODO Q: Why is the return value getting converted to an 'int' instead of remaining a c_int?
    # if intersect.value:
    if intersect:
        return t.value, u.value, v.value
    else:
        return None

def ray_mesh_intersections_slow( ray, mesh ):
    '''
    Given a ray (3d point, 3d direction) and 'mesh' (a structure containing a list of
    3d vertices 'vs' and a list of triangles (triplets of indices into 'vs') named 'faces'),
    returns a list of vectors ( t, fi, (u,v) ), where
        t is the distance to the plane in which the triangle lies,
        fi is the index of the face in mesh.faces, and
        (u, v) represents the coordinates inside the triangle.
    If the ray does not intersect the mesh, returns an empty list.
    
    used
    '''
    
    result = []
    for fi, face in enumerate( mesh.faces ):
        isect = ray_triangle_intersection(
            ray,
            [ mesh.vs[ vi ] for vi in face ]
            )
        if isect is None: continue
        result.append( ( isect[0], fi, ( isect[1], isect[2] ) ) )
    
    return result

class raymesh_intersection_t( ctypes.Structure ):
    _fields_ = [
        ( "t", ctypes.c_double ),
        ( "fi", ctypes.c_int ),
        ( "u", ctypes.c_double ),
        ( "v", ctypes.c_double )
        ]
__raytri_wrapper_ray_mesh_intersections = libtritri.ray_mesh_intersections
__raytri_wrapper_ray_mesh_intersections.argtypes = \
    [
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double ),
        ctypes.c_int,
        ctypes.POINTER( ctypes.c_double ),
        ctypes.c_int,
        ctypes.POINTER( ctypes.c_int ),
        ctypes.POINTER( ctypes.c_int )
    ]
__raytri_wrapper_ray_mesh_intersections.restype = ctypes.POINTER( raymesh_intersection_t )
__raytri_wrapper_ray_mesh_intersections_free = libtritri.ray_mesh_intersections_free
__raytri_wrapper_ray_mesh_intersections_free.argtypes = [ ctypes.POINTER( raymesh_intersection_t ) ]
__raytri_wrapper_ray_mesh_intersections_free.restype = None
def _raytri_wrapper_ray_mesh_intersections( ray, mesh ):
    '''
    Given a ray (3d point, 3d direction) and 'mesh' (a structure containing a list of
    3d vertices 'vs' and a list of triangles (triplets of indices into 'vs') named 'faces'),
    returns a list of vectors ( t, fi, (u,v) ), where
        t is the distance to the plane in which the triangle lies,
        fi is the index of the face in mesh.faces, and
        (u, v) represents the coordinates inside the triangle.
    If the ray does not intersect the mesh, returns an empty list.
    
    used
    '''
    
    ray = (
        to_ctypes_array( ray[0], ctypes.c_double ),
        to_ctypes_array( ray[1], ctypes.c_double )
        )
    vs = to_ctypes_array( mesh.vs, ctypes.c_double )
    faces = to_ctypes_array( mesh.faces, ctypes.c_int )
    
    num_intersections_out = ctypes.c_int()
    intersections = __raytri_wrapper_ray_mesh_intersections(
        ctypes_array_ref( ray[0], ctypes.c_double ),
        ctypes_array_ref( ray[1], ctypes.c_double ),
        len( vs ),
        ctypes_array_ref( vs, ctypes.c_double ),
        len( faces ),
        ctypes_array_ref( faces, ctypes.c_int ),
        ctypes.byref( num_intersections_out )
        )
    num_intersections_out = num_intersections_out.value
    result = [
        ( isect.t, isect.fi, ( isect.u, isect.v ) )
        for isect in intersections[:num_intersections_out]
        ]
    __raytri_wrapper_ray_mesh_intersections_free( intersections )
    
    ## These generate the same thing:
    #print 'ray_mesh_intersections fast', result
    #print 'ray_mesh_intersections slow', ray_mesh_intersections_slow( ray, mesh )
    #print 'ray_mesh_intersections difference:', set( result ).symmetric_difference( set( ray_mesh_intersections_slow( ray, mesh ) ) )
    
    return result

ray_mesh_intersections = _raytri_wrapper_ray_mesh_intersections


def mesh_silhouette_edges( mesh ):
    '''
    Given a TriMesh 'mesh' (a structure containing a list of
    3d vertices 'vs' and a list of triangles (triplets of indices into 'vs') named 'faces'
    and a list of 3d face normals 'face_normals'),
    returns a list of tuples (i,j), where (i,j) is an edge between vertices i and j,
    and the edge is a silhouette edge of the mesh when viewed along the z axis.
    
    untested
    '''
    
    result = []
    for edge in mesh.edges:
        hei0 = mesh.directed_edge2he_index( edge )
        if hei0 == -1:
            ## Boundary edges count as silhouette edges.
            result.append( edge )
            continue
        he0 = mesh.halfedges[ hei0 ]
        
        hei1 = he0.opposite_he
        if hei1 == -1:
            ## Boundary edges count as silhouette edges.
            result.append( edge )
            continue
        he1 = mesh.halfedges[ hei1 ]
        
        f0 = he0.face
        f1 = he1.face
        
        if mesh.face_normals[ f0 ][2] * mesh.face_normals[ f1 ][2] <= 0.:
            result.append( edge )
    
    return result
