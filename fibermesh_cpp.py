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

try:
    libfibermesh = ctypes.cdll.LoadLibrary( 'libfibermesh.' + platform_shared_library_suffix() )
except OSError:
    raise ImportError('Could not load "' + 'libfibermesh.' + platform_shared_library_suffix() + '"' )
int_t = ctypes.c_int
real_t = ctypes.c_double

def to_ctypes_array( sequence, ctypes_type ):
    return ascontiguousarray( asarray( sequence, dtype = ctypes_type ) )
def ctypes_array_ref( ctypes_array, ctypes_type ):
    return ctypes_array.ctypes.data_as( ctypes.POINTER( ctypes_type ) )


class optimize_surface2_params_t( ctypes.Structure ):
    _fields_ = [
        ( "max_iterations", int_t ),
        ( "absolute_convergence_threshold", real_t ),
        ( "graph_not_cotangent_weights", int_t ),
        ( "average_edge_lengths", int_t ),
        ( "vertex_normal_motion_only", int_t )
        ]
__initialize_optimize_surface2_params = libfibermesh.initialize_optimize_surface2_params
__initialize_optimize_surface2_params.argtypes = \
    [
        ctypes.POINTER( optimize_surface2_params_t )
    ]
__initialize_optimize_surface2_params.restype = None

__optimize_surface2 = libfibermesh.optimize_surface2
__optimize_surface2.argtypes = \
    [
        ctypes.c_void_p,
        
        ## x hard constraints
        int_t,
        ctypes.POINTER( int_t ),
        ctypes.POINTER( real_t ),
        ## x linear constraints
        int_t,
        ctypes.POINTER( int_t ),
        ctypes.POINTER( int_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        
        ## y hard constraints
        int_t,
        ctypes.POINTER( int_t ),
        ctypes.POINTER( real_t ),
        ## y linear constraints
        int_t,
        ctypes.POINTER( int_t ),
        ctypes.POINTER( int_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        
        ## z hard constraints
        int_t,
        ctypes.POINTER( int_t ),
        ctypes.POINTER( real_t ),
        ## z linear constraints
        int_t,
        ctypes.POINTER( int_t ),
        ctypes.POINTER( int_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        
        ## general linear constraints
        int_t,
        ctypes.POINTER( int_t ),
        ctypes.POINTER( int_t ),
        ctypes.POINTER( int_t ),
        ctypes.POINTER( real_t ),
        ctypes.POINTER( real_t ),
        
        ctypes.POINTER( optimize_surface2_params_t )
    ]
__optimize_surface2.restype = real_t
def optimize_surface2( mesh, constraints, **kwargs ):
    '''
    Given a TriMesh 'mesh', Constraints 'constraints',
    perform 'iterations' iterations of the non-linear inflation process.
    Returns the maximum displacement of any vertex in the last iteration.
    
    NOTE: Modifies 'mesh' in place.
    '''
    
    kwargs.setdefault( "max_iterations", None )
    kwargs.setdefault( "absolute_convergence_threshold", None )
    kwargs.setdefault( "graph_not_cotangent_weights", None )
    kwargs.setdefault( "average_edge_lengths", None )
    kwargs.setdefault( "vertex_normal_motion_only", None )
    kwargs.setdefault( "name", "" )
    
    optimize_surface_proxy( mesh, constraints, __optimize_surface2, **kwargs )

__optimize_surface_linear = libfibermesh.optimize_surface_linear
__optimize_surface_linear.argtypes = __optimize_surface2.argtypes
__optimize_surface_linear.restype = __optimize_surface2.restype
def optimize_surface_linear( mesh, constraints, **kwargs ):
    '''
    Given a TriMesh 'mesh', Constraints 'constraints',
    perform 'iterations' iterations of the non-linear inflation process.
    Returns the maximum displacement of any vertex in the last iteration.
    
    NOTE: Modifies 'mesh' in place.
    '''
    
    kwargs.setdefault( "max_iterations", None )
    kwargs.setdefault( "absolute_convergence_threshold", None )
    kwargs.setdefault( "graph_not_cotangent_weights", None )
    kwargs.setdefault( "average_edge_lengths", None )
    kwargs.setdefault( "vertex_normal_motion_only", None )
    kwargs.setdefault( "name", "" )
    
    optimize_surface_proxy( mesh, constraints, __optimize_surface_linear, **kwargs )

def optimize_surface_proxy( mesh, constraints, optimize_surface_ctypes_func, **kwargs ):
    '''
    Internal proxy func used to marshall the constraints and mesh for passing to the
    ctypes-wrapped function 'optimize_surface_ctypes_func'.
    '''
    
    #print 'constraints:'
    #pprint( constraints.__dict__ )
    constraints.assert_self( mesh )
    
    xyz_hard_and_linear_constraints = list( constraints.mymath_xy_hard_and_linear_constraints() )
    xyz_hard_and_linear_constraints.append( constraints.mymath_z_hard_and_linear_constraints() )
    from fibermesh import Constraints
    general_linear_constraints = list(map( Constraints.general_linear_constraint_collect_terms, constraints.general_linear_constraints ))
    
    #print 'len( mesh.vs ):', len( mesh.vs )
    #print 'len( mesh.faces ):', len( mesh.faces )
    openmesh = new_openmesh_from_data( mesh.vs, mesh.faces )
    
    #print 'len( openmesh.vs ):', openmesh_num_vertices( openmesh )
    #mesh.vs = asarray( mesh.vs )
    #mesh.vs[:] = extract_openmesh_vertices( openmesh )
    #print 'mesh.vs.shape:', mesh.vs.shape
    #delete_openmesh( openmesh )
    #return
    
    xyz_hard_constraints_num = []
    xyz_hard_constraints_indices = []
    xyz_hard_constraints_values = []
    xyz_linear_constraints_num = []
    xyz_linear_constraints_eq_num = []
    xyz_linear_constraints_indices = []
    xyz_linear_constraints_coeffs = []
    xyz_linear_constraints_rhss = []
    for c in (0,1,2):
        '''
        // These are a pair of equal-length 1d arrays.  The num parameter is the length.
        '''
        xyz_hard_constraints_num.append( len( xyz_hard_and_linear_constraints[c][0] ) )
        xyz_hard_constraints_indices.append( [ index for index, value in xyz_hard_and_linear_constraints[c][0] ] )
        xyz_hard_constraints_values.append( [ value for index, value in xyz_hard_and_linear_constraints[c][0] ] )
        
        '''
        // These are three arrays.  The indices and coeffs parameters represent an array of sub-arrays,
        // but have been flattened end to end into two 1d arrays.
        // There are x_linear_constraints_num arrays and the length of each sub-array is stored in
        // x_linear_constraints_eq_num.
        // The rhss parameter is a 1d array with x_linear_constraints_num values.
        '''
        xyz_linear_constraints_num.append( len( xyz_hard_and_linear_constraints[c][1] ) )
        xyz_linear_constraints_eq_num.append( [] )
        xyz_linear_constraints_indices.append( [] )
        xyz_linear_constraints_coeffs.append( [] )
        xyz_linear_constraints_rhss.append( [] )
        for indices, vals, rhs in xyz_hard_and_linear_constraints[c][1]:
            xyz_linear_constraints_eq_num[-1].append( len( indices ) )
            xyz_linear_constraints_indices[-1].extend( indices )
            xyz_linear_constraints_coeffs[-1].extend( vals )
            xyz_linear_constraints_rhss[-1].append( rhs )
    
    for c in (0,1,2):
        xyz_hard_constraints_indices[c] = to_ctypes_array( xyz_hard_constraints_indices[c], int_t )
        xyz_hard_constraints_values[c] = to_ctypes_array( xyz_hard_constraints_values[c], real_t )
        xyz_linear_constraints_eq_num[c] = to_ctypes_array( xyz_linear_constraints_eq_num[c], int_t )
        xyz_linear_constraints_indices[c] = to_ctypes_array( xyz_linear_constraints_indices[c], int_t )
        xyz_linear_constraints_coeffs[c] = to_ctypes_array( xyz_linear_constraints_coeffs[c], real_t )
        xyz_linear_constraints_rhss[c] = to_ctypes_array( xyz_linear_constraints_rhss[c], real_t )
    
    '''
    // These are four arrays.  The vertex_indices and coordinate_offsets and coeffs parameters
    // represent an array of sub-arrays, but have been flattened end to end into three 1d arrays.
    // There are general_linear_constraints_num arrays and the length of each sub-array is stored in
    // general_linear_constraints_eq_num.
    // The rhss parameter is a 1d array with general_linear_constraints_num values.
    '''
    general_linear_constraints_num = len( general_linear_constraints )
    general_linear_constraints_eq_num = []
    general_linear_constraints_vertex_indices = []
    general_linear_constraints_coordinate_offsets = []
    general_linear_constraints_coeffs = []
    general_linear_constraints_rhss = []
    for indices, vals, rhs in general_linear_constraints:
        general_linear_constraints_eq_num.append( len( indices ) )
        general_linear_constraints_vertex_indices.extend( [ index for index, coord in indices ] )
        general_linear_constraints_coordinate_offsets.extend( [ coord for index, coord in indices ] )
        general_linear_constraints_coeffs.extend( vals )
        general_linear_constraints_rhss.append( rhs )
    
    general_linear_constraints_eq_num = to_ctypes_array( general_linear_constraints_eq_num, int_t )
    general_linear_constraints_vertex_indices = to_ctypes_array( general_linear_constraints_vertex_indices, int_t )
    general_linear_constraints_coordinate_offsets = to_ctypes_array( general_linear_constraints_coordinate_offsets, int_t )
    general_linear_constraints_coeffs = to_ctypes_array( general_linear_constraints_coeffs, real_t )
    general_linear_constraints_rhss = to_ctypes_array( general_linear_constraints_rhss, real_t )
    
    params = optimize_surface2_params_t()
    __initialize_optimize_surface2_params( ctypes.byref( params ) )
    if kwargs.setdefault("max_iterations", None) is not None: params.max_iterations = kwargs["max_iterations"]
    if kwargs.setdefault("absolute_convergence_threshold", None) is not None: params.absolute_convergence_threshold = kwargs["absolute_convergence_threshold"]
    if kwargs.setdefault("graph_not_cotangent_weights", None) is not None: params.graph_not_cotangent_weights = kwargs["graph_not_cotangent_weights"]
    if kwargs.setdefault("average_edge_lengths", None) is not None: params.average_edge_lengths = kwargs["average_edge_lengths"]
    if kwargs.setdefault("vertex_normal_motion_only", None) is not None: params.vertex_normal_motion_only = kwargs["vertex_normal_motion_only"]
    
    result = optimize_surface_ctypes_func(
        openmesh,
        
        
        xyz_hard_constraints_num[0],
        ctypes_array_ref( xyz_hard_constraints_indices[0], int_t ),
        ctypes_array_ref( xyz_hard_constraints_values[0], real_t ),
        xyz_linear_constraints_num[0],
        ctypes_array_ref( xyz_linear_constraints_eq_num[0], int_t ),
        ctypes_array_ref( xyz_linear_constraints_indices[0], int_t ),
        ctypes_array_ref( xyz_linear_constraints_coeffs[0], real_t ),
        ctypes_array_ref( xyz_linear_constraints_rhss[0], real_t ),
        
        xyz_hard_constraints_num[1],
        ctypes_array_ref( xyz_hard_constraints_indices[1], int_t ),
        ctypes_array_ref( xyz_hard_constraints_values[1], real_t ),
        xyz_linear_constraints_num[1],
        ctypes_array_ref( xyz_linear_constraints_eq_num[1], int_t ),
        ctypes_array_ref( xyz_linear_constraints_indices[1], int_t ),
        ctypes_array_ref( xyz_linear_constraints_coeffs[1], real_t ),
        ctypes_array_ref( xyz_linear_constraints_rhss[1], real_t ),
        
        xyz_hard_constraints_num[2],
        ctypes_array_ref( xyz_hard_constraints_indices[2], int_t ),
        ctypes_array_ref( xyz_hard_constraints_values[2], real_t ),
        xyz_linear_constraints_num[2],
        ctypes_array_ref( xyz_linear_constraints_eq_num[2], int_t ),
        ctypes_array_ref( xyz_linear_constraints_indices[2], int_t ),
        ctypes_array_ref( xyz_linear_constraints_coeffs[2], real_t ),
        ctypes_array_ref( xyz_linear_constraints_rhss[2], real_t ),
        
        
        general_linear_constraints_num,
        ctypes_array_ref( general_linear_constraints_eq_num, int_t ),
        ctypes_array_ref( general_linear_constraints_vertex_indices, int_t ),
        ctypes_array_ref( general_linear_constraints_coordinate_offsets, int_t ),
        ctypes_array_ref( general_linear_constraints_coeffs, real_t ),
        ctypes_array_ref( general_linear_constraints_rhss, real_t ),
        
        ctypes.byref( params )
        )
    
    ## Update mesh.vs in-place.  If mesh.vs is not already an array, make it one.
    mesh.vs = asarray( mesh.vs )
    mesh.vs[:] = extract_openmesh_vertices( openmesh )
    delete_openmesh( openmesh )
    del openmesh
    
    mesh.positions_changed()
    return result

__new_openmesh_from_data = libfibermesh.new_openmesh_from_data
__new_openmesh_from_data.argtypes = \
    [
        int_t,
        ctypes.POINTER( real_t ),
        int_t,
        ctypes.POINTER( int_t )
    ]
__new_openmesh_from_data.restype = ctypes.c_void_p
def new_openmesh_from_data( vertices, faces ):
    '''
    Given a list of 3d points 'vertices' and a list 'faces' of triplets of indices into 'vertices',
    constructs and returns a handle to a corresponding OpenMesh mesh.
    
    tested
    '''
    
    n_vertices = len( vertices )
    n_faces = len( faces )
    vertices = ascontiguousarray( asarray( vertices, dtype = real_t ) )
    faces = ascontiguousarray( asarray( faces, dtype = int_t ) )
    return __new_openmesh_from_data(
        n_vertices,
        vertices.ctypes.data_as( ctypes.POINTER( real_t ) ),
        n_faces,
        faces.ctypes.data_as( ctypes.POINTER( int_t ) )
        )


__delete_openmesh = libfibermesh.delete_openmesh
__delete_openmesh.argtypes = \
    [
        ctypes.c_void_p
    ]
__delete_openmesh.restype = None
def delete_openmesh( openmesh ):
    '''
    Given a handle to an OpenMesh mesh (as returned by new_openmesh_from_data),
    deletes the mesh.
    After this call, the OpenMesh mesh handle is invalid.
    
    used
    '''
    
    __delete_openmesh( openmesh )


__update_openmesh_vertices = libfibermesh.update_openmesh_vertices
__update_openmesh_vertices.argtypes = \
    [
        ctypes.c_void_p,
        ctypes.POINTER( real_t )
    ]
__update_openmesh_vertices.restype = None
def update_openmesh_vertices( openmesh, vertices ):
    '''
    Given a list of 3d points 'vertices', one for each vertex in the OpenMesh mesh 'openmesh',
    replaces the positions of vertices in the mesh with those in vertices.
    
    tested
    '''
    
    vertices = ascontiguousarray( asarray( vertices, dtype = real_t ) )
    __update_openmesh_vertices(
        openmesh,
        vertices.ctypes.data_as( ctypes.POINTER( real_t ) )
        )


__extract_openmesh_vertices = libfibermesh.extract_openmesh_vertices
__extract_openmesh_vertices.argtypes = \
    [
        ctypes.c_void_p,
        ctypes.POINTER( real_t )
    ]
__extract_openmesh_vertices.restype = None
def extract_openmesh_vertices( openmesh ):
    '''
    Returns a list of 3d points 'vertices', one for each vertex in the OpenMesh mesh 'openmesh'.
    
    tested
    '''
    
    N = openmesh_num_vertices( openmesh )
    vertices = ascontiguousarray( empty( (N, 3), dtype = real_t ) )
    __extract_openmesh_vertices(
        openmesh,
        vertices.ctypes.data_as( ctypes.POINTER( real_t ) )
        )
    return vertices

__extract_openmesh_vertex_normals = libfibermesh.extract_openmesh_vertex_normals
__extract_openmesh_vertex_normals.argtypes = \
    [
        ctypes.c_void_p,
        ctypes.POINTER( real_t )
    ]
__extract_openmesh_vertex_normals.restype = None
def extract_openmesh_vertex_normals( openmesh ):
    '''
    Returns an array of 3d vectors, one for each vertex normal in the OpenMesh mesh 'openmesh'.
    
    tested
    '''
    
    N = openmesh_num_vertices( openmesh )
    vertex_normals = ascontiguousarray( empty( (N, 3), dtype = real_t ) )
    __extract_openmesh_vertex_normals(
        openmesh,
        vertex_normals.ctypes.data_as( ctypes.POINTER( real_t ) )
        )
    return vertex_normals

__extract_openmesh_face_normals = libfibermesh.extract_openmesh_face_normals
__extract_openmesh_face_normals.argtypes = \
    [
        ctypes.c_void_p,
        ctypes.POINTER( real_t )
    ]
__extract_openmesh_face_normals.restype = None
def extract_openmesh_face_normals( openmesh ):
    '''
    Returns an array of 3d vectors, one for each face normal in the OpenMesh mesh 'openmesh'.
    
    tested
    '''
    
    N = openmesh_num_faces( openmesh )
    face_normals = ascontiguousarray( empty( (N, 3), dtype = real_t ) )
    __extract_openmesh_face_normals(
        openmesh,
        face_normals.ctypes.data_as( ctypes.POINTER( real_t ) )
        )
    return face_normals

__extract_openmesh_face_areas = libfibermesh.extract_openmesh_face_areas
__extract_openmesh_face_areas.argtypes = \
    [
        ctypes.c_void_p,
        ctypes.POINTER( real_t )
    ]
__extract_openmesh_face_areas.restype = None
def extract_openmesh_face_areas( openmesh ):
    '''
    Returns an array of real_t's, one for each face area in the OpenMesh mesh 'openmesh'.
    
    tested
    '''
    
    N = openmesh_num_faces( openmesh )
    face_areas = ascontiguousarray( empty( N, dtype = real_t ) )
    __extract_openmesh_face_areas(
        openmesh,
        face_areas.ctypes.data_as( ctypes.POINTER( real_t ) )
        )
    return face_areas

__openmesh_num_vertices = libfibermesh.openmesh_num_vertices
__openmesh_num_vertices.argtypes = \
    [
        ctypes.c_void_p
    ]
__openmesh_num_vertices.restype = int_t
def openmesh_num_vertices( openmesh ):
    '''
    Returns the number of vertices in the OpenMesh mesh 'openmesh'.
    
    tested
    '''
    
    return __openmesh_num_vertices( openmesh )

__openmesh_num_faces = libfibermesh.openmesh_num_faces
__openmesh_num_faces.argtypes = \
    [
        ctypes.c_void_p
    ]
__openmesh_num_faces.restype = int_t
def openmesh_num_faces( openmesh ):
    '''
    Returns the number of faces in the OpenMesh mesh 'openmesh'.
    
    tested
    '''
    
    return __openmesh_num_faces( openmesh )
