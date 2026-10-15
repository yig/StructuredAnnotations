#ifndef __tritri_wrapper_h__
#define __tritri_wrapper_h__

typedef float real_t;

// Returns a pointer to a flattened sequence of line segments (pairs of points) on the intersection
// curve between two ellipsoids.
// The result has the layout
// [ line0.start.x, line0.start.y, line0.start.z, line0.end.x, line0.end.y, line0.end.z, line1.start.x, ... ]
// The number of lines is returned in the integer output parameter 'num_lines_out'.
// Ellipsoids 1 and 2 are each specified by an array of 3 floats (xyz) representing the center
// 'center1' and 'center2', and an array of 9 floats representing a row-major transform matrix
// (scaling and rotation, actually the three axes as columns) 'xform1' and 'xform2'.
real_t* tritri_wrapper_intersect_ellipsoids( real_t* center1, real_t* xform1, real_t* center2, real_t* xform2, int* num_lines_out );
// Returns a pointer to a flattened sequence of line segments (pairs of points) on the intersection
// curve between two meshes.
// The result has the layout
// [ line0.start.x, line0.start.y, line0.start.z, line0.end.x, line0.end.y, line0.end.z, line1.start.x, ... ]
// The number of lines is returned in the integer output parameter 'num_lines_out'.
// The two meshes are each specified by an array of vertices vs (triplets of xyz floats)
// and an array of triangular faces (triplets of indices into vs).
real_t* tritri_wrapper_intersect_meshes(
    real_t* vs0, int num_vertices0, int* faces0, int num_faces0,
    real_t* vs1, int num_vertices1, int* faces1, int num_faces1,
    int* num_lines_out );
// Call this to free the memory returned by tritri_wrapper_intersect_*()
void tritri_wrapper_free( real_t* memory );

// Same as tritri_wrapper_intersect_meshes(), but returns a pointer to an array of
// face_face_intersection_t structs.
typedef struct {
    // The xyz positions of the start and end of the line segment.
    real_t x0, y0, z0, x1, y1, z1;
    // The index of the intersected face in faces0 and faces1.
    int face0, face1;
} face_face_intersection_t;
face_face_intersection_t* tritri_wrapper_intersect_meshes_info(
    real_t* vs0, int num_vertices0, int* faces0, int num_faces0,
    real_t* vs1, int num_vertices1, int* faces1, int num_faces1,
    int* num_lines_out );
void tritri_wrapper_free_info( face_face_intersection_t* memory );

#endif /* __tritri_wrapper_h__ */
