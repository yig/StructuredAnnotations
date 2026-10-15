#ifndef __laplacian_editing_wrapper_h__
#define __laplacian_editing_wrapper_h__

typedef double Real;
typedef unsigned int uint;
typedef unsigned char BOOL;

/// The following was adapted from "meshopt/src/forces/LaplacianEditing.h".
/// Vectors are triplets of real_t's.

// Returns the cotangent weights and area at vertex 'center'.
// The one-ring is given by the array 'one_ring' with length 3*'valence'.
// If 'has_boundary' is true, the boundary is interpreted as being between one_ring[-1] and one_ring[0].
// Upon return, output array parameters 'cot_alpha_out' and 'cot_beta_out' will be filled with
// 'valence' values, corresponding to cotangent( alpha ) and cotangent( beta ) for
// the angles alpha and beta opposite each of the 'valence' edges ( one_ring[i] - center ).
// Upon return, output parameter 'area_out' will be set to the area owned by vertex 'center'.
// TODO Q: what about the reflected values in case of 'has_boundary'?
// TODO: return 0 for cot( alpha ) and cot( beta ) that lie beyond the boundary
void cotangentWeights(
    Real* center,
    Real* one_ring, uint valence,
    BOOL has_boundary,
    Real* cot_alpha_out, Real* cot_beta_out, Real* area_out
    );

// Returns in output parameters 'indices_out' and 'values_out' values such that
//      Hn = \sum_{ 0 <= i <= valence } mesh.vertex( indices_out[i] ) * values_out[i]
// Output parameters 'indices_out' and 'values_out' must be able to hold valence+1 values.
// Arguments are similar to those for cotangentWeights() and Hn(), except
// pass vertex indices instead of Vector3d's for 'center' and 'one_ring'.
// Also pass the output parameters from cotangentWeights(): 'cot_alpha', 'cot_beta', 'area'.
void Hn_coefficients(
    int center,
    int* one_ring, uint valence,
    BOOL has_boundary,
    Real* cot_alpha, Real* cot_beta, Real area,
    int* indices_out, Real* values_out
    );

#endif /* __laplacian_editing_wrapper_h__ */
