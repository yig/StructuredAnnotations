#if 0
c++ laplacian_editing_wrapper.cpp \
    -Iext \
    -dynamiclib -o liblaplacian_editing.dylib \
    -g -O1 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare
#endif


extern "C"
{
#include "laplacian_editing_wrapper.h"
}

#include "laplacian_helpers.h"

// #include "adt/cvec3t.h"
// typedef CVec3T<Real> Vector3d;
#include "minivec3_t.h"
typedef minivec3_t<Real> Vector3d;

extern "C"
{
#include "laplacian_editing_wrapper.h"

void cotangentWeights(
    Real* center,
    Real* one_ring, uint valence,
    BOOL has_boundary,
    Real* cot_alpha_out, Real* cot_beta_out, Real* area_out
    )
{
    // Convert the input parameters into Vector3ds.
    Vector3d centerr( center[0], center[1], center[2] );
    Vector3d* one_ringg = new Vector3d[ valence ];
    for( uint i = 0; i < valence; ++i )
    {
        one_ringg[i](0) = one_ring[ 3*i + 0 ];
        one_ringg[i](1) = one_ring[ 3*i + 1 ];
        one_ringg[i](2) = one_ring[ 3*i + 2 ];
    }
    
    laplacian_helpers::cotangentWeights(
        centerr, one_ringg, valence, has_boundary,
        cot_alpha_out, cot_beta_out, area_out
        );
    delete[] one_ringg;
}

void Hn_coefficients(
    int center,
    int* one_ring, uint valence,
    BOOL has_boundary,
    Real* cot_alpha, Real* cot_beta, Real area,
    int* indices_out, Real* values_out
    )
{
    laplacian_helpers::Hn_coefficients(
        center, one_ring, valence, has_boundary,
        cot_alpha, cot_beta, area,
        indices_out, values_out
        );
}

}
