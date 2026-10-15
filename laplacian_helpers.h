#ifndef __laplacian_helpers_h__
#define __laplacian_helpers_h__

// #include "adt/cvec3t.h"
#include "minivec3_t.h"
#define CVec3T minivec3_t
#include <valarray>
#include <cassert>

#define AREA_MIXED 1

namespace laplacian_helpers
{
typedef unsigned int uint;

// code lifted from Shells3d.
template< typename T >
T VoronoiArea( const CVec3T<T>& v_i1, const CVec3T<T>& v_i2, const CVec3T<T>& w_12 )
{
    const CVec3T<T> x_ji = -v_i1; // v_i - v_j;
    const CVec3T<T> x_ki = -v_i2; // v_i - v_k;

    const CVec3T<T> x_jk = w_12; // v_k - v_j;
    const CVec3T<T> x_kj = -x_jk;

    const T cot_alpha_ij = dot(x_ki,x_kj) / len(cross(x_ki,x_kj));
    const T cot_beta_ik = dot(x_ji,x_jk) / len(cross(x_ji,x_jk));

    const T areaVor = 1 / 8.0 * (cot_alpha_ij * dot(x_ji,x_ji) + cot_beta_ik * dot(x_ki,x_ki));
    
    return areaVor;
}

template< typename T >
T TriangleArea( const CVec3T<T>& v_i1, const CVec3T<T>& v_i2 )
{
    return 0.5 * len( cross( v_i1, v_i2 ) );
}

template< typename T >
T AreaMixed( const CVec3T<T>& v_i1, const CVec3T<T>& v_i2, const CVec3T<T>& w_12 )
{
    // if the angle at the center vertex is obtuse
    if( dot( v_i1, v_i2 ) < 0. )
    {
        return TriangleArea( v_i1, v_i2 ) / 2.;
    }
    // if the triangle is obtuse
    else if( dot( -v_i1, w_12 ) < 0. || dot( v_i2, w_12 ) < 0. )
    {
        return TriangleArea( v_i1, v_i2 ) / 4.;
    }
    // else return voronoi area
    else
    {
        return VoronoiArea( v_i1, v_i2, w_12 );
    }
}

template< typename vec_type, typename real_type >
void
cotangentWeightsT(
    const vec_type& center,
    const vec_type* one_ring, const uint valence,
    const bool has_boundary,
    real_type* cot_alpha_out, real_type* cot_beta_out, real_type* area_out
    )
{
    assert( one_ring );
    assert( cot_alpha_out && cot_beta_out && area_out );
    
    const vec_type& p_center = center;
    const vec_type* p = one_ring;
    
/*
           ->
       p0   w0    p1
        \       /   
         v0    v1  w1\
          \   /       _|
p5 <-v5- pcenter -v2-> p2
          /   \        
         v4    v3   w2/
        /       \    L
       p4    w3   p3
             <-
*/
    // w is difference between adjacent vertices on the one-ring
    // (around the tire), pointing in increasing index direction.
    std::valarray< vec_type > w( valence );
    // v is difference between vertex on the one-ring and vertex in the center
    // (spokes of the wheel).
    std::valarray< vec_type > v( valence );
    for( uint i = 0 ; i < valence ; ++i )
    {
        w[i] = p[ (i+1) % valence ] - p[i];
        v[i] = p[i] - p_center;
    }
    
    // A is total area of one-ring / 3.
    real_type A( 0. );
    // area forward is area between vertex at center, vertex on the one-ring,
    // and vertex on the one-ring in the increasing index direction.
    std::valarray< real_type > area_forward( valence );
    for( uint i = 0 ; i < valence ; ++i )
    {
        if( has_boundary && i + 1 == valence )
        {
            area_forward[i] = 0.;
            continue;
        }
        
        area_forward[i] = 0.5 * len( cross( w[i], v[i] ) );
        
#if AREA_MIXED
        A += AreaMixed( v[i], v[ (i+1) % valence ], w[i] );
#else
        A += 1./3. * area_forward[i];
#endif
    }
    
    // cot( alpha ) is the angle opposite the v[] vector in the increasing index direction.
    real_type* cot_alpha = cot_alpha_out;
    // cot( beta ) is the angle opposite the v[] vector in the decreasing index direction.
    real_type* cot_beta = cot_beta_out;
    for( uint i = 0 ; i < valence ; ++i )
    {
        const uint i_next( (i+1) % valence );
        const uint i_prev( (i+valence-1) % valence );
    	cot_alpha[i] = dot( w[i], v[ i_next ] ) / ( 2. * area_forward[i] );
	    cot_beta[i] = dot( w[ i_prev ], -v[ i_prev ] ) / ( 2. * area_forward[ i_prev ] );
    	// std::cout << "cot(alpha): " << cot_alpha[i] << ", cot(beta): " << cot_beta[i] << std::endl ;
    	
        if( has_boundary && i_next < i ) cot_alpha[i] = 0.;
        if( has_boundary && i_prev > i ) cot_beta[i] = 0.;
    }
    
    *area_out = A;
}

template< typename vec_type, typename real_type >
void
cotangentWeights(
    const vec_type& center,
    const vec_type* one_ring, const uint valence,
    const bool has_boundary,
    real_type* cot_alpha_out, real_type* cot_beta_out, real_type* area_out
    )
{
    cotangentWeightsT( center, one_ring, valence, has_boundary, cot_alpha_out, cot_beta_out, area_out );
}

/*
// Returns the normal at vertex 'center' scaled by the mean curvature H.
// The one-ring is given by the array 'one_ring' with length 'valence'.
// If 'has_boundary' is true, the boundary is interpreted as being between one_ring[-1] and one_ring[0].
Vector3d
LaplacianEditing::Hn(
    const Vector3d& center,
    const Vector3d* one_ring, const uint valence,
    const bool has_boundary
    )
{
    Real A;
    Real cot_alpha[ valence ], cot_beta[ valence ];
    cotangentWeights(
        center, one_ring, valence, has_boundary, // input parameters
        cot_alpha, cot_beta, &A // output parameters
        );
    
    // Hn = 1/(2*A) \sum ( cot alpha_i + cot beta_i )( v[i] ),
    // where H is mean curvature and n is the normal to the center vertex.
    // HnA = Hn * A
    Vector3d HnA( 0,0,0 );
    for( uint i = 0 ; i < valence ; ++i )
    {
        HnA += (cot_alpha[i] + cot_beta[i]) * (one_ring[i] - center);
    }
    return HnA / (-2*A);
}
*/

// Returns in output parameters 'indices_out' and 'values_out' values such that
//      Hn = \sum_{ 0 <= i <= valence } mesh.vertex( indices_out[i] ) * values_out[i]
// Output parameters 'indices_out' and 'values_out' must be able to hold valence+1 values.
// Arguments are similar to those for cotangentWeights() and Hn(), except
// pass vertex indices instead of Vector3d's for 'center' and 'one_ring'.
// Also pass the output parameters from cotangentWeights(): 'cot_alpha', 'cot_beta', 'area'.
template< typename real_type >
void
Hn_coefficients(
    const int center,
    const int* one_ring, const uint valence,
    const bool has_boundary,
    const real_type* cot_alpha, const real_type* cot_beta, const real_type area,
    int* indices_out, real_type* values_out
    )
{
    // This code follows Hn().
    // NOTE: We don't use 'has_boundary', but if we change boundary handling it might be
    //       necessary in the future, and in the meantime it keeps input uniformity
    //       with cotangentWeights() and Hn().
    
    assert( valence >= 0 ); // In case someone passed a negative value (signed).
    assert( valence < 1000 ); // In case someone passed a negative value (unsigned).
    assert( one_ring );
    assert( cot_alpha && cot_beta );
    assert( area >= 0. );
    assert( center >= 0 );
    for( uint i = 0; i < valence; ++i ) assert( one_ring[i] >= 0 );
    assert( indices_out && values_out );
    
    const real_type coeff = 1./(-2*area);
    
    // Initialize the center vertex.
    indices_out[ valence ] = center;
    values_out[ valence ] = 0.;
    for( uint i = 0; i < valence; ++i )
    {
        indices_out[i] = one_ring[i];
        
        const real_type val = coeff * ( cot_alpha[i] + cot_beta[i] );
        values_out[i] = val;
        values_out[valence] -= val;
    }
    
    
    /*
    // The following test doesn't fail.  The test is in
    //  void LaplacianTiltEditingOptimizerSimple::set_Hn_constraint( const VertexID vertex_id, const Vector3d& Hn, bool incremental )
#ifndef NDEBUG
        {
            Vector3d Hn( 0,0,0 );
            for( int i = 0; i < v_indices.size(); ++i )
            {
                Hn +=
                    v_coeffs.at(i)
                    * cast( ( incremental ? mMeshCopy : undeformed_mesh )
                        .point( TriMesh::VertexHandle( v_indices.at(i) ) ) )
                    ;
            }
            std::cout << "Hn (" << (incremental?"deformed":"undeformed") << "): " << Hn << std::endl;
            assert( len( Hn - get_Hn( vertex_id, incremental ) ) < kEps );
        }
#endif
    */
}

}

#endif /* __laplacian_helpers_h__ */
