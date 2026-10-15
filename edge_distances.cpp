/*
## The best way to turn this into a dynamic library is to run a custom setup.py.

## On MacOS X:
c++ edge_distances.cpp -dynamiclib -o edge_distances.dylib -g -O3 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare
*/


#include <cassert>
#include <iostream>

// #include "adt/cvec2t.h"
// #include "adt/cvec3t.h"
#include "minivec2_t.h"
#include "minivec3_t.h"
#define CVec2T minivec2_t
#define CVec3T minivec3_t

namespace edge_distances
{
double L_infinity = 3;

// linestrip is an array of vec_t's of length 'npts'.
// gradients_out must be an array of length 'npts-1' (the number of edges).
template< typename real_t, typename vec_t >
void point_distances_and_gradients_to_linestrip_edges( const vec_t& point, const vec_t* linestrip, int npts, real_t* distances_out, vec_t* gradients_out )
{
    assert( npts >= 2 );
    
    vec_t p_b = point - linestrip[0];
    for( int i = 0; i < npts-1; ++i )
    {
        vec_t p_a = p_b;
        p_b = point - linestrip[i+1];
        
        vec_t b_a = linestrip[i+1] - linestrip[i];
        vec_t d = b_a / lenSq( b_a );
        real_t cond = dot( p_a, d );
        
        if( cond < 0 )
        {
            distances_out[i] = len( p_a );
            gradients_out[i] = p_a;
        }
        else if( cond > 1 )
        {
            distances_out[i] = len( p_b );
            gradients_out[i] = p_b;
        }
        else
        {
            vec_t L = p_a - cond * b_a;
            distances_out[i] = len( L );
            gradients_out[i] = L - dot( L, d ) * b_a;
        }
        
        gradients_out[i] /= distances_out[i];
    }
}


// linestrip is an array of vec_t's of length 'npts'.
template< typename real_t, typename vec_t >
vec_t smooth_gradient_to_linestrip_edges( const vec_t& point, const vec_t* linestrip, int npts )
{
    assert( npts >= 2 );
    
    real_t* distances = new real_t[ npts-1 ];
    vec_t* gradients = new vec_t[ npts-1 ];
    
    point_distances_and_gradients_to_linestrip_edges( point, linestrip, npts, distances, gradients );
    
    real_t sum_distances = 0.;
    for( int i = 0; i < npts-1; ++i )
    {
        sum_distances += pow( distances[i], -L_infinity );
    }
    
    // I hope all vec_t's are 3-dimensional or less.
    real_t zeros[] = { 0., 0., 0. };
    vec_t smooth_gradient( zeros );
    real_t coeff = (-1./L_infinity) * pow( sum_distances, -1./L_infinity - 1 );
    for( int i = 0; i < npts-1; ++i )
    {
        smooth_gradient += coeff * ( ( -L_infinity * pow( distances[i], -L_infinity - 1 ) ) * gradients[i] );
    }
    
    delete [] gradients;
    delete [] distances;
    
    return smooth_gradient;
}

template< typename real_t, typename vec_t >
vec_t trace_point( const vec_t& point, const vec_t* linestrip, int npts, real_t distance, real_t dt )
{
    assert( npts >= 2 );
    assert( distance >= 0. );
    assert( dt > 0 );
    
    vec_t trace = point;
    // std::cout << "point: " << point << '\n';
    // int steps = 1;
    
    const bool kNormalizeOnce = false;
    real_t norm = -1.;
    while( distance > 0. )
    {
        const vec_t smooth_gradient = smooth_gradient_to_linestrip_edges< real_t, vec_t >( trace, linestrip, npts );
        
        if( norm == -1. || !kNormalizeOnce )
        {
            norm = 1. / len( smooth_gradient );
        }
        
        /*
        std::cout << "distance: " << distance << '\n';
        std::cout << "norm: " << norm << '\n';
        std::cout << "smooth_gradient: " << smooth_gradient << '\n';
        std::cout << "step: " << ( ( std::max( std::min( distance, dt ), real_t(0) ) * norm ) * smooth_gradient ) << '\n';
        */
        
        trace += ( std::max( std::min( distance, dt ), real_t(0) ) * norm ) * smooth_gradient;
        distance -= dt;
        
        // steps += 1;
        
        // std::cout << "trace: " << trace << '\n';
    }
    
    // std::cout << "steps: " << steps << '\n';
    
    return trace;
}
}

extern "C" {

typedef double real_t;

void point_distances_and_gradients_to_linestrip_edges(
    int dim, real_t* point, real_t* linestrip, int npts, real_t* distances_out, real_t* gradients_out
    )
{
    assert( point && linestrip && distances_out && gradients_out );
    assert( dim > 0 && npts > 2 );
    
    if( dim == 2 )
    {
        edge_distances::point_distances_and_gradients_to_linestrip_edges< real_t, CVec2T<real_t> >(
            *reinterpret_cast< CVec2T<real_t>* >( point ),
            reinterpret_cast< CVec2T<real_t>* >( linestrip ),
            npts,
            distances_out,
            reinterpret_cast< CVec2T<real_t>* >( gradients_out )
            );
    }
    else if( dim == 3 )
    {
        edge_distances::point_distances_and_gradients_to_linestrip_edges< real_t, CVec3T<real_t> >(
            *reinterpret_cast< CVec3T<real_t>* >( point ),
            reinterpret_cast< CVec3T<real_t>* >( linestrip ),
            npts,
            distances_out,
            reinterpret_cast< CVec3T<real_t>* >( gradients_out )
            );
    }
    else
    {
        assert( dim == 2 || dim == 3 );
    }
}

void trace_point( int dim, real_t* point, real_t* linestrip, int npts, real_t distance, real_t dt, real_t* trace_out )
{
    assert( point && linestrip && trace_out );
    assert( dim > 0 && npts > 2 );
    
    if( dim == 2 )
    {
        CVec2T<real_t> result = edge_distances::trace_point< real_t, CVec2T<real_t> >(
            *reinterpret_cast< CVec2T<real_t>* >( point ),
            reinterpret_cast< CVec2T<real_t>* >( linestrip ),
            npts,
            distance,
            dt
            );
        
        trace_out[0] = result(0);
        trace_out[1] = result(1);
    }
    else if( dim == 3 )
    {
        CVec3T<real_t> result = edge_distances::trace_point< real_t, CVec3T<real_t> >(
            *reinterpret_cast< CVec3T<real_t>* >( point ),
            reinterpret_cast< CVec3T<real_t>* >( linestrip ),
            npts,
            distance,
            dt
            );
        
        trace_out[0] = result(0);
        trace_out[1] = result(1);
        trace_out[2] = result(2);
    }
    else
    {
        assert( dim == 2 || dim == 3 );
    }
}
}
