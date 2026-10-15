/*
## The best way to turn this into a dynamic library is to run a custom setup.py.

## On MacOS X:
g++ Eline_ad.cpp -dynamiclib -o Eline_ad.dylib -arch i386 -arch x86_64 -g -O1 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare -Iext -Iext/autodiff -Iext/meta

## On Windows:
/cygdrive/c/MinGW/bin/migw32-g++ -c Eline_ad.cpp -Iext -DNDEBUG
/cygdrive/c/MinGW/bin/dllwrap --export-all-symbols --mno-cygwin --driver-name=c++ --dllname Eline_ad.dll Eline_ad.o
## Or is it?
g++ -c Eline_ad.cpp -Iext -DNDEBUG -shared -o Eline_ad.dll
*/

#include <cassert>
#include <iostream>

//#include "adt/adreal.h"
//#include "adt/advec.h"

#include "minivec3_t.h"
#include <autodiff.h>

typedef double real_t;

namespace
{
    bool kDebugging = false;
    
    template< typename T, typename U >
    void set_independent( T& advec, const U& vec3, int offset, int total )
    {
        assert( offset + 3 <= total );
        
        for( int i = 0; i < 3; ++i )
        {
            advec(i) = vec3(i);
            advec(i).diff( offset + i, total );
        }
    }
}

template< int DERIVS >
sl::Fad< real_t, DERIVS >
Emag_and_derivatives(
    real_t* p_0_in, real_t* p_1_in
    )
{
    typedef sl::Fad< real_t, DERIVS > adreal_t;
    typedef minivec3_t< real_t > vec_t;
    typedef minivec3_t< adreal_t > advec_t;
    
    advec_t p_0;
    advec_t p_1;
    set_independent( p_0, vec_t( p_0_in ), 0, 6 );
    set_independent( p_1, vec_t( p_1_in ), 3, 6 );
    
    if( kDebugging )
    {
        std::cout << "p_0: " << p_0 << '\n';
        std::cout << "p_1: " << p_1 << '\n';
    }
    
    // The energy will be the length of one point minus the other, squared.
    
    // Find the vector from p0 to p1.
    advec_t p01 = p_1 - p_0;
    
    // Set the energy to the magnitude, squared of the vector.
    adreal_t E = dot( p01, p01 );
    
    return E;
}


template< int DERIVS >
sl::Fad< real_t, DERIVS >
Ecurvature_and_derivatives(
    real_t* p_0_in, real_t* p_1_in, real_t* p_2_in
    )
{
    typedef sl::Fad< real_t, DERIVS > adreal_t;
    typedef minivec3_t< real_t > vec_t;
    typedef minivec3_t< adreal_t > advec_t;
    
    advec_t p_0;
    advec_t p_1;
    advec_t p_2;
    set_independent( p_0, vec_t( p_0_in ), 0, 9 );
    set_independent( p_1, vec_t( p_1_in ), 3, 9 );
    set_independent( p_2, vec_t( p_2_in ), 6, 9 );
    
    if( kDebugging )
    {
        std::cout << "p_0: " << p_0 << '\n';
        std::cout << "p_1: " << p_1 << '\n';
        std::cout << "p_2: " << p_2 << '\n';
    }
    
    // Q: Does setting lengths to 1 help this converge?
    // A: No.
    //adreal_t p_10mag = 1.;
    //adreal_t p_12mag = 1.;
    adreal_t p_10mag = len( p_0 - p_1 );
    adreal_t p_12mag = len( p_2 - p_1 );
    if( kDebugging )
    {
        std::cout << "p_10mag: " << p_10mag << '\n';
        std::cout << "p_12mag: " << p_12mag << '\n';
    }
    
    adreal_t A = .5 * ( p_10mag + p_12mag );
    if( kDebugging )
    {
        std::cout << "A: " << A << '\n';
    }
    
    adreal_t beta = dot( p_0 - p_1, p_2 - p_1 ) / ( p_10mag * p_12mag );
    if( kDebugging )
    {
        std::cout << "beta: " << beta << '\n';
    }
    
    /*
    // This expression has a sin of an acos(), which we can get around with some trig identities.
    adreal_t alpha = acos( beta );
    std::cout << "alpha: " << alpha << '\n';
    
    // Q: Why divide by A to get the curvature?
    // A: It looks like we do this to convert the integrated quantity into a point-wise quantity.
    //    Note that the energy is the curvature^2 times the area; multiplying by the area
    //    takes into account the area of integration (h).
    adreal_t curvature = 2. * sin( ( M_PI - alpha ) / 2. ) / A;
    std::cout << "curvature: " << curvature << '\n';
    */
    
    /*
    // Using the half angle formula sin( a / 2 ) = +- sqrt( ( 1 - cos(a) ) / 2 ):
    //   curvature = 2 * (+-) sqrt( ( 1 - cos( M_PI - alpha ) ) / 2 ) / A
    //             = 2 * (+-) sqrt( ( 1 + cos( alpha ) ) / 2 ) / A
    //             = 2 * (+-) sqrt( ( 1 + beta ) / 2 ) / A
    // Since E = curvature^2 * A, we don't care about the +-.  We also don't care about the sqrt().
    // std::cout << "sqrt( 1. + beta ): " << sqrt( 1. + beta ) << '\n';
    adreal_t curvature = 2. * sqrt( ( 1. + beta ) / 2. ) / A;
    std::cout << "curvature: " << curvature << '\n';
    // If curvature is 0, the sqrt must be of 0, and the gradient has 1 / sqrt() in it,
    // so the gradient is nan.
    adreal_t E = A*curvature*curvature;
    std::cout << "E: " << E << '\n';
    */
    
    // Using the half angle formula, and the fact that we're squaring the curvature,
    //   E = A * curvature * curvature
    //     = A * ( 2 * (+-) sqrt( ( 1 + beta ) / 2 ) / A ) * ( 2 * (+-) sqrt( ( 1 + beta ) / 2 ) / A )
    //     = 4 * ( ( 1 + beta ) / 2 ) / A
    //     = 2 * ( 1 + beta ) / A
    adreal_t E = 2. * ( 1. + beta ) / A;
    
    // If it's curvature = 2 * tan( ( pi - alpha ) / 2 ) / A, as in [Bergou et al 2008] then it's
    //   E = A * curvature * curvature
    //     = A * ( 2 * (+-) sqrt( ( 1 + beta ) / ( 1 - beta ) ) / A ) * ( 2 * (+-) sqrt( ( 1 + beta ) / ( 1 - beta ) ) / A )
    //     = 4 * ( ( 1 + beta ) / ( 1 - beta ) ) / A
    // adreal_t E = 4. * ( ( 1. + beta ) / ( 1. - beta ) ) / A;
    
    return E;
}


template< int DERIVS >
sl::Fad< real_t, DERIVS >
Ecurvature_and_derivatives_undeformed(
    real_t* p_0_in, real_t* p_1_in, real_t* p_2_in,
    real_t* u_0_in, real_t* u_1_in, real_t* u_2_in
    )
{
    typedef sl::Fad< real_t, DERIVS > adreal_t;
    typedef minivec3_t< real_t > vec_t;
    typedef minivec3_t< adreal_t > advec_t;
    
    advec_t p_0;
    advec_t p_1;
    advec_t p_2;
    set_independent( p_0, vec_t( p_0_in ), 0, 9 );
    set_independent( p_1, vec_t( p_1_in ), 3, 9 );
    set_independent( p_2, vec_t( p_2_in ), 6, 9 );
    
    vec_t u_0( u_0_in );
    vec_t u_1( u_1_in );
    vec_t u_2( u_2_in );
    
    if( kDebugging )
    {
        std::cout << "p_0: " << p_0 << '\n';
        std::cout << "p_1: " << p_1 << '\n';
        std::cout << "p_2: " << p_2 << '\n';
        
        std::cout << "u_0: " << u_0 << '\n';
        std::cout << "u_1: " << u_1 << '\n';
        std::cout << "u_2: " << u_2 << '\n';
    }
    
    real_t u_10mag = len( u_0 - u_1 );
    real_t u_12mag = len( u_2 - u_1 );
    if( kDebugging )
    {
        std::cout << "u_10mag: " << u_10mag << '\n';
        std::cout << "u_12mag: " << u_12mag << '\n';
    }
    
    real_t A = .5 * ( u_10mag + u_12mag );
    if( kDebugging )
    {
        std::cout << "A: " << A << '\n';
    }
    
    adreal_t beta = dot( p_0 - p_1, p_2 - p_1 ) / ( u_10mag * u_12mag );
    if( kDebugging )
    {
        std::cout << "beta: " << beta << '\n';
    }
    
    /*
    // This expression has a sin of an acos(), which we can get around with some trig identities.
    adreal_t alpha = acos( beta );
    std::cout << "alpha: " << alpha << '\n';
    
    // Q: Why divide by A to get the curvature?
    // A: It looks like we do this to convert the integrated quantity into a point-wise quantity.
    //    Note that the energy is the curvature^2 times the area; multiplying by the area
    //    takes into account the area of integration (h).
    adreal_t curvature = 2. * sin( ( M_PI - alpha ) / 2. ) / A;
    std::cout << "curvature: " << curvature << '\n';
    */
    
    /*
    // Using the half angle formula sin( a / 2 ) = +- sqrt( ( 1 - cos(a) ) / 2 ):
    //   curvature = 2 * (+-) sqrt( ( 1 - cos( M_PI - alpha ) ) / 2 ) / A
    //             = 2 * (+-) sqrt( ( 1 + cos( alpha ) ) / 2 ) / A
    //             = 2 * (+-) sqrt( ( 1 + beta ) / 2 ) / A
    // Since E = curvature^2 * A, we don't care about the +-.  We also don't care about the sqrt().
    // std::cout << "sqrt( 1. + beta ): " << sqrt( 1. + beta ) << '\n';
    adreal_t curvature = 2. * sqrt( ( 1. + beta ) / 2. ) / A;
    std::cout << "curvature: " << curvature << '\n';
    // If curvature is 0, the sqrt must be of 0, and the gradient has 1 / sqrt() in it,
    // so the gradient is nan.
    adreal_t E = A*curvature*curvature;
    std::cout << "E: " << E << '\n';
    */
    
    // Using the half angle formula, and the fact that we're squaring the curvature,
    //   E = A * curvature * curvature
    //     = A * ( 2 * (+-) sqrt( ( 1 + beta ) / 2 ) / A ) * ( 2 * (+-) sqrt( ( 1 + beta ) / 2 ) / A )
    //     = 4 * ( ( 1 + beta ) / 2 ) / A
    //     = 2 * ( 1 + beta ) / A
    adreal_t E = 2. * ( 1. + beta ) / A;
    
    // If it's curvature = 2 * tan( ( pi - alpha ) / 2 ) / A, as in [Bergou et al 2008] then it's
    //   E = A * curvature * curvature
    //     = A * ( 2 * (+-) sqrt( ( 1 + beta ) / ( 1 - beta ) ) / A ) * ( 2 * (+-) sqrt( ( 1 + beta ) / ( 1 - beta ) ) / A )
    //     = 4 * ( ( 1 + beta ) / ( 1 - beta ) ) / A
    // adreal_t E = 4. * ( ( 1. + beta ) / ( 1. - beta ) ) / A;
    
    return E;
}


template< int DERIVS >
sl::Fad< real_t, DERIVS >
Elaplacian_editing_and_derivatives(
    real_t* p_0_in, real_t* p_1_in, real_t* p_2_in
    )
{
    typedef sl::Fad< real_t, DERIVS > adreal_t;
    typedef minivec3_t< real_t > vec_t;
    typedef minivec3_t< adreal_t > advec_t;
    
    advec_t p_0;
    advec_t p_1;
    advec_t p_2;
    set_independent( p_0, vec_t( p_0_in ), 0, 9 );
    set_independent( p_1, vec_t( p_1_in ), 3, 9 );
    set_independent( p_2, vec_t( p_2_in ), 6, 9 );
    
    if( kDebugging )
    {
        std::cout << "p_0: " << p_0 << '\n';
        std::cout << "p_1: " << p_1 << '\n';
        std::cout << "p_2: " << p_2 << '\n';
    }
    
    adreal_t p_10mag = len( p_0 - p_1 );
    adreal_t p_12mag = len( p_2 - p_1 );
    if( kDebugging )
    {
        std::cout << "p_10mag: " << p_10mag << '\n';
        std::cout << "p_12mag: " << p_12mag << '\n';
    }
    
    // Q: Does setting lengths to 1 help this converge?
    // A: Yes.
    p_10mag = 1.;
    p_12mag = 1.;
    
    adreal_t a = p_12mag;
    adreal_t b = p_10mag;
    adreal_t Area = a+b;
    
    // This comes from Nealen's PhD thesis.
    advec_t LArea = Area*p_1 - a*p_0 - b*p_2;
    
    adreal_t E = dot( LArea, LArea ) / Area;
    
    /*
    advec_t grad12 = (p_2 - p_1) / p_12mag;
    advec_t grad01 = (p_1 - p_0) / p_10mag;
    advec_t L = (grad12 - grad01) / (.5*Area);
    adreal_t E = dot( L,L ) * .5 * Area;
    */
    
    return E;
}


template< int DERIVS >
sl::Fad< real_t, DERIVS >
Elaplacian_editing_and_derivatives_undeformed(
    real_t* p_0_in, real_t* p_1_in, real_t* p_2_in,
    real_t* u_0_in, real_t* u_1_in, real_t* u_2_in
    )
{
    typedef sl::Fad< real_t, DERIVS > adreal_t;
    typedef minivec3_t< real_t > vec_t;
    typedef minivec3_t< adreal_t > advec_t;
    
    advec_t p_0;
    advec_t p_1;
    advec_t p_2;
    if( p_0_in ) set_independent( p_0, vec_t( p_0_in ), 0, 9 );
    set_independent( p_1, vec_t( p_1_in ), 3, 9 );
    if( p_2_in ) set_independent( p_2, vec_t( p_2_in ), 6, 9 );
    
    vec_t u_0;
    if( u_0_in ) u_0 = vec_t( u_0_in );
    vec_t u_1( u_1_in );
    vec_t u_2;
    if( u_2_in ) u_2 = vec_t( u_2_in );
    
    if( kDebugging )
    {
        if( p_0_in ) std::cout << "p_0: " << p_0 << '\n';
        std::cout << "p_1: " << p_1 << '\n';
        if( p_2_in ) std::cout << "p_2: " << p_2 << '\n';
        
        if( u_0_in ) std::cout << "u_0: " << u_0 << '\n';
        std::cout << "u_1: " << u_1 << '\n';
        if( u_2_in ) std::cout << "u_2: " << u_2 << '\n';
    }
    
    real_t u_10mag( 0. );
    real_t u_12mag( 0. );
    if( u_0_in ) u_10mag = len( u_0 - u_1 );
    if( u_2_in ) u_12mag = len( u_2 - u_1 );
    if( kDebugging )
    {
        if( u_0_in ) std::cout << "u_10mag: " << u_10mag << '\n';
        if( u_1_in ) std::cout << "u_12mag: " << u_12mag << '\n';
    }
    
    adreal_t E;
    if( !p_0_in )
    {
        assert( !u_0_in );
        
        adreal_t a = u_12mag;
        adreal_t Area = a;
        
        // This comes from Nealen's PhD thesis.
        // TODO Q: Is this the right thing to do at endpoints?
        //         It's not quite the generalization, because that would have b*p_2, which is zero.
        advec_t LArea = Area*p_1 - a*p_2;
        
        E = dot( LArea, LArea ) / Area;
    }
    else if( !p_2_in )
    {
        assert( !u_2_in );
        
        adreal_t b = u_10mag;
        adreal_t Area = b;
        
        // This comes from Nealen's PhD thesis.
        // TODO Q: Is this the right thing to do at endpoints?
        //         It's not quite the generalization, because that would have a*p_0, which is zero.
        advec_t LArea = Area*p_1 - b*p_0;
        
        E = dot( LArea, LArea ) / Area;
    }
    else
    {
        adreal_t a = u_12mag;
        adreal_t b = u_10mag;
        adreal_t Area = a+b;
        
        // This comes from Nealen's PhD thesis.
        advec_t LArea = Area*p_1 - a*p_0 - b*p_2;
        
        E = dot( LArea, LArea ) / Area;
    }
    
    /*
    advec_t grad12 = (p_2 - p_1) / p_12mag;
    advec_t grad01 = (p_1 - p_0) / p_10mag;
    advec_t L = (grad12 - grad01) / (.5*Area);
    adreal_t E = dot( L,L ) * .5 * Area;
    */
    
    return E;
}


extern "C" {

void Emag(
    real_t* p_0, real_t* p_1,
    real_t* E_out,
    real_t* dE_dp_out,
    real_t* d2E_dp2_out
    )
{
    assert( p_0 );
    assert( p_1 );
    assert( E_out );
    assert( dE_dp_out );
    
    static const int N = 6;
    
    if( d2E_dp2_out )
    {
        sl::Fad< real_t, 2 > e = Emag_and_derivatives<2>( p_0, p_1 );
        *E_out = e.getV();
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = e.getD(i);
        }
        for( int i = 0; i < N; ++i )
        for( int j = 0; j < N; ++j )
        {
            d2E_dp2_out[ N*i + j ] = e.getH(i,j);
        }
    }
    else
    {
        sl::Fad< real_t, 1 > e = Emag_and_derivatives<1>( p_0, p_1 );
        *E_out = e.getV();
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = e.getD(i);
        }
    }
}

void Ecurvature(
    real_t* p_0, real_t* p_1, real_t* p_2,
    real_t* E_out,
    real_t* dE_dp_out,
    real_t* d2E_dp2_out
    )
{
    assert( p_0 );
    assert( p_1 );
    assert( p_2 );
    assert( E_out );
    assert( dE_dp_out );
    
    static const int N = 9;
    
    if( d2E_dp2_out )
    {
        sl::Fad< real_t, 2 > e = Ecurvature_and_derivatives<2>( p_0, p_1, p_2 );
        *E_out = e.getV();
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = e.getD(i);
        }
        for( int i = 0; i < N; ++i )
        for( int j = 0; j < N; ++j )
        {
            d2E_dp2_out[ N*i + j ] = e.getH(i,j);
        }
    }
    else
    {
        sl::Fad< real_t, 1 > e = Ecurvature_and_derivatives<1>( p_0, p_1, p_2 );
        *E_out = e.getV();
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = e.getD(i);
        }
    }
}

void Ecurvature_undeformed(
    real_t* p_0, real_t* p_1, real_t* p_2,
    real_t* u_0, real_t* u_1, real_t* u_2,
    real_t* E_out,
    real_t* dE_dp_out,
    real_t* d2E_dp2_out
    )
{
    assert( p_0 );
    assert( p_1 );
    assert( p_2 );
    assert( E_out );
    assert( dE_dp_out );
    
    static const int N = 9;
    
    if( d2E_dp2_out )
    {
        sl::Fad< real_t, 2 > e = Ecurvature_and_derivatives_undeformed<2>( p_0, p_1, p_2, u_0, u_1, u_2 );
        *E_out = e.getV();
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = e.getD(i);
        }
        for( int i = 0; i < N; ++i )
        for( int j = 0; j < N; ++j )
        {
            d2E_dp2_out[ N*i + j ] = e.getH(i,j);
        }
    }
    else
    {
        sl::Fad< real_t, 1 > e = Ecurvature_and_derivatives_undeformed<1>( p_0, p_1, p_2, u_0, u_1, u_2 );
        *E_out = e.getV();
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = e.getD(i);
        }
    }
}

void Elaplacian_editing(
    real_t* p_0, real_t* p_1, real_t* p_2,
    real_t* E_out,
    real_t* dE_dp_out,
    real_t* d2E_dp2_out
    )
{
    assert( p_0 );
    assert( p_1 );
    assert( p_2 );
    assert( E_out );
    assert( dE_dp_out );
    
    static const int N = 9;
    
    if( d2E_dp2_out )
    {
        sl::Fad< real_t, 2 > e = Elaplacian_editing_and_derivatives<2>( p_0, p_1, p_2 );
        *E_out = e.getV();
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = e.getD(i);
        }
        for( int i = 0; i < N; ++i )
        for( int j = 0; j < N; ++j )
        {
            d2E_dp2_out[ N*i + j ] = e.getH(i,j);
        }
    }
    else
    {
        sl::Fad< real_t, 1 > e = Elaplacian_editing_and_derivatives<1>( p_0, p_1, p_2 );
        *E_out = e.getV();
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = e.getD(i);
        }
    }
}

void Elaplacian_editing_undeformed(
    real_t* p_0, real_t* p_1, real_t* p_2,
    real_t* u_0, real_t* u_1, real_t* u_2,
    real_t* E_out,
    real_t* dE_dp_out,
    real_t* d2E_dp2_out
    )
{
    assert( (p_0 == NULL) == (u_0 == NULL) );
    assert( (p_2 == NULL) == (u_2 == NULL) );
    
    // u/p_0 and u/p_2 can be NULL, which means we are at an endpoint of the line.
    // assert( p_0 );
    if( kDebugging ) std::cout << "p_0: " << p_0 << '\n';
    assert( p_1 );
    // assert( p_2 );
    if( kDebugging ) std::cout << "p_2: " << p_2 << '\n';
    
    // assert( u_0 );
    if( kDebugging ) std::cout << "u_0: " << u_0 << '\n';
    assert( u_1 );
    // assert( u_2 );
    if( kDebugging ) std::cout << "u_2: " << u_2 << '\n';
    
    assert( E_out );
    assert( dE_dp_out );
    
    static const int N = 9;
    
    if( d2E_dp2_out )
    {
        sl::Fad< real_t, 2 > e = Elaplacian_editing_and_derivatives_undeformed<2>( p_0, p_1, p_2, u_0, u_1, u_2 );
        *E_out = e.getV();
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = e.getD(i);
        }
        for( int i = 0; i < N; ++i )
        for( int j = 0; j < N; ++j )
        {
            d2E_dp2_out[ N*i + j ] = e.getH(i,j);
        }
    }
    else
    {
        sl::Fad< real_t, 1 > e = Elaplacian_editing_and_derivatives_undeformed<1>( p_0, p_1, p_2, u_0, u_1, u_2 );
        *E_out = e.getV();
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = e.getD(i);
        }
    }
}

void Etest(
    real_t* p_0, real_t* p_1, real_t* p_2,
    real_t* E_out,
    real_t* dE_dp_out,
    real_t* d2E_dp2_out
    )
{
    assert( p_0 );
    assert( p_1 );
    assert( p_2 );
    assert( E_out );
    assert( dE_dp_out );
    
    static const int N = 9;
    
    if( d2E_dp2_out )
    {
        *E_out = -31337;
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = i;
        }
        for( int i = 0; i < N; ++i )
        for( int j = 0; j < N; ++j )
        {
            d2E_dp2_out[ N*i + j ] = N*i + j;
        }
    }
    else
    {
        *E_out = -31337;
        for( int i = 0; i < N; ++i )
        {
            dE_dp_out[i] = i;
        }
    }
}

}
