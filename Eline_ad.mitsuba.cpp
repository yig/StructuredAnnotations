/*
## The best way to turn this into a dynamic library is to run a custom setup.py.

## On MacOS X:
c++ Eline_ad.mitsuba.cpp -dynamiclib -o Eline_ad.dylib -g -O3 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare -Iext -Iext/mitsuba -I/usr/local/include/eigen3

## On Windows:
/cygdrive/c/MinGW/bin/migw32-g++ -c Eline_ad.cpp -Iext -DNDEBUG
/cygdrive/c/MinGW/bin/dllwrap --export-all-symbols --mno-cygwin --driver-name=c++ --dllname Eline_ad.dll Eline_ad.o
## Or is it?
g++ -c Eline_ad.cpp -Iext -DNDEBUG -shared -o Eline_ad.dll
*/

#include <cassert>
#include <iostream>

#include "minivec3_t.h"
#include "autodiff.h"

DECLARE_DIFFSCALAR_BASE();

typedef double real_t;
typedef minivec3_t< real_t > vec_t;
typedef DScalar1< real_t, Eigen::Matrix< real_t, 6, 1 > > adrealG6_t;
typedef DScalar1< real_t, Eigen::Matrix< real_t, 9, 1 > > adrealG9_t;
typedef DScalar2< real_t, Eigen::Matrix< real_t, 6, 1 >, Eigen::Matrix< real_t, 6, 6 > > adrealH6_t;
typedef DScalar2< real_t, Eigen::Matrix< real_t, 9, 1 >, Eigen::Matrix< real_t, 9, 9 > > adrealH9_t;

namespace
{
    bool kDebugging = false;
    
    template< typename adreal_t >
    void fill_energy_gradient( const adreal_t& e, real_t* E_out, real_t* dE_dp_out )
    {
        assert( E_out );
        assert( dE_dp_out );
        
        // Value
        *E_out = e.getValue();
        // Gradient
        Eigen::Map< typename adreal_t::Gradient >( dE_dp_out, DiffScalarBase::getVariableCount(), 1 ) = e.getGradient();
    }
    
    template< typename adreal_t >
    void fill_energy_gradient_hessian(
        const adreal_t& e,
        real_t* E_out,
        real_t* dE_dp_out,
        real_t* d2E_dp2_out
        )
    {
        fill_energy_gradient( e, E_out, dE_dp_out );
        
        // Hessian
        assert( d2E_dp2_out );
        Eigen::Map< typename adreal_t::Hessian >( d2E_dp2_out, DiffScalarBase::getVariableCount(), DiffScalarBase::getVariableCount() ) = e.getHessian();
    }
    
    template< typename adreal_t >
    void fill_vec3( minivec3_t< adreal_t >& dst, const real_t* src, const int dof_off )
    {
        assert( src );
        for( int i = 0; i < 3; ++i ) dst(i) = adreal_t( dof_off + i, src[i] );
    }
    
    // It is annoying to have to do this, but I get an error:
    //     implicit instantiation of undefined template 'Eigen::internal::traits ...'
    // if I don't force temporaries for complex expressions.
    template< typename scalar_t, typename advec_t >
    advec_t
    compute_LArea( const scalar_t& s1, const advec_t& v1, const scalar_t& s2, const advec_t& v2 )
    {
        advec_t tmp1 = v1;
        tmp1 *= s1;
        
        advec_t tmp2 = v2;
        tmp2 *= s2;
        
        advec_t result = tmp1;
        result += tmp2;
        
        return result;
    }
    
    template< typename scalar_t, typename advec_t >
    advec_t
    compute_LArea( const scalar_t& s1, const advec_t& v1, const scalar_t& s2, const advec_t& v2, const scalar_t& s3, const advec_t& v3 )
    {
        advec_t tmp1 = v1;
        tmp1 *= s1;
        
        advec_t tmp2 = v2;
        tmp2 *= s2;
        
        advec_t tmp3 = v3;
        tmp3 *= s3;
        
        advec_t result = tmp1;
        result += tmp2;
        result += tmp3;
        
        return result;
    }
}

template< typename adreal_t >
adreal_t
Emag_and_derivatives(
    const real_t* p_0_in, const real_t* p_1_in
    )
{
    typedef minivec3_t< adreal_t > advec_t;
    advec_t p_0, p_1;
    fill_vec3( p_0, p_0_in, 0 );
    fill_vec3( p_1, p_1_in, 3 );
    
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


template< typename adreal_t >
adreal_t
Ecurvature_and_derivatives(
    const real_t* p_0_in, const real_t* p_1_in, const real_t* p_2_in
    )
{
    typedef minivec3_t< adreal_t > advec_t;
    advec_t p_0, p_1, p_2;
    fill_vec3( p_0, p_0_in, 0 );
    fill_vec3( p_1, p_1_in, 3 );
    fill_vec3( p_2, p_2_in, 6 );
    
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


template< typename adreal_t >
adreal_t
Ecurvature_and_derivatives_undeformed(
    const real_t* p_0_in, const real_t* p_1_in, const real_t* p_2_in,
    const vec_t& u_0, const vec_t& u_1, const vec_t& u_2
    )
{
    typedef minivec3_t< adreal_t > advec_t;
    advec_t p_0, p_1, p_2;
    fill_vec3( p_0, p_0_in, 0 );
    fill_vec3( p_1, p_1_in, 3 );
    fill_vec3( p_2, p_2_in, 6 );
    
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


template< typename adreal_t >
adreal_t
Elaplacian_editing_and_derivatives(
    const real_t* p_0_in, const real_t* p_1_in, const real_t* p_2_in
    )
{
    typedef minivec3_t< adreal_t > advec_t;
    advec_t p_0, p_1, p_2;
    fill_vec3( p_0, p_0_in, 0 );
    fill_vec3( p_1, p_1_in, 3 );
    fill_vec3( p_2, p_2_in, 6 );
    
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
    
    const adreal_t a = p_12mag;
    const adreal_t b = p_10mag;
    const adreal_t Area = a+b;
    
    // This comes from Nealen's PhD thesis.
    // const advec_t LArea = Area*p_1 - a*p_0 - b*p_2;
    const advec_t LArea = compute_LArea( Area,p_1, -a,p_0, -b,p_2 );
    
    adreal_t E = dot( LArea, LArea ) / Area;
    
    /*
    advec_t grad12 = (p_2 - p_1) / p_12mag;
    advec_t grad01 = (p_1 - p_0) / p_10mag;
    advec_t L = (grad12 - grad01) / (.5*Area);
    adreal_t E = dot( L,L ) * .5 * Area;
    */
    
    return E;
}


template< typename adreal_t >
adreal_t
Elaplacian_editing_and_derivatives_undeformed(
    const real_t* p_0_in, const real_t* p_1_in, const real_t* p_2_in,
    const real_t* u_0_in, const real_t* u_1_in, const real_t* u_2_in
    )
{
    typedef minivec3_t< adreal_t > advec_t;
    advec_t p_0, p_1, p_2;
    if( p_0_in ) fill_vec3( p_0, p_0_in, 0 );
    fill_vec3( p_1, p_1_in, 3 );
    if( p_2_in ) fill_vec3( p_2, p_2_in, 6 );
    
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
        
        const real_t a = u_12mag;
        const real_t Area = a;
        
        // This comes from Nealen's PhD thesis.
        // TODO Q: Is this the right thing to do at endpoints?
        //         It's not quite the generalization, because that would have b*p_2, which is zero.
        // const advec_t LArea = Area*p_1 - a*p_2;
        const advec_t LArea = compute_LArea( Area,p_1, -a,p_2 );
        
        E = dot( LArea, LArea ) / Area;
    }
    else if( !p_2_in )
    {
        assert( !u_2_in );
        
        const real_t b = u_10mag;
        const real_t Area = b;
        
        // This comes from Nealen's PhD thesis.
        // TODO Q: Is this the right thing to do at endpoints?
        //         It's not quite the generalization, because that would have a*p_0, which is zero.
        // const advec_t LArea = Area*p_1 - b*p_0;
        const advec_t LArea = compute_LArea( Area,p_1, -b,p_0 );
        
        E = dot( LArea, LArea ) / Area;
    }
    else
    {
        const real_t a = u_12mag;
        const real_t b = u_10mag;
        const real_t Area = a+b;
        
        // This comes from Nealen's PhD thesis.
        // const advec_t LArea = Area*p_1 - a*p_0 - b*p_2;
        const advec_t LArea = compute_LArea( Area,p_1, -a,p_0, -b,p_2 );
        
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
    
    DiffScalarBase::setVariableCount(6);
    
    if( d2E_dp2_out )
    {
        adrealH6_t e = Emag_and_derivatives< adrealH6_t >( p_0, p_1 );
        fill_energy_gradient_hessian( e, E_out, dE_dp_out, d2E_dp2_out );
    }
    else
    {
        adrealG6_t e = Emag_and_derivatives< adrealG6_t >( p_0, p_1 );
        fill_energy_gradient( e, E_out, dE_dp_out );
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
    
    DiffScalarBase::setVariableCount(9);
    
    if( d2E_dp2_out )
    {
        adrealH9_t e = Ecurvature_and_derivatives< adrealH9_t >( p_0, p_1, p_2 );
        fill_energy_gradient_hessian( e, E_out, dE_dp_out, d2E_dp2_out );
    }
    else
    {
        adrealG9_t e = Ecurvature_and_derivatives< adrealG9_t >( p_0, p_1, p_2 );
        fill_energy_gradient( e, E_out, dE_dp_out );
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
    
    DiffScalarBase::setVariableCount(9);
    
    if( d2E_dp2_out )
    {
        adrealH9_t e = Ecurvature_and_derivatives_undeformed< adrealH9_t >( p_0, p_1, p_2, u_0, u_1, u_2 );
        fill_energy_gradient_hessian( e, E_out, dE_dp_out, d2E_dp2_out );
    }
    else
    {
        adrealG9_t e = Ecurvature_and_derivatives_undeformed< adrealG9_t >( p_0, p_1, p_2, u_0, u_1, u_2 );
        fill_energy_gradient( e, E_out, dE_dp_out );
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
    
    DiffScalarBase::setVariableCount(9);
    
    if( d2E_dp2_out )
    {
        adrealH9_t e = Elaplacian_editing_and_derivatives< adrealH9_t >( p_0, p_1, p_2 );
        fill_energy_gradient_hessian( e, E_out, dE_dp_out, d2E_dp2_out );
    }
    else
    {
        adrealG9_t e = Elaplacian_editing_and_derivatives< adrealG9_t >( p_0, p_1, p_2 );
        fill_energy_gradient( e, E_out, dE_dp_out );
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
    
    DiffScalarBase::setVariableCount(9);
    
    if( d2E_dp2_out )
    {
        adrealH9_t e = Elaplacian_editing_and_derivatives_undeformed< adrealH9_t >( p_0, p_1, p_2, u_0, u_1, u_2 );
        fill_energy_gradient_hessian( e, E_out, dE_dp_out, d2E_dp2_out );
    }
    else
    {
        adrealG9_t e = Elaplacian_editing_and_derivatives_undeformed< adrealG9_t >( p_0, p_1, p_2, u_0, u_1, u_2 );
        fill_energy_gradient( e, E_out, dE_dp_out );
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
