#ifndef __minivec2_t_h__
#define __minivec2_t_h__

/*
Author: Yotam Gingold <yotam (strudel) yotamgingold.com>
License: Public Domain.  (I, Yotam Gingold, the author, dedicate any copyright to the Public Domain.)
         http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <cmath>
#include <iostream>

template< typename T >
class minivec2_t
{
public:
    // Initialize to zero.
    minivec2_t() { vals[0] = vals[1] = 0; }
    // Copy constructor.
    template< typename U > minivec2_t( const minivec2_t<U>& cc ) { vals[0] = cc.vals[0]; vals[1] = cc.vals[1];}
    // Initialize from xy.
    minivec2_t( const T& x, const T& y ) { vals[0] = x; vals[1] = y; }
    // Initialize from an arrray.
    minivec2_t( const T* xy ) { vals[0] = xy[0]; vals[1] = xy[1]; }
    
    // Silent cast to a pointer to the first element.
    // inline operator const T*() const { return &vals[0]; }
    
    // Get elements with parentheses:
    inline T& operator()( const unsigned int i ) { return vals[i]; }
    inline const T& operator()( const unsigned int i ) const { return vals[i]; };
    
    // Assignment operator.
    template< typename U > inline minivec2_t<T>& operator=( const minivec2_t<U>& rhs ) { vals[0] = rhs.vals[0]; vals[1] = rhs.vals[1]; return *this; }
    
private:
    T vals[2];
};

template< typename T > std::ostream& operator<<( std::ostream& out, const minivec2_t<T>& v ) { return out << v(0) << ' ' << v(1); }
template< typename T > std::istream& operator>>( std::istream& in, minivec2_t<T>& v ) { return in >> v(0) >> v(1); }

template< typename T >
inline minivec2_t<T> operator-( const minivec2_t<T>& lhs )
{
    return minivec2_t<T>(  - lhs(0),  - lhs(1) );
}

template< typename T, typename U, typename V >
inline minivec2_t<T> operator-( const minivec2_t<U>& lhs, const minivec2_t<V>& rhs )
{
    return minivec2_t<T>( lhs(0) - rhs(0), lhs(1) - rhs(1) );
}

template< typename T, typename U, typename V >
inline minivec2_t<T> operator+( const minivec2_t<U>& lhs, const minivec2_t<V>& rhs )
{
    return minivec2_t<T>( lhs(0) + rhs(0), lhs(1) + rhs(1) );
}

template< typename T, typename U, typename V >
inline minivec2_t<T>& operator+=( minivec2_t<U>& lhs, const minivec2_t<V>& rhs )
{
    lhs(0) += rhs(0); lhs(1) += rhs(1);
    return lhs;
}

template< typename T, typename U, typename V >
inline minivec2_t<T>& operator-=( minivec2_t<U>& lhs, const minivec2_t<V>& rhs )
{
    lhs(0) -= rhs(0); lhs(1) -= rhs(1);
    return lhs;
}

template< typename T, typename U, typename V >
inline minivec2_t<T> operator*( const minivec2_t<U>& v, const V& s )
{
    return minivec2_t<T>( v(0) * s, v(1) * s );
}

template< typename T, typename U, typename V >
inline minivec2_t<T> operator*( const U s, const minivec2_t<V>& v )
{
    return v * s;
}

template< typename T, typename U, typename V >
inline minivec2_t<T> operator/( const minivec2_t<U>& v, const V& s )
{
    const V oos = 1./s;
    return v * oos;
}

template< typename T, typename U, typename V >
inline minivec2_t<T>& operator*=( minivec2_t<U>& v, const V& s )
{
    v(0) *= s; v(1) *= s;
    return v;
}

template< typename T, typename U, typename V >
inline minivec2_t<T>& operator/=( minivec2_t<U>& v, const V& s )
{
    const V oos = 1./s;
    return v *= oos;
}

template< typename T >
inline T len( const minivec2_t<T>& v )
{
    return sqrt( v(0)*v(0) + v(1)*v(1) );
}

template< typename T >
inline T lenSq( const minivec2_t<T>& v )
{
    return v(0)*v(0) + v(1)*v(1);
}

template< typename T, typename U, typename V >
inline T dot( const minivec2_t<U>& v1, const minivec2_t<V>& v2 )
{
    return v1(0)*v2(0) + v1(1)*v2(1);
}

// The compiler needs help inferring template types sometimes.
template< typename T > inline minivec2_t<T> operator-( const minivec2_t<T>& lhs, const minivec2_t<T>& rhs ) { return operator-<T,T,T>( lhs, rhs ); }
template< typename T > inline minivec2_t<T> operator+( const minivec2_t<T>& lhs, const minivec2_t<T>& rhs ) { return operator+<T,T,T>( lhs, rhs ); }
template< typename T > inline minivec2_t<T>& operator+=( minivec2_t<T>& lhs, const minivec2_t<T>& rhs ) { return operator+=<T,T,T>( lhs, rhs ); }
template< typename T > inline minivec2_t<T>& operator-=( minivec2_t<T>& lhs, const minivec2_t<T>& rhs ) { return operator-=<T,T,T>( lhs, rhs ); }
template< typename T > inline minivec2_t<T> operator*( const minivec2_t<T>& v, const T& s ) { return operator*<T,T,T>( v, s ); }
template< typename T > inline minivec2_t<T> operator*( const T& s, const minivec2_t<T>& v ) { return operator*<T,T,T>( s, v ); }
template< typename T, typename U > inline minivec2_t<T> operator*( const minivec2_t<T>& v, const U& s ) { return operator*<T,T,U>( v, s ); }
template< typename T, typename U > inline minivec2_t<T> operator*( const U& s, const minivec2_t<T>& v ) { return operator*<T,U,T>( s, v ); }
template< typename T > inline minivec2_t<T> operator/( const minivec2_t<T>& v, const T& s ) { return operator/<T,T,T>( v, s ); }
template< typename T > inline minivec2_t<T>& operator*=( minivec2_t<T>& v, const T& s ) { return operator*=<T,T,T>( v, s ); }
template< typename T, typename U > inline minivec2_t<T>& operator*=( minivec2_t<T>& v, const U& s ) { return operator*=<T,T,U>( v, s ); }
template< typename T > inline minivec2_t<T>& operator/=( minivec2_t<T>& v, const T& s ) { return operator/=<T,T,T>( v, s ); }
template< typename T > inline T dot( const minivec2_t<T>& v1, const minivec2_t<T>& v2 ) { return dot<T,T,T>( v1, v2 ); }

#endif /* __minivec2_t_h__ */

