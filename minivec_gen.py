#!/usr/bin/env python

'''
Author: Yotam Gingold <yotam (strudel) yotamgingold.com>
License: Public Domain.  (I, Yotam Gingold, the author, dedicate any copyright to the Public Domain.)
         http://creativecommons.org/publicdomain/zero/1.0/
'''

def valsN_equals_chain( N ):
    base = 'vals[%d] ='
    return ' '.join( [ ( base % i ) for i in range(N) ] )

def base_repeat_like_array( base, N ):
    return [ base + '[' + str(i) + ']' for i in range(N) ]

def base_repeat_like_parens( base, N ):
    return [ base + '(' + str(i) + ')' for i in range(N) ]

def equals_sequence( lhs, rhs ):
    return operation_sequence( lhs, rhs, ' = ' )

def operation_sequence( lhs, rhs, op ):
    result = '; '.join( [ l + op + r for l,r in zip( lhs, rhs ) ] )
    result = result + ';'
    return result

def operation_sequence_joined( lhs, rhs, op, joinstr ):
    result = joinstr.join( [ l + op + r for l,r in zip( lhs, rhs ) ] )
    return result

def alphabet_sequence( N ):
    base = 'x' if N <= 3 else 'a'
    
    return ', '.join( [ 'const T& %s' % chr(ord(base) + i) for i in range(N) ] )

def alphabet_letters( N ):
    base = 'x' if N <= 3 else 'a'
    return ''.join( [ chr(ord(base) + i) for i in range(N) ] )

def gen( N ):
    result = '''#ifndef __minivecN_t_h__
#define __minivecN_t_h__

/*
Author: Yotam Gingold <yotam (strudel) yotamgingold.com>
License: Public Domain.  (I, Yotam Gingold, the author, dedicate any copyright to the Public Domain.)
         http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <cmath>
#include <iostream>

template< typename T >
class minivecN_t
{
public:
    // Initialize to zero.
    minivecN_t() { ''' + valsN_equals_chain( N ) + ''' 0; }
    // Copy constructor.
    template< typename U > minivecN_t( const minivecN_t<U>& cc ) { ''' + equals_sequence( base_repeat_like_array( 'vals', N ), base_repeat_like_array( 'cc.vals', N ) ) + '''}
    // Initialize from ''' + alphabet_letters( N ) + '''.
    minivecN_t( ''' + alphabet_sequence( N ) + ''' ) { ''' + equals_sequence( base_repeat_like_array( 'vals', N ), alphabet_letters( N ) ) + ''' }
    // Initialize from an arrray.
    minivecN_t( const T* ''' + alphabet_letters( N ) + ''' ) { ''' + equals_sequence( base_repeat_like_array( 'vals', N ), base_repeat_like_array( alphabet_letters( N ), N ) ) + ''' }
    
    // Silent cast to a pointer to the first element.
    // inline operator const T*() const { return &vals[0]; }
    
    // Get elements with parentheses:
    inline T& operator()( const unsigned int i ) { return vals[i]; }
    inline const T& operator()( const unsigned int i ) const { return vals[i]; };
    
    // Assignment operator.
    template< typename U > inline minivecN_t<T>& operator=( const minivecN_t<U>& rhs ) { ''' + equals_sequence( base_repeat_like_array( 'vals', N ), base_repeat_like_array( 'rhs.vals', N ) ) + ''' return *this; }
    
private:
    T vals[''' + str(N) + '''];
};

template< typename T > std::ostream& operator<<( std::ostream& out, const minivecN_t<T>& v ) { return out << ''' + operation_sequence_joined( ['']*N, base_repeat_like_parens( 'v', N ), '', " << ' ' << " ) + '''; }
template< typename T > std::istream& operator>>( std::istream& in, minivecN_t<T>& v ) { return in >> ''' + operation_sequence_joined( ['']*N, base_repeat_like_parens( 'v', N ), '', " >> " ) + '''; }

template< typename T >
inline minivecN_t<T> operator-( const minivecN_t<T>& lhs )
{
    return minivecN_t<T>( ''' + operation_sequence_joined( ['']*N, base_repeat_like_parens( 'lhs', N ), ' - ', ', ' ) + ''' );
}

template< typename T, typename U, typename V >
inline minivecN_t<T> operator-( const minivecN_t<U>& lhs, const minivecN_t<V>& rhs )
{
    return minivecN_t<T>( ''' + operation_sequence_joined( base_repeat_like_parens( 'lhs', N ), base_repeat_like_parens( 'rhs', N ), ' - ', ', ' ) + ''' );
}

template< typename T, typename U, typename V >
inline minivecN_t<T> operator+( const minivecN_t<U>& lhs, const minivecN_t<V>& rhs )
{
    return minivecN_t<T>( ''' + operation_sequence_joined( base_repeat_like_parens( 'lhs', N ), base_repeat_like_parens( 'rhs', N ), ' + ', ', ' ) + ''' );
}

template< typename T, typename U, typename V >
inline minivecN_t<T>& operator+=( minivecN_t<U>& lhs, const minivecN_t<V>& rhs )
{
    ''' + operation_sequence( base_repeat_like_parens( 'lhs', N ), base_repeat_like_parens( 'rhs', N ), ' += ' ) + '''
    return lhs;
}

template< typename T, typename U, typename V >
inline minivecN_t<T>& operator-=( minivecN_t<U>& lhs, const minivecN_t<V>& rhs )
{
    ''' + operation_sequence( base_repeat_like_parens( 'lhs', N ), base_repeat_like_parens( 'rhs', N ), ' -= ' ) + '''
    return lhs;
}

template< typename T, typename U, typename V >
inline minivecN_t<T> operator*( const minivecN_t<U>& v, const V& s )
{
    return minivecN_t<T>( ''' + operation_sequence_joined( base_repeat_like_parens( 'v', N ), ['s']*N, ' * ', ', ' ) + ''' );
}

template< typename T, typename U, typename V >
inline minivecN_t<T> operator*( const U s, const minivecN_t<V>& v )
{
    return v * s;
}

template< typename T, typename U, typename V >
inline minivecN_t<T> operator/( const minivecN_t<U>& v, const V& s )
{
    const V oos = 1./s;
    return v * oos;
}

template< typename T, typename U, typename V >
inline minivecN_t<T>& operator*=( minivecN_t<U>& v, const V& s )
{
    ''' + operation_sequence( base_repeat_like_parens( 'v', N ), ['s']*N, ' *= ' ) + '''
    return v;
}

template< typename T, typename U, typename V >
inline minivecN_t<T>& operator/=( minivecN_t<U>& v, const V& s )
{
    const V oos = 1./s;
    return v *= oos;
}

template< typename T >
inline T len( const minivecN_t<T>& v )
{
    return sqrt( ''' + operation_sequence_joined( base_repeat_like_parens( 'v', N ), base_repeat_like_parens( 'v', N ), '*', ' + ' ) + ''' );
}

template< typename T >
inline T lenSq( const minivecN_t<T>& v )
{
    return ''' + operation_sequence_joined( base_repeat_like_parens( 'v', N ), base_repeat_like_parens( 'v', N ), '*', ' + ' ) + ''';
}

template< typename T, typename U, typename V >
inline T dot( const minivecN_t<U>& v1, const minivecN_t<V>& v2 )
{
    return ''' + operation_sequence_joined( base_repeat_like_parens( 'v1', N ), base_repeat_like_parens( 'v2', N ), '*', ' + ' ) + ''';
}
''' + (
# There is a cross product in 3D
'''
template< typename T, typename U, typename V >
inline minivecN_t<T> cross( const minivecN_t<U>& v1, const minivecN_t<V>& v2 )
{
    return minivecN_t<T>(
        v1(1)*v2(2) - v1(2)*v2(1),
        v1(2)*v2(0) - v1(0)*v2(2),
        v1(0)*v2(1) - v1(1)*v2(0)
        );
}
'''
if N == 3
else
''
) + '''
// The compiler needs help inferring template types sometimes.
template< typename T > inline minivecN_t<T> operator-( const minivecN_t<T>& lhs, const minivecN_t<T>& rhs ) { return operator-<T,T,T>( lhs, rhs ); }
template< typename T > inline minivecN_t<T> operator+( const minivecN_t<T>& lhs, const minivecN_t<T>& rhs ) { return operator+<T,T,T>( lhs, rhs ); }
template< typename T > inline minivecN_t<T>& operator+=( minivecN_t<T>& lhs, const minivecN_t<T>& rhs ) { return operator+=<T,T,T>( lhs, rhs ); }
template< typename T > inline minivecN_t<T>& operator-=( minivecN_t<T>& lhs, const minivecN_t<T>& rhs ) { return operator-=<T,T,T>( lhs, rhs ); }
template< typename T > inline minivecN_t<T> operator*( const minivecN_t<T>& v, const T& s ) { return operator*<T,T,T>( v, s ); }
template< typename T > inline minivecN_t<T> operator*( const T& s, const minivecN_t<T>& v ) { return operator*<T,T,T>( s, v ); }
template< typename T, typename U > inline minivecN_t<T> operator*( const minivecN_t<T>& v, const U& s ) { return operator*<T,T,U>( v, s ); }
template< typename T, typename U > inline minivecN_t<T> operator*( const U& s, const minivecN_t<T>& v ) { return operator*<T,U,T>( s, v ); }
template< typename T > inline minivecN_t<T> operator/( const minivecN_t<T>& v, const T& s ) { return operator/<T,T,T>( v, s ); }
template< typename T > inline minivecN_t<T>& operator*=( minivecN_t<T>& v, const T& s ) { return operator*=<T,T,T>( v, s ); }
template< typename T, typename U > inline minivecN_t<T>& operator*=( minivecN_t<T>& v, const U& s ) { return operator*=<T,T,U>( v, s ); }
template< typename T > inline minivecN_t<T>& operator/=( minivecN_t<T>& v, const T& s ) { return operator/=<T,T,T>( v, s ); }
template< typename T > inline T dot( const minivecN_t<T>& v1, const minivecN_t<T>& v2 ) { return dot<T,T,T>( v1, v2 ); }
''' + (
# There is a cross product in 3D
'''template< typename T > inline minivecN_t<T> cross( const minivecN_t<T>& v1, const minivecN_t<T>& v2 ) { return cross<T,T,T>( v1, v2 ); }
'''
if N == 3
else
''
) + '''
#endif /* __minivecN_t_h__ */
'''
    result = result.replace( 'minivecN_t', 'minivec%d_t' % N )
    return result

def main():
    import sys
    
    def usage():
        print('Usage:', sys.argv[0], 'N', file=sys.stderr)
        print('Example:', sys.argv[0], '2 > minivec2_t.h', file=sys.stderr)
        print('Example:', sys.argv[0], '3 > minivec3_t.h', file=sys.stderr)
        sys.exit(-1)
    
    try:
        assert len( sys.argv ) == 2
        N = int( sys.argv[1] )
    except:
        usage()
    
    print(gen( N ))

if __name__ == '__main__': main()
