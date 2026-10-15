#if 0
cc tritri.c tritri_wrapper.cpp raytri.c raytri_wrapper.cpp -lstdc++ \
    -dynamiclib -o libtritri.dylib \
    -g -O3 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare
#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cassert>

extern "C"
{
#include "tritri_wrapper.h"
}

namespace
{
std::vector<real_t> vs;
std::vector<int> faces;

void tritri_prepare()
{
    assert( vs.empty() == faces.empty() );
    // We're done if we've already loaded something.
    if( !vs.empty() ) return;
    
    // Load the OBJ file and parse it.
    std::ifstream sphere_obj( "./ext/sphere.obj" );
    
    std::string line, type;
    int v0, v1, v2;
    real_t x, y, z;
    
    while( !(sphere_obj >> std::ws).eof() )
    {
        getline( sphere_obj, line );
        std::istringstream linestr( line );
        
        linestr >> type;
        // We have to have some characters because we're not at eof().
        // Reading a word into a string should never fail.
        assert( linestr );
        assert( !type.empty() );
        if( type[0] == 'v' )
        {
            linestr >> x >> y >> z;
            if( !linestr ) throw "Invalid vertex";
            
            vs.push_back( x );
            vs.push_back( y );
            vs.push_back( z );
        }
        else if( type[0] == 'f' )
        {
            linestr >> v0 >> v1 >> v2;
            if( !linestr ) throw "Invalid face";
            
            if( v0 < 0 ) v0 = vs.size() + v0;
            else v0 -= 1;
            if( v1 < 0 ) v1 = vs.size() + v1;
            else v1 -= 1;
            if( v2 < 0 ) v2 = vs.size() + v2;
            else v2 -= 1;
            
            faces.push_back( v0 );
            faces.push_back( v1 );
            faces.push_back( v2 );
        }
    }
    
    assert( ( faces.size() % 3 ) == 0 );
    assert( ( vs.size() % 3 ) == 0 );
    for( int i = 0; i < faces.size(); ++i ) assert( faces[i] >= 0 && ( 3*faces[i] + 2 < vs.size() ) );
}

void transform_vertex( real_t center[3], real_t xform[9], real_t v[3], real_t v_out[3] )
{
    v_out[0] = center[0] + xform[0] * v[0] + xform[1] * v[1] + xform[2] * v[2];
    v_out[1] = center[1] + xform[3] * v[0] + xform[4] * v[1] + xform[5] * v[2];
    v_out[2] = center[2] + xform[6] * v[0] + xform[7] * v[1] + xform[8] * v[2];
}

}

extern "C"
{
int tri_tri_intersect_with_isectline(real_t V0[3],real_t V1[3],real_t V2[3],
				     real_t U0[3],real_t U1[3],real_t U2[3],int *coplanar,
				     real_t isectpt1[3],real_t isectpt2[3]);

real_t* tritri_wrapper_intersect_ellipsoids( real_t* center1, real_t* xform1, real_t* center2, real_t* xform2, int* num_lines_out )
{
    assert( center1 );
    assert( xform1 );
    assert( center2 );
    assert( xform2 );
    assert( num_lines_out );
    
    std::vector<real_t> line_segments;
    
    // Adapted from the python:
    /*
    for i in xrange( len(fs) ):
        for j in xrange( len(fs) ):
            intersect, coplanar, line_segment = _tri_tri_intersect_with_isectline(
                dot( xform1, vs[ fs[i][0] ] ) + center1,
                dot( xform1, vs[ fs[i][1] ] ) + center1,
                dot( xform1, vs[ fs[i][2] ] ) + center1,
                dot( xform2, vs[ fs[j][0] ] ) + center2,
                dot( xform2, vs[ fs[j][1] ] ) + center2,
                dot( xform2, vs[ fs[j][2] ] ) + center2
                )
            if intersect and not coplanar:
                result.append( line_segment )
    */
    
    tritri_prepare();
    assert( ( faces.size() % 3 ) == 0 );
    assert( ( vs.size() % 3 ) == 0 );
    const int n_faces = faces.size()/3;
    
    int intersect;
    int coplanar;
    real_t v0[3], v1[3], v2[3];
    real_t u0[3], u1[3], u2[3];
    real_t isectpt1[3], isectpt2[3];
    for( int i = 0; i < n_faces; ++i )
    for( int j = 0; j < n_faces; ++j )
    {
        transform_vertex( center1, xform1, &vs[ 3*faces[ i*3 + 0 ] ], v0 );
        transform_vertex( center1, xform1, &vs[ 3*faces[ i*3 + 1 ] ], v1 );
        transform_vertex( center1, xform1, &vs[ 3*faces[ i*3 + 2 ] ], v2 );
        
        transform_vertex( center2, xform2, &vs[ 3*faces[ j*3 + 0 ] ], u0 );
        transform_vertex( center2, xform2, &vs[ 3*faces[ j*3 + 1 ] ], u1 );
        transform_vertex( center2, xform2, &vs[ 3*faces[ j*3 + 2 ] ], u2 );
        
        intersect = tri_tri_intersect_with_isectline(
            v0, v1, v2,
            u0, u1, u2,
            &coplanar,
            isectpt1, isectpt2
            );
        if( intersect && !coplanar )
        {
            line_segments.push_back( isectpt1[0] );
            line_segments.push_back( isectpt1[1] );
            line_segments.push_back( isectpt1[2] );
            
            line_segments.push_back( isectpt2[0] );
            line_segments.push_back( isectpt2[1] );
            line_segments.push_back( isectpt2[2] );
        }
    }
    
    real_t* result = new real_t[ line_segments.size() ];
    std::copy( line_segments.begin(), line_segments.end(), result );
    assert( ( line_segments.size() % 6 ) == 0 );
    *num_lines_out = line_segments.size() / 6;
    return result;
}

real_t* tritri_wrapper_intersect_meshes(
    real_t* vs0, int num_vertices0, int* faces0, int num_faces0,
    real_t* vs1, int num_vertices1, int* faces1, int num_faces1,
    int* num_lines_out
    )
{
    assert( num_vertices0 > 0 );
    assert( num_vertices1 > 0 );
    assert( num_faces0 > 0 );
    assert( num_faces1 > 0 );
    assert( vs0 && faces0 && vs1 && faces1 );
    assert( num_lines_out );
    
    std::vector<real_t> line_segments;
    
    // In python (untested):
    /*
    for i in xrange( len(faces0) ):
        for j in xrange( len(faces1) ):
            intersect, coplanar, line_segment = _tri_tri_intersect_with_isectline(
                vs0[ faces0[i][0] ] ),
                vs0[ faces0[i][1] ] ),
                vs0[ faces0[i][2] ] ),
                
                vs1[ faces1[j][0] ] ),
                vs1[ faces1[j][1] ] ),
                vs1[ faces1[j][2] ] )
                )
            if intersect and not coplanar:
                result.append( line_segment )
    */
    
    int intersect;
    int coplanar;
    real_t isectpt1[3], isectpt2[3];
    for( int i = 0; i < num_faces0; ++i )
    for( int j = 0; j < num_faces1; ++j )
    {
        intersect = tri_tri_intersect_with_isectline(
            vs0 + 3*faces0[ 3*i + 0 ],
            vs0 + 3*faces0[ 3*i + 1 ],
            vs0 + 3*faces0[ 3*i + 2 ],
            
            vs1 + 3*faces1[ 3*j + 0 ],
            vs1 + 3*faces1[ 3*j + 1 ],
            vs1 + 3*faces1[ 3*j + 2 ],
            
            &coplanar,
            isectpt1, isectpt2
            );
        if( intersect && !coplanar )
        {
            line_segments.push_back( isectpt1[0] );
            line_segments.push_back( isectpt1[1] );
            line_segments.push_back( isectpt1[2] );
            
            line_segments.push_back( isectpt2[0] );
            line_segments.push_back( isectpt2[1] );
            line_segments.push_back( isectpt2[2] );
        }
    }
    
    real_t* result = new real_t[ line_segments.size() ];
    std::copy( line_segments.begin(), line_segments.end(), result );
    assert( ( line_segments.size() % 6 ) == 0 );
    *num_lines_out = line_segments.size() / 6;
    return result;
}

face_face_intersection_t* tritri_wrapper_intersect_meshes_info(
    real_t* vs0, int num_vertices0, int* faces0, int num_faces0,
    real_t* vs1, int num_vertices1, int* faces1, int num_faces1,
    int* num_lines_out
    )
{
    assert( num_vertices0 > 0 );
    assert( num_vertices1 > 0 );
    assert( num_faces0 > 0 );
    assert( num_faces1 > 0 );
    assert( vs0 && faces0 && vs1 && faces1 );
    assert( num_lines_out );
    
    std::vector<real_t> line_segments;
    std::vector<int> line_segment_faces;
    
    // In python (untested):
    /*
    for i in xrange( len(faces0) ):
        for j in xrange( len(faces1) ):
            intersect, coplanar, line_segment = _tri_tri_intersect_with_isectline(
                vs0[ faces0[i][0] ] ),
                vs0[ faces0[i][1] ] ),
                vs0[ faces0[i][2] ] ),
                
                vs1[ faces1[j][0] ] ),
                vs1[ faces1[j][1] ] ),
                vs1[ faces1[j][2] ] )
                )
            if intersect and not coplanar:
                result_lines.append( line_segment )
                result_faces.append( ( i, j ) )
    */
    
    int intersect;
    int coplanar;
    real_t isectpt1[3], isectpt2[3];
    for( int i = 0; i < num_faces0; ++i )
    for( int j = 0; j < num_faces1; ++j )
    {
        intersect = tri_tri_intersect_with_isectline(
            vs0 + 3*faces0[ 3*i + 0 ],
            vs0 + 3*faces0[ 3*i + 1 ],
            vs0 + 3*faces0[ 3*i + 2 ],
            
            vs1 + 3*faces1[ 3*j + 0 ],
            vs1 + 3*faces1[ 3*j + 1 ],
            vs1 + 3*faces1[ 3*j + 2 ],
            
            &coplanar,
            isectpt1, isectpt2
            );
        if( intersect && !coplanar )
        {
            line_segments.push_back( isectpt1[0] );
            line_segments.push_back( isectpt1[1] );
            line_segments.push_back( isectpt1[2] );
            
            line_segments.push_back( isectpt2[0] );
            line_segments.push_back( isectpt2[1] );
            line_segments.push_back( isectpt2[2] );
            
            line_segment_faces.push_back( i );
            line_segment_faces.push_back( j );
        }
    }
    
    assert( ( line_segments.size() % 6 ) == 0 );
    assert( ( line_segment_faces.size() % 2 ) == 0 );
    assert( line_segments.size() / 6 == line_segment_faces.size() / 2 );
    *num_lines_out = line_segments.size() / 6;
    
    face_face_intersection_t* result = new face_face_intersection_t[ *num_lines_out ];
    for( int i = 0; i < *num_lines_out; ++i )
    {
        result[i].x0 = line_segments[ 6*i + 0 ];
        result[i].y0 = line_segments[ 6*i + 1 ];
        result[i].z0 = line_segments[ 6*i + 2 ];
        
        result[i].x1 = line_segments[ 6*i + 3 ];
        result[i].y1 = line_segments[ 6*i + 4 ];
        result[i].z1 = line_segments[ 6*i + 5 ];
        
        result[i].face0 = line_segment_faces[ 2*i + 0 ];
        result[i].face1 = line_segment_faces[ 2*i + 1 ];
    }
    
    return result;
}

void tritri_wrapper_free( real_t* memory )
{
    delete [] memory;
}

void tritri_wrapper_free_info( face_face_intersection_t* memory )
{
    delete [] memory;
}

}
