/*
## TODO: Is there a setup.py I should be updating somewhere?

## On MacOS X:
cc face_stitch_pattern.c -dynamiclib -o face_stitch_pattern.dylib -std=c99 -g -O3 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare
*/

#include <assert.h>

// Output variable 'out' must be large enough to hold
// ( num_cross_section_verts ) * ( num_cross_sections - 1 ) * 6 ints.
void face_stitch_pattern( int num_cross_sections, int num_cross_section_verts, int* out )
{
    assert( out );
    assert( num_cross_sections >= 0 );
    assert( num_cross_section_verts >= 0 );
    
    for( int csi = 0; csi < num_cross_sections - 1; ++csi )
    for( int vi = 0; vi < num_cross_section_verts; ++vi )
    {
        const int v00 = csi * num_cross_section_verts + vi;
        const int v10 = csi * num_cross_section_verts + ( (vi + 1) % num_cross_section_verts );
        const int v01 = (csi + 1) * num_cross_section_verts + vi;
        const int v11 = (csi + 1) * num_cross_section_verts + ( (vi + 1) % num_cross_section_verts );
        
        *out++ = v00;
        *out++ = v10;
        *out++ = v11;
        
        *out++ = v00;
        *out++ = v11;
        *out++ = v01;
    }
}
