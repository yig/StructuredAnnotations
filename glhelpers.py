from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

from myarray import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.

def draw_lineloop( points ):
    '''
    Given a list of 2d or 3d tuples 'points', draw them as a line loop (an implicit segment between
    the last and first points).
    '''
    if len( points ) < 2: return
    
    ## Q: Do I need to do this?
    # points = asarray( points )
    ## A1: If I convert all input into a list of tuples,
    # points = [ tuple( pt ) for pt in points ]
    ##     everything still runs fine.  I'll let PyOpenGL do the conversion since it knows what's
    ##     best for it.
    ## A2: The latest PyOpenGL doesn't take lists, so I always convert to arrays.
    
    glEnableClientState( GL_VERTEX_ARRAY )
    glVertexPointer( len(points[0]), GL_DOUBLE, 0, asarray( points ) )
    glDrawArrays( GL_LINE_LOOP, 0, len(points) )
    glDisableClientState( GL_VERTEX_ARRAY )

def draw_linestrip( points ):
    '''
    Given a list of 2d or 3d tuples 'points', draw them as a line strip.
    '''
    if len( points ) < 2: return
    
    ## Q: Do I need to do this?
    # points = asarray( points )
    ## A1: If I convert all input into a list of tuples,
    # points = [ tuple( pt ) for pt in points ]
    ##     everything still runs fine.  I'll let PyOpenGL do the conversion since it knows what's
    ##     best for it.
    ## A2: The latest PyOpenGL doesn't take lists, so I always convert to arrays.
    
    glEnableClientState( GL_VERTEX_ARRAY )
    glVertexPointer( len(points[0]), GL_DOUBLE, 0, asarray( points ) )
    glDrawArrays( GL_LINE_STRIP, 0, len(points) )
    glDisableClientState( GL_VERTEX_ARRAY )

def draw_thick_linestrip2d( points, thickness, colors = None ):
    '''
    Given a list of 2d 'points', draw them as a line strip with given diameter 'thickness'.
    
    Also takes an optional array of floating point color values 'colors', where
    'colors' is an array of the same length as 'points', each element in 'colors'
    has either 3 (RGB) or 4 (RGBA) elements, and each element is in the range [0..1].
    If 'colors' is not specified, the OpenGL color state is not modified and you can
    call glColor*() before calling this function.
    '''
    if len( points ) < 2: return
    
    ## Q: Do I need to do this?
    # points = asarray( points )
    ## A1: If I convert all input into a list of tuples,
    # points = [ tuple( pt ) for pt in points ]
    ##     everything still runs fine.  I'll let PyOpenGL do the conversion since it knows what's
    ##     best for it.
    
    #assert len( points[0] ) == 2
    points = asarrayf( points )
    
    from helpers import tangents_and_normals_along_line_strip_2d, normalize_vectors
    ## Note that these edge normals are a right turn along the path.
    unit_edges, enormals = tangents_and_normals_along_line_strip_2d( points )
    ## Average normals (and repeat the first and last).
    point_normals = concatenate( ( enormals[:1], enormals[:-1] + enormals[1:], enormals[-1:] ) )
    ## Re-Normalize.
    point_normals[:,2:] = 0.
    point_normals = normalize_vectors( point_normals )
    
    offsets = .5*thickness * point_normals
    thick_points = list(zip( points + offsets, points - offsets ))
    
    if colors is not None:
        colors = asarray( colors, dtype = float32 )
        thick_colors = list(zip( colors, colors ))
        assert len( colors.shape ) == 2
        assert colors.shape[0] == len( points )
        assert colors.shape[1] in (3,4)
        
        glEnableClientState( GL_COLOR_ARRAY )
        glColorPointer( len(colors[0]), GL_FLOAT, 0, asarray( thick_colors ) )
    
    glEnableClientState( GL_VERTEX_ARRAY )
    glVertexPointer( len(points[0]), GL_DOUBLE, 0, asarray( thick_points ) )
    glDrawArrays( GL_TRIANGLE_STRIP, 0, 2*len(thick_points) )
    glDisableClientState( GL_VERTEX_ARRAY )
    
    if colors is not None:
        glDisableClientState( GL_COLOR_ARRAY )

def draw_thick_lineloop2d( points, thickness, colors = None ):
    '''
    Given a list of 2d 'points', draw them as a line loop with given diameter 'thickness'.
    
    Also takes an optional array of floating point color values 'colors', where
    'colors' is an array of the same length as 'points', each element in 'colors'
    has either 3 (RGB) or 4 (RGBA) elements, and each element is in the range [0..1].
    If 'colors' is not specified, the OpenGL color state is not modified and you can
    call glColor*() before calling this function.
    '''
    if len( points ) < 2: return
    
    ## Q: Do I need to do this?
    # points = asarray( points )
    ## A1: If I convert all input into a list of tuples,
    # points = [ tuple( pt ) for pt in points ]
    ##     everything still runs fine.  I'll let PyOpenGL do the conversion since it knows what's
    ##     best for it.
    
    #assert len( points[0] ) == 2
    points = append( points, points[:1], 0 )
    
    from helpers import tangents_and_normals_along_line_strip_2d, normalize_vectors
    ## Note that these edge normals are a right turn along the path.
    unit_edges, enormals = tangents_and_normals_along_line_strip_2d( points[:-1] )
    ## Average normals, including the first with the last.
    ## Repeat the first/last average once at the beginning and twice at the end, since this
    ## is a line loop.
    middle_normals = enormals[:-1] + enormals[1:]
    first_last_normal = enormals[:1] + enormals[-1:]
    point_normals = concatenate( ( first_last_normal, middle_normals, first_last_normal, first_last_normal ) )
    ## Re-Normalize.
    point_normals[:,2:] = 0.
    point_normals = normalize_vectors( point_normals )
    
    offsets = .5*thickness * point_normals
    thick_points = list(zip( points + offsets, points - offsets ))
    
    if colors is not None:
        colors = asarray( colors, dtype = float32 )
        colors = append( colors, colors[:1], 0 )
        thick_colors = list(zip( colors, colors ))
        assert len( colors.shape ) == 2
        assert colors.shape[0] == len( points )
        assert colors.shape[1] in (3,4)
        
        glEnableClientState( GL_COLOR_ARRAY )
        glColorPointer( len(colors[0]), GL_FLOAT, 0, asarray( thick_colors ) )
    
    glEnableClientState( GL_VERTEX_ARRAY )
    glVertexPointer( len(points[0]), GL_DOUBLE, 0, asarray( thick_points ) )
    glDrawArrays( GL_TRIANGLE_STRIP, 0, 2*len(thick_points) )
    glDisableClientState( GL_VERTEX_ARRAY )
    
    if colors is not None:
        glDisableClientState( GL_COLOR_ARRAY )

def draw_lines( points ):
    '''
    Given a list of 2d or 3d tuples 'points', draw them as a series of lines.
    Each successive pair of points is a line.
    '''
    if len( points ) < 2: return
    
    ## Q: Do I need to do this?
    # points = asarray( points )
    ## A1: If I convert all input into a list of tuples,
    # points = [ tuple( pt ) for pt in points ]
    ##     everything still runs fine.  I'll let PyOpenGL do the conversion since it knows what's
    ##     best for it.
    
    glEnableClientState( GL_VERTEX_ARRAY )
    glVertexPointer( len(points[0]), GL_DOUBLE, 0, asarray( points ) )
    glDrawArrays( GL_LINES, 0, len(points) )
    glDisableClientState( GL_VERTEX_ARRAY )

def set_linear_alpha_for_points( points_values, base_color3f, min_alpha = 0., max_alpha = 1., min_value = None, max_value = None ):
    '''
    This function is for use with the draw_line*() functions.
    Call it before calling one of them.
    Then call the matching unset_linear_alpha_for_points() function.
    
    It sets up N colors for N points, based on the N values in 'points_values'.
    Specifically, the minimum value in 'points_values' (or 'min_value', if set) maps to 'min_alpha'
    and the maximum value (or 'max_value', if set) to 'max_alpha'.
    The rest of the values in 'points_values' map linearly.
    The alpha is used with the base color 'base_color3f' (3 floating point values between 0 and 1).
    '''
    real_t = float32
    
    colors = colors_from_linear_alpha_for_points( points_values, base_color3f, min_alpha, max_alpha, min_value, max_value )
    colors = asarray( colors, dtype = real_t )
    assert colors.shape[0] == len( points_values )
    assert colors.shape[1] == 4
    assert len( colors.shape ) == 2
    
    glEnableClientState( GL_COLOR_ARRAY )
    glColorPointer( len(colors[0]), GL_FLOAT, 0, asarray( colors ) )

def colors_from_linear_alpha_for_points( points_values, base_color3f, min_alpha = 0., max_alpha = 1., min_value = None, max_value = None ):
    '''
    This function is called inside of set_linear_alpha_for_points().
    
    It returns N colors for N points, based on the N values in 'points_values'.
    Specifically, the minimum value in 'points_values' (or 'min_value', if set) maps to 'min_alpha'
    and the maximum value (or 'max_value', if set) to 'max_alpha'.
    The rest of the values in 'points_values' map linearly.
    The alpha is used with the base color 'base_color3f' (3 floating point values between 0 and 1).
    '''
    real_t = float32
    
    assert len( base_color3f ) == 3
    base_color3f = list(map( real_t, base_color3f ))
    
    colors = empty( ( len(points_values), 4 ), dtype = real_t )
    colors[:,0:3] = base_color3f
    
    points_values = asarray( points_values, dtype = real_t )
    min_value = min_value if ( min_value is not None ) else points_values.min()
    max_value = max_value if ( max_value is not None ) else points_values.max()
    if max_value == min_value: colors[:,3] = max_alpha
    else: colors[:,3] = min_alpha + max_alpha * ( (points_values - min_value) / (max_value - min_value) )
    
    return colors

def unset_linear_alpha_for_points():
    '''
    See set_linear_alpha_for_points().
    '''
    glDisableClientState( GL_COLOR_ARRAY )

_draw_unit_circle_cache = {}
def _draw_unit_circle( start_radians, end_radians, NUM_SLICES ):
    key = ( start_radians, end_radians, NUM_SLICES )
    try:
        glCallList( _draw_unit_circle_cache[ key ] )
    except KeyError:
        list_id = glGenLists( 1 )
        _draw_unit_circle_cache[ key ] = list_id
        glNewList( list_id, GL_COMPILE_AND_EXECUTE )
        
        samples = linspace( start_radians, end_radians, NUM_SLICES )
        points = arrayf( ( cos( samples ), sin( samples ) ) ).T
        
        glEnableClientState( GL_VERTEX_ARRAY )
        glVertexPointer( 2, GL_DOUBLE, 0, asarray( points ) )
        glDrawArrays( GL_LINE_STRIP, 0, len(points) )
        glDisableClientState( GL_VERTEX_ARRAY )
        
        glEndList()

_draw_unit_circle_interior_cache = {}
def _draw_unit_circle_interior( start_radians, end_radians, NUM_SLICES ):
    key = ( start_radians, end_radians, NUM_SLICES )
    try:
        glCallList( _draw_unit_circle_interior_cache[ key ] )
    except KeyError:
        list_id = glGenLists( 1 )
        _draw_unit_circle_interior_cache[ key ] = list_id
        glNewList( list_id, GL_COMPILE_AND_EXECUTE )
        
        samples = linspace( start_radians, end_radians, NUM_SLICES )
        points = emptyf( ( NUM_SLICES + 1, 2 ) )
        ## The first point of the triangle fan
        points[0,:] = 0.
        ## The rest of the points
        points[1:,0] = cos( samples )
        points[1:,1] = sin( samples )
        
        glEnableClientState( GL_VERTEX_ARRAY )
        glVertexPointer( 2, GL_DOUBLE, 0, asarray( points ) )
        glDrawArrays( GL_TRIANGLE_FAN, 0, len(points) )
        glDisableClientState( GL_VERTEX_ARRAY )
        
        glEndList()

def draw_circle( center, radius, NUM_SLICES = 40, start_radians = 0, end_radians = 2*pi ):
    '''
    Given a center and radius for a circle, draw the circle.
    NOTE: This does not draw a disc, it draws a circle, which is the outline of a disc.
    
    Optional parameters are:
        NUM_SLICES, which determines the number of points to sample along the ellipse (default 40)
        start_radians, which determines the starting angle in radians (default 0)
        end_radians, which determines the ending angle in radians (default 2*pi)
    NOTE: Angles 0, pi/2, pi, ... correspond to +axis1, +axis2, -axis1, ...
    '''
    
    mmode = glGetInteger( GL_MATRIX_MODE )
    if mmode != GL_MODELVIEW: glMatrixMode( GL_MODELVIEW )
    glPushMatrix()
    
    center3 = zeros( 3 )
    center3[:len(center)] = center
    glTranslated( *center3 )
    
    glScaled( radius, radius, 1. )
    
    _draw_unit_circle( start_radians, end_radians, NUM_SLICES )
    
    glPopMatrix()
    if mmode != GL_MODELVIEW: glMatrixMode( mmode )

def draw_circle_interior( center, radius, NUM_SLICES = 40, start_radians = 0, end_radians = 2*pi ):
    '''
    Given a center and radius for a circle, draw the interior disc.
    NOTE: This does not draw a circle, it draws a disc, which is the interior of a circle.
    
    Optional parameters are:
        NUM_SLICES, which determines the number of points to sample along the ellipse (default 40)
        start_radians, which determines the starting angle in radians (default 0)
        end_radians, which determines the ending angle in radians (default 2*pi)
    NOTE: Angles 0, pi/2, pi, ... correspond to +axis1, +axis2, -axis1, ...
    '''
    
    mmode = glGetInteger( GL_MATRIX_MODE )
    if mmode != GL_MODELVIEW: glMatrixMode( GL_MODELVIEW )
    glPushMatrix()
    
    center3 = zeros( 3 )
    center3[:len(center)] = center
    glTranslated( *center3 )
    
    glScaled( radius, radius, 1. )
    
    _draw_unit_circle_interior( start_radians, end_radians, NUM_SLICES )
    
    glPopMatrix()
    if mmode != GL_MODELVIEW: glMatrixMode( mmode )

def draw_ellipse( center, axis1, axis2, NUM_SLICES = 40, start_radians = 0, end_radians = 2*pi ):
    '''
    Given a center point 'center' and two axes, 'axis1' and 'axis2', draw the
    ellipse which has 'axis1' and 'axis2' as its axes.
    NOTE: This draws a closed loop in space, not a filled disc.
    
    Optional parameters are:
        NUM_SLICES, which determines the number of points to sample along the ellipse (default 40)
        start_radians, which determines the starting angle in radians (default 0)
        end_radians, which determines the ending angle in radians (default 2*pi)
    NOTE: Angles 0, pi/2, pi, ... correspond to +axis1, +axis2, -axis1, ...
    '''
    
    mmode = glGetInteger( GL_MATRIX_MODE )
    if mmode != GL_MODELVIEW: glMatrixMode( GL_MODELVIEW )
    glPushMatrix()
    
    center3 = zeros( 3 )
    center3[:len(center)] = center
    glTranslated( *center3 )
    
    rot_scale_mat = identityf( 4 )
    rot_scale_mat[ :len(axis1), 0 ] = axis1
    rot_scale_mat[ :len(axis2), 1 ] = axis2
    glMultMatrixd( rot_scale_mat.T.ravel() )
    
    _draw_unit_circle( start_radians, end_radians, NUM_SLICES )
    
    glPopMatrix()
    if mmode != GL_MODELVIEW: glMatrixMode( mmode )

def draw_ellipse_interior( center, axis1, axis2, NUM_SLICES = 40, start_radians = 0, end_radians = 2*pi ):
    '''
    Given a center point 'center' and two axes, 'axis1' and 'axis2', draw the
    filled interior of the ellipse which has 'axis1' and 'axis2' as its axes.
    
    Optional parameters are:
        NUM_SLICES, which determines the number of points to sample along the ellipse (default 40)
        start_radians, which determines the starting angle in radians (default 0)
        end_radians, which determines the ending angle in radians (default 2*pi)
    NOTE: Angles 0, pi/2, pi, ... correspond to +axis1, +axis2, -axis1, ...
    '''
    
    mmode = glGetInteger( GL_MATRIX_MODE )
    if mmode != GL_MODELVIEW: glMatrixMode( GL_MODELVIEW )
    glPushMatrix()
    
    center3 = zeros( 3 )
    center3[:len(center)] = center
    glTranslated( *center3 )
    
    rot_scale_mat = identityf( 4 )
    rot_scale_mat[ :len(axis1), 0 ] = axis1
    rot_scale_mat[ :len(axis2), 1 ] = axis2
    glMultMatrixd( rot_scale_mat.T.ravel() )
    
    _draw_unit_circle_interior( start_radians, end_radians, NUM_SLICES )
    
    glPopMatrix()
    if mmode != GL_MODELVIEW: glMatrixMode( mmode )

def draw_ellipse_crust(
    center, axis1, axis2, radius1_larger, radius2_larger, radius1_smaller, radius2_smaller,
    larger_color = None,
    smaller_color = None,
    NUM_SLICES = 40, start_radians = 0, end_radians = 2*pi
    ):
    '''
    Given a center point 'center' and two axes, 'axis1' and 'axis2', draw a filled
    crust of ellipse which has larger radii 'radius1/2_larger' and smaller radii 'radius1/2_smaller'.
    
    Optional parameters are:
        larger_color, which passes a color for points on the larger ellipse (default no color is set)
        smaller_color, which passes a color for points on the smaller ellipse (default no color is set)
        NUM_SLICES, which determines the number of points to sample along the ellipse (default 40)
        start_radians, which determines the starting angle in radians (default 0)
        end_radians, which determines the ending angle in radians (default 2*pi)
    NOTE: Angles 0, pi/2, pi, ... correspond to +axis1, +axis2, -axis1, ...
    '''
    
    assert (larger_color is None) == (smaller_color is None)
    assert radius1_larger > radius1_smaller
    assert radius2_larger > radius2_smaller
    assert radius1_larger >= 0
    assert radius2_larger >= 0
    assert radius1_smaller >= 0
    assert radius2_smaller >= 0
    
    mmode = glGetInteger( GL_MATRIX_MODE )
    if mmode != GL_MODELVIEW: glMatrixMode( GL_MODELVIEW )
    glPushMatrix()
    
    center3 = zeros( 3 )
    center3[:len(center)] = center
    glTranslated( *center3 )
    
    rot_scale_mat = identityf( 4 )
    rot_scale_mat[ :len(axis1), 0 ] = dir( axis1 )
    rot_scale_mat[ :len(axis2), 1 ] = dir( axis2 )
    glMultMatrixd( rot_scale_mat.T.ravel() )
    
    samples = linspace( start_radians, end_radians, NUM_SLICES )
    cos_samples = cos( samples )
    sin_samples = sin( samples )
    
    points = emptyf( ( NUM_SLICES, 4 ) )
    ## Points on the larger ellipse.
    points[:,0] = radius1_larger * cos_samples
    points[:,1] = radius2_larger * sin_samples
    ## Points on the smaller ellipse.
    points[:,2] = radius1_smaller * cos_samples
    points[:,3] = radius2_smaller * sin_samples
    
    if larger_color is not None:
        colors = onesf( ( NUM_SLICES, 8 ) )
        colors[ :, :len(larger_color) ] = larger_color
        colors[ :, 4:4+len(smaller_color) ] = smaller_color
        
        glEnableClientState( GL_COLOR_ARRAY )
        glColorPointer( 4, GL_DOUBLE, 0, asarray( colors ) )
    
    glEnableClientState( GL_VERTEX_ARRAY )
    glVertexPointer( 2, GL_DOUBLE, 0, asarray( points ) )
    glDrawArrays( GL_TRIANGLE_STRIP, 0, 2*len(points) )
    glDisableClientState( GL_VERTEX_ARRAY )
    
    if larger_color is not None: glDisableClientState( GL_COLOR_ARRAY )
    
    glPopMatrix()
    if mmode != GL_MODELVIEW: glMatrixMode( mmode )

def draw_text( pos, text, font = GLUT_BITMAP_HELVETICA_18 ):
    '''
    Given a 2d or 3d position 'pos' and string of text to draw 'text',
    draws the string at the position using the given 'font' (defaults to GLUT_BITMAP_HELVETICA_18).
    '''
    glRasterPos = glRasterPos2f if len(pos) == 2 else glRasterPos3f
    
    glRasterPos( *pos )
    for c in text: glutBitmapCharacter( font, ord(c) )

## It seems glutBitmapLength is broken, but it is defined to be:
def glutBitmapLength( font, text ): return sum( [ glutBitmapWidth( font, ord(c) ) for c in text ] )

def draw_text_pixel_size( text, font = GLUT_BITMAP_HELVETICA_18 ):
    '''
    Given a string of text to draw 'text', returns the 2D (width,height) in pixels
    of the string using the given 'font' (defaults to GLUT_BITMAP_HELVETICA_18).
    '''
    
    heights = {
        GLUT_BITMAP_8_BY_13.value: 13,
        GLUT_BITMAP_9_BY_15.value: 15,
        GLUT_BITMAP_TIMES_ROMAN_10.value: 10,
        GLUT_BITMAP_TIMES_ROMAN_24.value: 24,
        GLUT_BITMAP_HELVETICA_10.value: 10,
        GLUT_BITMAP_HELVETICA_12.value: 12,
        GLUT_BITMAP_HELVETICA_18.value: 18
        }
    
    return ( heights[ font.value ], glutBitmapLength( font, text ) )

def draw_text_world_size( world2pixels, pixels2world, text, font = GLUT_BITMAP_HELVETICA_18 ):
    '''
    Given a string of text to draw 'text', returns the 2D (width,height) in world coordinates
    of the string using the given 'font' (defaults to GLUT_BITMAP_HELVETICA_18).
    '''
    
    pixel_size = draw_text_pixel_size( text, font )
    
    pixels_vec = list( world2pixels( (0,0) ) )
    pixels_vec[0] += pixel_size[0]
    pixels_vec[1] += pixel_size[1]
    return pixels2world( pixels_vec )

def draw_gltexture( gltex, tex_pos, tex_size ):
    '''
    Given a GLTexture gltex and a 2d position and 2d size, draws the texture modulated
    with the current color.
    '''
    
    tex_pos = asarray( tex_pos )
    tex_size = asarray( tex_size )
    
    glBindTexture( GL_TEXTURE_2D, gltex.id )
    
    tex_was_off = not glIsEnabled( GL_TEXTURE_2D )
    if tex_was_off: glEnable( GL_TEXTURE_2D )
    
    env_mode = glGetTexEnviv( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE )
    if env_mode != GL_MODULATE: glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE )
    
    xone = gltex.actual_width
    yone = gltex.actual_height
    
    glBegin( GL_QUADS )
    glTexCoord2f( 0,0 )
    glVertex2f( *tex_pos )
    
    glTexCoord2f( xone,0 )
    glVertex2f( *( tex_pos + array((tex_size[0],0)) ) )
    
    glTexCoord2f( xone,yone )
    glVertex2f( *( tex_pos + tex_size ) )
    
    glTexCoord2f( 0,yone )
    glVertex2f( *( tex_pos + array((0,tex_size[1])) ) )
    glEnd()
    
    if env_mode != GL_MODULATE: glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, env_mode )
    if tex_was_off: glDisable( GL_TEXTURE_2D )

'''
def set_blend_mode( mode ):
    if mode == 'over':
        ## see
        ##   image_overlap.cpp,
        ##   GPU ImageProcessing.pdf from: <ftp://download.nvidia.com/developer/presentations/2004/SIGGRAPH/ImageProcessing.pdf>
        ##   http://en.wikipedia.org/wiki/Alpha_compositing
        ## I need to assume glColor4 premultiplies, and that textures are premultiplied, and that
        ## the result in the framebuffer is composited with black (un-premultiplied) before showing
        ## on-screen.
        glBlendFuncSeparate( ... )
    else:
        raise NotImplementedError, 'this function is a work-in-progress'
'''

class NoShaderProgram( object ):
    def __init__( self ):
        self.program = 0
    def UseProgram( self ): pass
    def UnUseProgram( self ): pass

class ShaderProgram( object ):
    def __init__( self ):
        extensions = glGetString( GL_EXTENSIONS ).decode("utf-8", "strict").split(" ")
        version = glGetString( GL_VERSION ).decode("utf-8", "strict").split(".")
        ## TODO Q: Does this break if the GL version is less than 2 but the GL_ARB extension is
        ##         present?  Should I use alternates as in "OpenGLContext/tests/shaderobject.py"?
        if not( int(version[0]) >= 2 or 'GL_ARB_shader_objects' in extensions ):
            raise NotImplementedError("Shaders are not implemented on this platform.")
        
        self.program = glCreateProgram()
        ## This doesn't work, it takes lots of arguments like the C version.
        ## Fortunately, someone (the python bindings?) are printing out lots of information
        ## when a shader fails.
        #info_log = glGetProgramInfoLog( self.program )
        #print info_log
        
        self.vertex_shader = None
        self.fragment_shader = None
    
    def Destroy( self ):
        self.UnloadVertexShader()
        self.UnloadFragmentShader()
        glDeleteProgram( self.program )
        self.program = None
    
    def LoadVertexShaderFromFilename( self, filename ):
        self.LoadVertexShaderFromString( open( filename ).read() )
    def LoadVertexShaderFromString( self, string ):
        assert self.program is not None
        
        self.vertex_shader = self.shader_from_string( GL_VERTEX_SHADER, string, self.vertex_shader )
        ## This doesn't work, it takes lots of arguments like the C version.
        ## Fortunately, someone (the python bindings?) are printing out lots of information
        ## when a shader fails.
        #info_log = glGetShaderInfoLog( self.vertex_shader )
        #print info_log
        
        glAttachShader( self.program, self.vertex_shader )
        glLinkProgram( self.program )
    
    def UnloadVertexShader( self ):
        assert self.program is not None
        
        if self.vertex_shader is None:
            print("No vertex program loaded.")
            return
        
        glDetachShader( self.program, self.vertex_program )
        glLinkProgram( self.program )
        
        glDeleteShader( self.vertex_program )
        self.vertex_shader = None
    
    def LoadFragmentShaderFromFile( self, filename ):
        self.LoadFragmentShaderFromString( open( filename ).read() )
    def LoadFragmentShaderFromString( self, string ):
        assert self.program is not None
        
        self.fragment_shader = self.shader_from_string( GL_FRAGMENT_SHADER, string, self.fragment_shader )
        ## This doesn't work, it takes lots of arguments like the C version.
        ## Fortunately, someone (the python bindings?) are printing out lots of information
        ## when a shader fails.
        #info_log = glGetShaderInfoLog( self.fragment_shader )
        #print info_log
        
        glAttachShader( self.program, self.fragment_shader )
        glLinkProgram( self.program )
    
    def UnloadFragmentShader( self ):
        assert self.program is not None
        
        if self.fragment_shader is None:
            print("No fragment program loaded.")
            return
        
        glDetachShader( self.program, self.fragment_shader )
        glLinkProgram( self.program )
        
        glDeleteShader( self.fragment_shader )
        self.fragment_shader = None
    
    def UseProgram( self ):
        assert self.program is not None
        glUseProgram( self.program )
    
    def UnUseProgram( self ):
        assert self.program is not None
        glUseProgram( 0 )
    
    @staticmethod
    def shader_from_string( shader_type, string, shader = None ):
        if shader is None:
            shader = glCreateShader( shader_type )
        
        glShaderSource( shader, string )
        glCompileShader( shader )
        
        return shader
