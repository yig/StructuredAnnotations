from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *
from helpers import *
from vecutil import *
kEps = 1e-5

from OpenGL.GL import *
from glhelpers import *
from GLTexture import GLTexture

def draw_widget( widget_name, pos, yvec, size ):
    '''
    Given a widget's name as a string,
    a 2+d center position 'pos',
    a 2d vector orienting the widget by specifying its image's y direction, and
    the edge length of its square,
    draws the widget.
    '''
    
    draw_widget_tex( get_texture(widget_name), pos, yvec, size )

def draw_widget_tex( tex, pos, yvec, size ):
    '''
    Given a widget's GLTexture object 'tex',
    a 2+d center position 'pos',
    a 2d vector orienting the widget by specifying its image's y direction, and
    the edge length of its square,
    draws the widget.
    '''
    
    pos3d = zerosf(3)
    pos3d[ :len(pos) ] = pos
    
    ## We want a 3d vector, and we only want the 2d y-vector direction to affect the widget.
    yvec3d = append( yvec[:2], 0. )
    
    tex_was_enabled = glIsEnabled( GL_TEXTURE_2D )
    
    #glColor3f( 1,1,1 )
    glEnable( GL_TEXTURE_2D )
    
    glMatrixMode( GL_MODELVIEW )
    glPushMatrix()
    glTranslate( *pos3d )
    axis, radians = get_axis_angle_rotation( (0,1,0), yvec3d )
    glRotate( (180./pi)*radians, *axis )
    glScale( size, size, 1. )
    
    dims = arrayf( ( tex.width, tex.height ) )
    dims /= max( dims )
    draw_gltexture( tex, -.5*dims, dims )
    
    glPopMatrix()
    
    if not tex_was_enabled: glDisable( GL_TEXTURE_2D )

_textures = {}
def get_texture( tex_name ):
    if tex_name not in _textures:
        _textures[ tex_name ] = GLTexture( os.path.join( '.', 'widgets', tex_name ) + '.png' )
    
    return _textures[ tex_name ]
