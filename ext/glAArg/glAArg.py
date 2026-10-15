try:
    import ctypes
except:
    print('''
ERROR: ctypes not installed properly. (ctypes needed for laplacian_editing_1d())
        ''')
    import sys
    sys.exit()

from OpenGL import GL as _GL

def platform_shared_library_suffix():
    import sys
    result = 'so'
    if 'win' in sys.platform.lower(): result = 'dll'
    ## No else if, because we want darwin to override win (which is a substring of darwin)
    if 'darwin' in sys.platform.lower(): result = 'dylib'
    return result

import os.path
lib = ctypes.cdll.LoadLibrary( os.path.join( '.', 'libglAArg.' + platform_shared_library_suffix() ) )

phf = 32
glAA_texture = ( _GL.GLuint * 8 ).in_dll( lib, 'glAA_texture' )
GLAA_VERTEX_ARRAY = 0x851D

class GLCoord2( ctypes.Structure ):
    _fields_ = [
        ( "x", _GL.GLfloat ),
        ( "y", _GL.GLfloat )
        ]

## Use an init function to prevent my locals from leaking into the module.
def cimport():
    ## meta-programming to create entries for every function in glAArg
    types2names = [
        ( None, ( 'glAAInit', 'glAAExit', 'glAAEnd', 'glAAFlush' ) ),
        ( ( _GL.GLenum, ), ( 'glAABegin', 'glAAEnable', 'glAADisable' ) ),
        ( ( ctypes.c_float, _GL.GLuint, ctypes.c_float ), ( 'glAAGenerateAATex', ) ),
        ( ( _GL.GLfloat, ), ( 'glAAPointSize', 'glAALineWidth' ) ),
        ( ( GLCoord2, ), ( 'glAAPointLight', ) ),
        ( ( _GL.GLint, _GL.GLushort ), ( 'glAALineStipple', ) ),
        ( ( _GL.GLfloat, _GL.GLushort, _GL.GLfloat ), ( 'glAALineStipplePhase', ) ),
        ( ( _GL.GLfloat, _GL.GLfloat, _GL.GLfloat ), ( 'glAAColor3f', ) ),
        ( ( _GL.GLfloat, _GL.GLfloat, _GL.GLfloat, _GL.GLfloat ), ( 'glAAColor4f', ) ),
        ( ( ctypes.POINTER( _GL.GLubyte ), ), ( 'glAAColor4ubv', ) ),
        ( ( _GL.GLuint, ), ( 'glAAColor1ui', ) ),
        ( ( _GL.GLfloat, _GL.GLfloat ), ( 'glAAVertex2f', ) )
        ]
    for argtypes, funcnames in types2names:
        for funcname in funcnames:
            func = getattr( lib, funcname )
            func.restype = None
            func.argtypes = argtypes
            globals()[ funcname ] = func
cimport()
del cimport

_glAAColor4ubv = glAAColor4ubv
def glAAColor4ubv( four_unsigned_bytes ):
    four_unsigned_bytes_array = ( _GL.GLubyte * 4 )()
    four_unsigned_bytes_array[:] = four_unsigned_bytes
    _glAAColor4ubv( four_unsigned_bytes_array )


## Helper functions.
def glAAPrepareBlendingState():
    '''
    Helper function to set up the required blending and alpha state, and turns off lighting.
    '''
    _GL.glDisable( _GL.GL_LIGHTING )
    _GL.glDisable( _GL.GL_DITHER )
    _GL.glShadeModel( _GL.GL_SMOOTH )
    _GL.glEnable( _GL.GL_BLEND )
    _GL.glBlendFunc( _GL.GL_SRC_ALPHA, _GL.GL_ONE_MINUS_SRC_ALPHA )
    _GL.glEnable( _GL.GL_ALPHA_TEST )
    _GL.glAlphaFunc( _GL.GL_GREATER, 1/255.0 )

__texture_2D_was_enabled = None
def glAAPrepareTextureState():
    '''
    Helper function to set up the required texture state.
    
    NOTE: This does not bind to one of the glAA_texture[]'s.
    NOTE: This function affects the current matrix mode and the texture matrix.
    '''
    global __texture_2D_was_enabled
    __texture_2D_was_enabled = _GL.glIsEnabled( _GL.GL_TEXTURE_2D )
    
    #_GL.glPixelStorei( _GL.GL_UNPACK_ALIGNMENT, 1 )
    _GL.glMatrixMode( _GL.GL_TEXTURE )
    _GL.glLoadIdentity()
    _GL.glScalef( 0.25/phf, 0.25/phf, 1.0 )
    _GL.glEnable( _GL.GL_TEXTURE_2D )

def glAARestorePriorTextureState():
    '''
    Helper function which resets the texture state and restores the prior
    value of glIsEnabled( GL_TEXTURE_2D ).
    
    NOTE: This function affects the current matrix mode and the texture matrix.
    '''
    global __texture_2D_was_enabled
    assert __texture_2D_was_enabled is not None
    if not __texture_2D_was_enabled: _GL.glDisable( _GL.GL_TEXTURE_2D )
    __texture_2D_was_enabled = None
    
    _GL.glMatrixMode( _GL.GL_TEXTURE )
    _GL.glLoadIdentity()
