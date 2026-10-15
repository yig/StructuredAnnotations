from PIL import Image
from OpenGL.GL import *

'''
A class to load and display images.
'''
class GLTexture:
    
    def __init__( self, filename ):
        self.id = 0
        self.width = 0
        self.height = 0
        self.actual_width = 1.0
        self.actual_height = 1.0
        
        self.LoadImage( Image.open( filename ) )
    
    def Unload( self ):
        if( glIsTexture( self.id ) ):
            assert( 0 != self.id )
            assert( self.width > 0 and self.height > 0 )
            assert( self.actual_width > 0.0 and self.actual_height > 0.0 )
            assert( self.actual_width <= 1.0 and self.actual_height <= 1.0 )
            glDeleteTextures( 1, self.id )
        
        self.id = 0
        self.width = 0
        self.height = 0
        self.actual_width = 1.0
        self.actual_height = 1.0
    
    def Load( self, filename ):
        self.LoadImage( Image.open( filename ) )
    
    def LoadImage( self, image ):
        '''
        Given a PIL Image object, convert it to a texture usable with OpenGL.
        '''
        if non_power_of_two_textures_allowed():
            self.width, self.height = image.size
            
            self.actual_width = 1.
            self.actual_height = 1.
        else:
            newsize = ( next_power_of_two( image.size[0] ), next_power_of_two( image.size[1] ) )
            oldsize = image.size
            
            self.width, self.height = image.size
            self.actual_width = float( oldsize[0] ) / float( newsize[0] )
            self.actual_height = float( oldsize[1] ) / float( newsize[1] )
            
            newimg = Image.new( "RGBA", newsize, (1,1,1,0) )
            newimg.paste( image, (0,newsize[1] - oldsize[1] - 1) )
            image = newimg
        
        image = image.convert( "RGBA" )
        rawimg = image.tobytes( "raw", "RGBA", 0, -1 )
        
        self.id = glGenTextures( 1 )
        glBindTexture( GL_TEXTURE_2D, self.id )
        
        glPixelStorei( GL_UNPACK_ALIGNMENT, 1 )
        glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, image.size[0], image.size[1], 0, GL_RGBA, GL_UNSIGNED_BYTE, rawimg )
        
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE )
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE )
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR )
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR )
        
        ## I can't set this here, because it's not attached to the texture object,
        ## but these are good values:
        ## multiply the texture by the color (tint)
        glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE )
        
        ## if we just want the texture unmodified by the current glColor
        # glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE ) 
    
    def Draw( self, pt, size, opacity = 1.0, upside_down = False, tint = None ):
        from glhelpers import draw_gltexture
        
        if tint is not None:
            glColor4f( tint[0], tint[1], tint[2], opacity )
        else: glColor4f( 1.0, 1.0, 1.0, opacity )
        draw_gltexture( self, list(pt), size, upside_down )

def next_power_of_two( n ):
    x = 1
    while x < n:
        x = x * 2
    return x

def non_power_of_two_textures_allowed():
    #print glGetString( GL_VERSION )
    #print glGetString( GL_EXTENSIONS )
    #print "GL_ARB_texture_non_power_of_two:", "GL_ARB_texture_non_power_of_two" in glGetString( GL_EXTENSIONS ).split(" ")
    
    if int( glGetString( GL_VERSION ).decode("utf-8", "strict").split(".")[0] ) >= 2: return True
    if "GL_ARB_texture_non_power_of_two" in glGetString( GL_EXTENSIONS ).decode("utf-8", "strict").split(" "): return True
    
    #print "non-power-of-two textures not supported."
    return False
