## Following http://wiki.cython.org/WrappingCPlusPlus
## To build:
# python GLTexture-build.py build_ext --inplace

cdef extern from "GLTexture.hpp":
    ctypedef struct c_GLTexture "GLTexture":
        unsigned int id
        int width
        int height
        float actual_width
        float actual_height
        
        void Load( char* filename ) except +IOError
        void Unload()
        
    c_GLTexture *new_GLTexture "new GLTexture" ( char* filename ) except +IOError
    void del_GLTexture "delete" (c_GLTexture *gltex)

cdef class GLTexture:
    ## hold a C++ instance which we're wrapping
    cdef c_GLTexture *thisptr
    def __cinit__( self, char* filename ):
        self.thisptr = new_GLTexture( filename )
    def __dealloc__( self ):
        del_GLTexture( self.thisptr )
    
    def Load( self, char* filename ):
        self.thisptr.Load( filename )
    def Unload( self ):
        self.thisptr.Unload()
    
    property id:
        def __get__(self): return self.thisptr.id
    property width:
        def __get__(self): return self.thisptr.width
    property height:
        def __get__(self): return self.thisptr.height
    property actual_width:
        def __get__(self): return self.thisptr.actual_width
    property actual_height:
        def __get__(self): return self.thisptr.actual_height
