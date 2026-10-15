from distutils.core import setup
from distutils.extension import Extension
try:
    from Cython.Distutils import build_ext
except ImportError:
    try:
        from Pyrex.Distutils import build_ext
    except:
        raise

ext = Extension(
    "GLTexture",                 # name of extension
    ["GLTexture.pyx"],           # filename of our Pyrex/Cython source
    language="c++",              # this causes Pyrex/Cython to create C++ source
    include_dirs=[
        ".",
        "/Users/yotam/Work/meshopt/git-devel/meshopt/src/new_viewer/yig"
        ],
    library_dirs=["/usr/X11/lib","/Users/yotam/Work/meshopt/git-devel/meshopt/src/new_viewer/glew"],
    libraries=["png","glew"],
    # if needed
    extra_link_args=[
        "/Users/yotam/Work/meshopt/git-devel/meshopt/src/new_viewer/yig/GLTexture.o",
        "/Users/yotam/Work/meshopt/git-devel/meshopt/src/new_viewer/yig/YImage.o",
        "-framework", "OpenGL"
        ]
    )
setup(
    ext_modules = [ext],
    cmdclass = {'build_ext': build_ext}
    )
