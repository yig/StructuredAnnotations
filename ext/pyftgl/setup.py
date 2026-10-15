#!/usr/bin/env python

from distutils.core import setup, Extension

__version__ = "0.5b"
__author__ = "Anders Dahnielson"
__doc__ = """Python FTGL binding

FTGL is a free, open source library to enable developers to use
arbitrary fonts in their OpenGL applications. Unlike other OpenGL font
libraries FTGL uses standard font file formats and the Freetype font
library to open the fonts. PyFTGL in turn wraps the functionality of
FTGL into a Python module so that it can be used in conjunction with
PyOpenGL."""

doclines = __doc__.split("\n")

##########################################################################
## FTGL module

module_ftgl_src = [
    'ftgl.cpp',
    ]

module_ftgl_include_dirs = [
    '/usr/include',
    '/usr/local/include',
    '/usr/include/freetype2',
    '/usr/local/include/freetype2',
    '/opt/homebrew/include/',
    '/opt/homebrew/include/freetype2',
    ]

import sys
module_ftgl_libs = [
    # 'GLU',
    # 'GL',
    'freetype',
    'z',
    'ftgl',
    'boost_python' + ( '%s%s' % ( sys.version_info[0], sys.version_info[1] ) ),
    ]

module_ftgl = Extension(
    'FTGL',
    module_ftgl_src,
    include_dirs=module_ftgl_include_dirs,
    libraries=module_ftgl_libs,
    library_dirs = ['/usr/local/lib','/opt/homebrew/lib'],
    extra_link_args=['-framework', 'OpenGL']
    )

##########################################################################
## setup

setup(
    name = 'PyFTGL',
    version = __version__,
    author = __author__,
    author_email = 'anders@dahnielson.com',
    maintainer = __author__,
    maintainer_email = 'anders@dahnielson.com',
    url = 'http://code.google.com/p/pyftgl/',
    description = doclines[0],
    long_description = "\n".join(doclines[2:]),
    download_url = 'http://code.google.com/p/pyftgl/downloads/list',
    classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: GNU General Public License (GPL)',
    'Operating System :: POSIX',
    'Programming Language :: C++',
    'Programming Language :: Python',
    'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    license = "GPLv2",
    platforms = "POSIX",
    ext_modules = [module_ftgl],
    )
