# Structured Annotations for 2D-to-3D Modeling

`Structured Annotations for 2D-to-3D Modeling` is a program for creating 3D models from 2D artwork, or guide images.

This software implements the [2009 SIGGRAPH Asia paper of the same name by Yotam Gingold, Takeo Igarashi, and Denis Zorin](http://cs.gmu.edu/~ygingold/annotations/).

Yotam Gingold is the author of this software. Please email <ygingold@gmu.edu> or <yotam@yotamgingold.com> with feedback or comments.

Except the third-party software listed in the `LICENSE` file, all source code and graphics have been placed into the public domain (CC0). Please see `LICENSE` for details.


## Troubleshooting

`Structured Annotations for 2D-to-3D Modeling` has been tested on Mac OS X versions 10.8–10.14.


## Running from source

The software can be run from the included source distribution. Follow the instructions below under dependencies. Then run it via:

    python3 wx_annotate_collection.py

or

    PYTHONPATH=$(pipenv --venv) python3 wx_annotate_collection.py

(We have to use PYTHONPATH instead of `pipenv run` or `pipenv shell` because [wxPython needs Python to think it is in a Framework](https://github.com/pypa/pipenv/issues/1416).
Alternatively, we can install the blank virtual environment using `venv` and then use `pipenv` for dependencies:

    python3 -m venv .venv
    pipenv install

Then `pipenv run python wx_annotate_collection.py` works.

## Creating the double-clickable application

    pipenv install --dev py2app
    pipenv run python setup.py py2app

Then open it from the `dist` folder.
See `release.HOW` for instructions on packaging the source and binary distributions.


## Dependencies

There are some system library dependencies, some Python dependencies, and some C/C++ files to turn into shared libraries. See `INSTALL` for details.

### System libraries

    brew install ftgl eigen

### Python dependencies

To install the Python dependencies, either use the `Pipfile` by running `pipenv install` or use `pip3 install -r requirements.txt`.
To avoid the `PYTHONHOME` workaround mentioned above,
run `python3 -m venv .venv` before `pipenv install`.

### Compiled code

To compile the compiled code:

    cythonize -i geometry.pyx
    # If the cythonize command isn't in the path, it is accessible like so:
    # python -m Cython.Build.Cythonize -i geometry.pyx

    cc tritri.c tritri_wrapper.cpp raytri.c raytri_wrapper.cpp -lstdc++ -fkeep-inline-functions -dynamiclib -o libtritri.dylib -g -O3 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare

    c++ Eline_ad.mitsuba.cpp -dynamiclib -o Eline_ad.dylib -g -O3 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare -Iext -Iext/mitsuba -I/usr/local/include/eigen3

    c++ edge_distances.cpp -dynamiclib -o edge_distances.dylib -g -O3 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare

    cc face_stitch_pattern.c -dynamiclib -o face_stitch_pattern.dylib -std=c99 -g -O3 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare

    c++ laplacian_editing_wrapper.cpp -Iext -dynamiclib -o liblaplacian_editing.dylib -g -O1 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare
