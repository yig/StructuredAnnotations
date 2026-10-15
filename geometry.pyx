## With cythonize
# cythonize -i geometry.pyx
## If the cythonize command isn't in the path, it is accessible like so:
# python -m Cython.Build.Cythonize -i geometry.pyx

## With my built-in Leopard python with Cython installed:
## python geometry-build.py build_ext --inplace

## pyrexc geometry.pyx
## gcc -O3 -c -fPIC geometry.c -I/System/Library/Frameworks/Python.framework/Headers
## gcc -dynamiclib -o geometry.so geometry.o /System/Library/Frameworks/Python.framework/Python

## cp geometry.pyx geometry25.pyx
## python2.5 geometry-build.py build_ext --inplace

from __future__ import print_function

cdef extern from "math.h":
    float fabs( float )

cdef struct pair:
    float one
    float two

#import animate_path
cdef float eps
#eps = animate_path.kEps
eps = 1e-5

## Erm, I wish I could declare this as returning a tuple of floats: cdef (float, float)
def line_line_intersection( line1, line2 ):
    cdef pair p
    p = __line_line_intersection( line1, line2 )
    return p.one, p.two
cdef pair __line_line_intersection( line1, line2 ) except *:
    ## '''
    ## Takes two 2d lines 'line1' = (start1, end1) and 'line2' = (start2, end2)
    ## and returns t1, t2 such that lerp( start1, end1, t1 ) == lerp( start2, end2, t2 ).
    ## 
    ## tested
    ## '''
    
    ## from http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
    
    cdef pair p
    cdef float x1, x2, x3, x4, y1, y2, y3, y4, t_denom, t1, t2
    #eps = animate_path.kEps
    
    #print 'line1[0]:', line1[0]
    
    ## Unpack line coordinates
    #x1, y1 = line1[0]
    #x2, y2 = line1[1]
    #x3, y3 = line2[0]
    #x4, y4 = line2[1]
    
    x1 = line1[0][0]
    y1 = line1[0][1]
    x2 = line1[1][0]
    y2 = line1[1][1]
    x3 = line2[0][0]
    y3 = line2[0][1]
    x4 = line2[1][0]
    y4 = line2[1][1]
    
    #print 'line_segments_cross( %s, %s )' % (str(line1), str(line2))
    #print 'line_segments_cross (x1,y1), (x2,y2), (x3,y3), (x4,y4):', (x1,y1), (x2,y2), (x3,y3), (x4,y4)
    
    t_denom = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
    
    t1 = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)
    t2 = (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)
    
    if fabs( t_denom ) <= eps:
        ## line are coincident
        if fabs( t1 ) <= eps and fabs( t2 ) <= eps:
            #set_trace()
            ## shouldn't happen!
            #assert False and "lines are coincident"
            raise RuntimeError, "lines are coincident"
            p.one = .5
            p.two = .5
            return p
        
        ## lines are parallel
        else:
            p.one = -1.
            p.two = -1.
            return p
    
    t1 = t1 / t_denom
    t2 = t2 / t_denom
    
    p.one = t1
    p.two = t2
    return p

def line_segments_cross( line1, line2 ):
    ## '''
    ## Takes two 2d line segments 'line1' = (start1, end1) and 'line2' = (start2, end2)
    ## and returns True if the line segments cross, false otherwise.
    ## 
    ## tested
    ## '''
    
    cdef pair p
    cdef float t1, t2
    
    p = __line_line_intersection( line1, line2 )
    t1 = p.one
    t2 = p.two
    
    if(
        t1 > -eps and
        t1 < 1 + eps and
        t2 > -eps and
        t2 < 1 + eps
        ):
        return True
    else:
        return False

## Erm, I wish I could declare this as returning a tuple of floats: cdef (float, float)
def line_line_closest_points2d( line1, line2 ):
    assert len( line1[0] ) == 2
    
    line1 = ( tuple( line1[0] ) + (0,), tuple( line1[1] ) + (0,) )
    line2 = ( tuple( line2[0] ) + (0,), tuple( line2[1] ) + (0,) )
    return line_line_closest_points3d( line1, line2 )

## Erm, I wish I could declare this as returning a tuple of floats: cdef (float, float)
def line_line_closest_points3d( line1, line2 ):
    assert len( line1[0] ) == 3
    
    cdef pair p
    p = __line_line_closest_points3d( line1, line2 )
    return p.one, p.two

cdef pair __line_line_closest_points3d( line1, line2 ) except *:
    ## '''
    ## Takes two 3d lines 'line1' = (start1, end1) and 'line2' = (start2, end2)
    ## and returns t1, t2 such that lerp( start1, end1, t1 ) and lerp( start2, end2, t2 )
    ## are the closest points to each other on the two line segments.
    ## 
    ## untested
    ## '''
    
    ## from http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline3d/
    
    cdef pair p
    cdef float p1_x, p2_x, p3_x, p4_x, p1_y, p2_y, p3_y, p4_y, p1_z, p2_z, p3_z, p4_z
    
    cdef float p13_x, p13_y, p13_z, p43_x, p43_y, p43_z, p21_x, p21_y, p21_z
    cdef float d1343, d4321, d1321, d4343, d2121
    cdef float numer, denom
    cdef float mua, mub
    
    #print 'line1[0]:', line1[0]
    
    ## Unpack line coordinates
    #p1_x, p1_y, p1_z = line1[0]
    #p2_x, p2_y, p2_z = line1[1]
    #p3_x, p3_y, p3_z = line2[0]
    #p4_x, p4_y, p4_z = line2[1]
    
    p1_x = line1[0][0]
    p1_y = line1[0][1]
    p1_z = line1[0][2]
    p2_x = line1[1][0]
    p2_y = line1[1][1]
    p2_z = line1[1][2]
    p3_x = line2[0][0]
    p3_y = line2[0][1]
    p3_z = line2[0][2]
    p4_x = line2[1][0]
    p4_y = line2[1][1]
    p4_z = line2[1][2]
    
    p13_x = p1_x - p3_x
    p13_y = p1_y - p3_y
    p13_z = p1_z - p3_z
    p43_x = p4_x - p3_x
    p43_y = p4_y - p3_y
    p43_z = p4_z - p3_z
    if (fabs(p43_x)  < eps and fabs(p43_y)  < eps and fabs(p43_z)  < eps):
        raise RuntimeError, "both line points are the same"
    p21_x = p2_x - p1_x
    p21_y = p2_y - p1_y
    p21_z = p2_z - p1_z
    if (fabs(p21_x)  < eps and fabs(p21_y)  < eps and fabs(p21_z)  < eps):
        raise RuntimeError, "both line points are the same"
    
    d1343 = p13_x * p43_x + p13_y * p43_y + p13_z * p43_z
    d4321 = p43_x * p21_x + p43_y * p21_y + p43_z * p21_z
    d1321 = p13_x * p21_x + p13_y * p21_y + p13_z * p21_z
    d4343 = p43_x * p43_x + p43_y * p43_y + p43_z * p43_z
    d2121 = p21_x * p21_x + p21_y * p21_y + p21_z * p21_z
    
    denom = d2121 * d4343 - d4321 * d4321
    if (fabs(denom) < eps):
        raise RuntimeError, "lines are parallel"
    numer = d1343 * d4321 - d1321 * d4343
    
    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343
    
    p.one = mua
    p.two = mub
    return p

def test_timing():
    import timeit
    setup = setup = '''
import geometry
l1 = ((0,0), (1,0))
l2 = ((.5,-1), (.5,1))

z1 = ((0,0,0), (1,0,0))
z2 = ((.5,-1,0), (.5,1,0))
'''
    
    N = 1000000
    print( "timing test iterations:", N )
    print( "line_line_intersection():", timeit.Timer( 'geometry.line_line_intersection( l1, l2 )', setup ).timeit( N ) )
    print( "line_segments_cross():", timeit.Timer( 'geometry.line_segments_cross( l1, l2 )', setup ).timeit( N ) )
    print( "line_line_closest_points3d():", timeit.Timer( 'geometry.line_line_closest_points3d( z1, z2 )', setup ).timeit( N ) )
    print( "line_line_closest_points2d():", timeit.Timer( 'geometry.line_line_closest_points2d( l1, l2 )', setup ).timeit( N ) )
