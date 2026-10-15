'''
I used values of npts (1,10,100,1000, and even 10000 which took almost 10 seconds for each run).

I tried with vecutil.axes_inverse() set to use linalg.inv and the custom code, which
made a factor of 2 difference for testing 1 or 10 points inside the ellipse,
a 50% improvement for 100 points, and a 20% improvement for 1000 points.

Passing arrays makes a 50% speed difference, no matter the number of points (meaning the
work necessary to create the array is proportional to the number of points, just as the algorithm
is, and coefficient of proportion is 1.
'''

npts = 10

setup2d = '''
from vecutil import ellipse_contains_points
ellipse = ( (0,0), [(1,0), (0,1)] )
points = [(0,0)] * %s
''' % (npts,)

setup2darr = '''
from vecutil import ellipse_contains_points
ellipse = ( (0,0), [(1,0), (0,1)] )
points = [(0,0)] * %s

from myarray import asarrayf
ellipse = ( ellipse[0], asarrayf( ellipse[1] ) )
points = asarrayf( points )
''' % (npts,)

setup3d = '''
from vecutil import ellipse_contains_points
ellipse = ( (0,0,0), [(1,0,0), (0,1,0), (0,0,1)] )
points = [(0,0,0)] * %s
''' % (npts,)
setup3darr = '''
from vecutil import ellipse_contains_points
ellipse = ( (0,0,0), [(1,0,0), (0,1,0), (0,0,1)] )
points = [(0,0,0)] * %s

from myarray import asarrayf
ellipse = ( ellipse[0], asarrayf( ellipse[1] ) )
points = asarrayf( points )
''' % (npts,)

import timeit
N = 1000

print('num points:', npts)
print('num repititions:', N)

print('2d:', end=' ')
print(timeit.Timer( 'ellipse_contains_points( ellipse, points )', setup2d ).timeit(N))
print('3d:', end=' ')
print(timeit.Timer( 'ellipse_contains_points( ellipse, points )', setup3d ).timeit(N))

print('2d array:', end=' ')
print(timeit.Timer( 'ellipse_contains_points( ellipse, points )', setup2darr ).timeit(N))
print('3d array:', end=' ')
print(timeit.Timer( 'ellipse_contains_points( ellipse, points )', setup3darr ).timeit(N))
