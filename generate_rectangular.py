#!/usr/bin/env python

# This script generates either a regular, a 4-8, or an equilateral
# flat rectangular mesh and saves it as an OBJ.  This script also
# allows a user to specify the size of a boundary, functions for
# initial positions of boundary and nonboundary vertices, and a fixed
# dof pattern for vertices not on the boundary.
# [jasonr_at_nyu.edu]

import sys, os, getopt
import math
import random

from math import *

# this constant is useful for switching which rows 
# in equilateral meshes contain the extra vertex
ODD_EVEN = 1

# generates points for regular and 4-8 rectangular meshes
# u_function is the undeformed value function for non-boundary vertices
# i_function is the initial value function for non-boundary vertices
# b_function is the initial value function for boundary vertices
# b is the size of the fixed boundary
def generate_points( x0, y0, width, height, num_rows, num_cols, u_function, i_function, b_function, b=0 ):
    undeformed_points = []
    initial_points = []    
    x_step = float( width ) / num_cols
    y_step = float( height ) / num_rows
    points_this_row = num_cols + 1
    y = y0
    for i in range( 0, num_rows + 1 ):
        x = x0
        for j in range( 0, points_this_row ):
            undeformed_points.append( u_function( x,y ) )

            if i < b or i > num_rows-b or j < b or j >= points_this_row-b:
                if b_function != None:
                    initial_points.append( b_function( x,y ) )
                else:
                    initial_points.append( (x,y,0) )
            else:
                initial_points.append( i_function( x,y ) )                
            x += x_step
        y += y_step
    return ( undeformed_points, initial_points )

# generates points for equilateral meshes
def generate_equilateral_points( x0, y0, side_length, num_rows, triangles_per_row, u_function, i_function, b_function, b=0 ):
    grid_points = []
    initial_points = []
    half_side_length = 0.5 * side_length    
    hypotensus = half_side_length / math.sin( math.radians( 30 ) )
    height = math.sqrt( hypotensus**2 - half_side_length**2 )
    y = y0
    for i in range( 0, num_rows + 1 ):
        if i % 2 == ODD_EVEN:
            points_per_row = triangles_per_row + 1
            x = x0
            for j in range( 0, points_per_row ):
                grid_points.append( u_function( x,y ) )
                if i < b or i > num_rows-b or j < b or j >= points_per_row-b:
                    initial_points.append( b_function( x,y ) )
                else:
                    initial_points.append( i_function( x,y ) )                    
                x += side_length
        else:
            points_per_row = triangles_per_row + 2
            x = x0+half_side_length            
            grid_points.append( u_function( x0,y ) )
            #if i < b or i > points_per_row-b:
            if b:
                initial_points.append( b_function( x0,y ) )
            else:
                initial_points.append( i_function( x0,y ) )                
            for k in range( 2, points_per_row ):
                grid_points.append( u_function( x,y ) )
                if i < b or i > num_rows-b or k < b+1 or k >= points_per_row-b+1:
                    initial_points.append( b_function( x,y ) )
                else:
                    initial_points.append( i_function( x,y ) )                                    
                x += side_length
            x -= half_side_length
            grid_points.append( u_function( x,y ) )
            #if i < b or i > points_per_row-b:
            if b:
                initial_points.append( b_function( x,y ) )
            else:
                initial_points.append( i_function( x,y ) )                          
        y += height
    return ( grid_points, initial_points )

# generates regular faces (to be used with regular grid points)
def generate_regular_triangles( num_rows, num_cols ):
    triangles = []
    cur = 0
    points_this_row = num_cols + 1
    triangles_this_row = 2 * num_cols
    for i in range( 0, num_rows ):
        for j in range( 0, triangles_this_row ):
            if j % 2 == 0:
                triangles.append( ( cur, cur+1, cur+points_this_row ) )
                cur += 1
            else:
                triangles.append( ( cur, cur+points_this_row, cur+points_this_row-1 ) )
        cur += 1
    return triangles

# generates 4-8 faces (to be used with regular grid points)
def generate_four_eight_triangles( num_rows, num_cols ):
    triangles = []
    cur = 0
    points_this_row = num_cols + 1
    triangles_this_row = 2 * num_cols
    count = 0
    down = True
    for i in range( 0, num_rows ):
        down = i % 2 == 0
        count = 0
        for j in range( 0, triangles_this_row ):
            if down:
                if j % 2 == 0:
                    triangles.append( ( cur, cur+1, cur+points_this_row ) )
                    cur += 1
                else:
                    triangles.append( ( cur, cur+points_this_row, cur+points_this_row-1 ) )
            else:
                if j % 2 == 0:
                    triangles.append( ( cur, cur+points_this_row+1, cur+points_this_row ) )
                else:
                    triangles.append( ( cur, cur+1, cur+points_this_row+1 ) )
                    cur += 1                    
            count += 1
            if count == 2:
                count = 0
                down = not down
        cur += 1
    return triangles

# generates triangles for use with the equilateral grid
def generate_equilateral_triangles( num_rows, triangles_per_row ):
    triangles = []
    cur = 0
    for i in range( 0, num_rows ):
        if i % 2 == ODD_EVEN:
            points_per_row = triangles_per_row + 1
        else:
            points_per_row = triangles_per_row + 2            
        for j in range( 0, 2*triangles_per_row+1 ):
            if i % 2 == ODD_EVEN:
                if j % 2 == 0:
                    triangles.append( (cur,cur+points_per_row,cur+points_per_row+1) )
                else:
                    triangles.append( (cur,cur+points_per_row+1,cur+1) )
                    cur += 1
            else:
                if j % 2 == 0:
                    triangles.append( (cur,cur+points_per_row,cur+1) )
                    cur += 1
                else:
                    triangles.append( (cur,cur+points_per_row-1,cur+points_per_row) )
        cur += 1
    return triangles

# fixes the boundary dofs for a regular grid
# b is the size of the boundary (e.g. 1, 2)
def fix_dofs( rows, cols, b, nonfixed_pattern ):
    fixed_dofs = []
    v = 1
    for r in range( 0, rows ):
        for c in range( 0, cols ):
            if r < b or r >= rows-b or c < b or c >= cols-b:
                fixed_dofs.append( "b " + str( v ) + " 1 1 1" )
            else:
                fixed_dofs.append( "b " + str( v ) + " " + nonfixed_pattern )                
            v += 1
    return fixed_dofs

# fixes the boundary dofs for an equilateral grid
# b is the size of the boundary (e.g. 1, 2)
def fix_dofs_equilateral( rows, cols, b, nonfixed_pattern ):
    fixed_dofs = []
    v = 1
    for r in range( 0, rows ):
        points_this_row = cols
        if r % 2 == 0:
            points_this_row +=1
        for c in range( 0, points_this_row ):
            if r < b or r >= rows-b or c < b or c >= points_this_row-b:
                fixed_dofs.append( "b " + str( v ) + " 1 1 1" )
            else:
                fixed_dofs.append( "b " + str( v ) + " " + nonfixed_pattern )
            v += 1
    return fixed_dofs

# writes an obj file given points and triangles (faces)
def write_obj( filename, options, points, initial, triangles, fixed ):
    f = open( filename, "w" )
    f.write( "# generated by generate_rectangular.py\n" )
    f.write( "# invoked with: ")
    for o in options:
        f.write( o + " " )
    f.write( "\n\n" )
    f.write( "# vertices\n" )
    for p in points:
        f.write( "v " + str( p[0] ) + " " + str( p[1] ) + " " + str( p[2] ) + "\n" )
    f.write( "\n# triangles\n" )    
    for t in triangles:
        f.write( "f " + str( t[0] + 1 ) + " " + str( t[1] + 1 ) + " " + str( t[2] + 1 ) + "\n" )
    f.write( "\n# initial positions\n" )
    k = 1
    for i in initial:
        f.write( "i " + str( k ) + " " + str( i[0] ) + " " + str( i[1] ) + " " + str( i[2] ) + "\n" )
        k += 1
    f.write( "\n# fixed dofs\n" )
    for d in fixed:
        f.write( d + "\n" )
    f.close()

def print_options():
    print("\nOptions:")
    print("\t--type\t\t\t the type of mesh to generation (regular, 4-8, or eq) [required]")    
    print("\t--width\t\t\t the width of the generated mesh [required]")
    print("\t--height\t\t the height of the generated mesh [required for regular and 4-8 but not equilateral meshes]")
    print("\t--rows\t\t\t the number of rows of vertices [required]")
    print("\t--cols\t\t\t the number of columns of vertices [required]")
    print("\t--x0\t\t\t x coordinate of lower left corner [optional, default 0]")
    print("\t--y0\t\t\t y coordinate of lower left corner [optional, default 0]")    
    print("\t--undeformed\t\t a function (of x and y) describing the undeformed positions of the non-boundary vertices[optional, default (x,y,0)]")
    print("\t--initial\t\t a function (of x and y) describing the initial positions of the non-boundary vertices [optional]")
    print("\t--boundary-size\t\t the depth of fixed vertices along the boundary (defaults to 0) [optional]")
    print("\t--boundary\t\t a function (of x and y) describing the initial positions of the boundary vertices [optional]")
    print("\t--nonboundary-pattern\t a string describing which dofs are fixed for non-boundary vertices [optional]")
    print("\t--filename\t\t path to file to save the resulting mesh in [required]")    

def print_examples():
    print("\nExamples:")
    print("\tgenerate_rectangular.py --type=regular --width=11.6 --height=7.0 --rows=5 --cols=3 --initial=\"(x,y,x+y)\" --filename=rect.obj")
    print("\tgenerate_rectangular.py --type=regular --width=2.0 --height=2.0 --rows=5 --cols=5 --initial=\"(x,y,x+y)\" --boundary-size=1 --nonboundary-pattern=\"0 1 0\" --filename=rect.obj")    
    print("\tgenerate_rectangular.py --type=4-8 --width=11.6 --height=7.0 --rows=5 --cols=4 --initial=\"(x,y,x*y)\" --filename=rect.obj")
    print("\tgenerate_rectangular.py --type=eq --width=10.0 --rows=5 --cols=6 --initial=\"(x,y,x+y)\" --boundary-size=1 --filename=rect.obj")

def usage():
    print("generate_rectangular.py")
    print("\nUsage: generate_rectangular.py --type=<regular|4-8|eq> --width=<width> --height=<height> --rows=<num rows> --cols=<num columns> [--undeformed=<f(x,y)->(x,y,z)>] [--initial=<f(x,y)->(x,y,z)>] [--boundary-size=<boundary size>] [--nonboundary-pattern=<1|0 1|0 1|0>] [--boundary=<f(x,y)->(x,y,z)] --filename=<obj output filename>")
    print_options()
    print_examples()
    sys.exit( -1 )

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt( sys.argv[1:], "", ["type=", "width=", "height=", "rows=", "cols=", "undeformed=", "initial=", "filename=", "boundary=", "nonboundary-pattern=", "boundary-size=", "x0=", "y0="] )
    except getopt.GetoptError:
        usage()

    # default values
    type = "regular"
    x0 = 0.0
    y0 = 0.0    
    width = 1.0 # width of mesh
    height = 1.0 # height of mesh
    rows = 1 # number of rows
    cols = 1 # number of columns
    filename = None # output filename
    u_function = lambda x,y : (x,y,0)  # function for undeformed position (defaults to flat)
    i_function = lambda x,y : (x,y,0)  # function for initial position (defaults to flat)
    b_function = lambda x,y : (x,y,0)  # boundary function (like i_function but for fixed boundary)
    fixed = 0  # number of fixed points on boundary
    nonfixed_pattern="0 0 0" # the fixed dof pattern for non-boundary vertices


    # read options
    optsMap = {}
    for n, v in opts:
        optsMap[n[2:]] = v

    if "type" in optsMap:
        type = optsMap["type"]
    if "x0" in optsMap:
        x0 = float( optsMap["x0"] )
    if "y0" in optsMap:
        y0 = float( optsMap["y0"] )                
    if "width" in optsMap:
        width = float( optsMap["width"] )
    if "height" in optsMap:
        height = float( optsMap["height"] )
    if "rows" in optsMap:
        rows = int( optsMap["rows"] )
    if "cols" in optsMap:
        cols = int( optsMap["cols"] )
    if "undeformed"  in optsMap:
        u_function = eval( "lambda x,y : ( " + optsMap["undeformed"] + " )" )        
    if "initial"  in optsMap:
        i_function = eval( "lambda x,y : ( " + optsMap["initial"] + " )" )        
    if "boundary" in optsMap:
        b_function = eval( "lambda x,y : ( " + optsMap["boundary"] + " )" )
    if "filename" in optsMap:
        filename = optsMap["filename"]
    if "boundary-size" in optsMap:
        fixed = int( optsMap["boundary-size"] )
    if "nonboundary-pattern" in optsMap:
        nonfixed_pattern = optsMap["nonboundary-pattern"]
        print(nonfixed_pattern)
        
    # make sure we have a real type
    if type != "regular" and type != "4-8" and type != "eq":
        usage()

    # make sure we have an output filename
    if filename == None or (i_function == None and u_function == None):
        usage()        

    fixed_dofs = []
    if type == "regular" or type == "4-8":
        (grid_points,initial_points) = generate_points( x0, y0, width, height, rows-1, cols-1, u_function, i_function, b_function, fixed )
        if type == "regular":
            triangles = generate_regular_triangles( rows-1, cols-1 )
        else:
            triangles = generate_four_eight_triangles( rows-1, cols-1 )
        if fixed:
            fixed_dofs = fix_dofs( rows, cols, fixed, nonfixed_pattern )
    elif type == "eq":
        side_length = width/(cols-1)
        (grid_points,initial_points) = generate_equilateral_points( x0, y0, side_length, rows-1, cols-1, u_function, i_function, b_function, fixed )
        triangles = generate_equilateral_triangles( rows-1, cols-1 )
        if fixed:
            fixed_dofs = fix_dofs_equilateral( rows, cols, fixed, nonfixed_pattern )
    write_obj( filename, sys.argv[1:], grid_points, initial_points, triangles, fixed_dofs  )
