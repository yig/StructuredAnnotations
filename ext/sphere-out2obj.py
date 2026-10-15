#!/usr/bin/env python

import os, sys

if len( sys.argv ) < 2 or len( sys.argv ) > 3:
    print('Usage:', sys.argv[0], 'input.sphere.out', '[output.obj]')
    sys.exit(1)

in_fname = sys.argv[1]
out_fname = ''
if len( sys.argv ) == 3:
    out_fname = sys.argv[2]
else:
    out_fname = os.path.basename( in_fname )
    if in_fname == '-': out_fname = 'stdin'
    out_fname = out_fname + '.obj'

if os.path.exists( out_fname ):
    print('[Error: output file \"' + out_fname + '\" already exists]', file=sys.stderr)
    sys.exit( -1 )

print('[Writing output OBJ to \"' + out_fname + '\"]', file=sys.stderr)

if in_fname == '-':
    in_file = sys.stdin
else:
    in_file = file( in_fname )

fs = []
vs = []
for line in in_file:
    sline = line.split()
    if len( sline ) == 0: continue
    if sline[0] == 'triangle':
        fs.append( [] )
        continue
    else:
        ## Turn the vertex string into a tuple
        v = tuple( sline )
        ## Attempt to find the index of the vertex in the list.
        try:
            vindex = vs.index( v )
        ## If it's not, add it to the end (and take the index of the end).
        except ValueError:
            vindex = len( vs )
            vs.append( v )
        ## Add the vertex index to the face.
        fs[-1].append( vindex )

in_file.close()
del in_file

## Write out OBJ format.
if out_fname == '-':
    out_file = sys.stdout
else:
    out_file = file( out_fname, 'w' )
out_file.write( '# %s %s %s\n\n' % (sys.argv[0],in_fname,out_fname) )
## Write out vertices
for v in vs: out_file.write( 'v %s %s %s\n' % v )
out_file.write( '\n\n' )
## Write out faces
for f in fs: out_file.write( 'f %s %s %s\n' % tuple( [ vindex + 1 for vindex in f ] ) )

out_file.close()
del out_file
