import direct.directbase.DirectStart
from pandac.PandaModules import *

def create_GeomVertexData( vs, name = '' ):
    '''
    Returns a Panda3D GeomVertexData given a list of 3d points 'vs'.
    This is useful to pass to other function such as
    create_Linestrip_with_GeomVertexData()
    and
    create_Triangles_with_GeomVertexData().
    '''
    
    ## Make the vertex data object and fill it with our vertices.
    format = GeomVertexFormat.getV3()
    vdata = GeomVertexData( name, format, Geom.UHStatic )
    vertex = GeomVertexWriter( vdata, 'vertex' )
    for v in vs: vertex.addData3f( *v )
    
    return vdata

def create_GeomTriangles( triangles ):
    '''
    Returns a Panda3D GeomTriangles data representing the Triangles 'triangles',
    a list of CCW 0-indexed triplets.
    '''
    
    ## Make the geometry primitive object and fill it with our triangles.
    prim = GeomTriangles( Geom.UHStatic )
    for tri in triangles:
        prim.addVertices( *tri )
        prim.closePrimitive()
    
    return prim

def create_Linestrip( strip ):
    '''
    Returns a Panda3D GeomLinestrips representing the line strip 'strip,
    a list of (0-indexed) indices.
    '''
    
    ## Make the geometry primitive object and fill it with our triangles.
    prim = GeomLinestrips( Geom.UHStatic )
    for index in strip: prim.addVertex( index )
    prim.closePrimitive()
    
    return prim

def create_Geom( vdata, primitives ):
    '''
    Returns a Panda3D Geom representing the GeomVertexData 'vdata'
    and the list of GeomPrimitives 'primitives'.
    '''
    
    ## Make the Geom object
    geom = Geom( vdata )
    for prim in primitives: geom.addPrimitive( prim )
    
    return geom

def create_GeomNode( geom, name = '' ):
    '''
    Returns a Panda3D GeomNode containing the given Geom object.
    '''
    
    ## Make the GeomNode
    node = GeomNode( name )
    node.addGeom( geom )
    
    return node

def add_GeomPrimitive_to_GeomNode( prim, node ):
    '''
    Adds the given GeomPrimitive 'prim' to the GeomNode 'node' so that they can share vertex data.
    '''
    assert node.getNumGeoms() == 1
    node.getGeom(0).addPrimitive( prim )

def create_Triangles_with_GeomVertexData( triangles, vdata, name = '' ):
    '''
    Returns a Panda3D GeomNode representing the Triangles 'triangles',
    a list of CCW 0-indexed triplets into 'vdata', the GeomVertexData.
    (To create a GeomVertexData call create_GeomVertexData() or use the convenience
    function create_Triangles_for_TriMesh().)
    '''
    
    prim = create_GeomTriangles( triangles )
    
    ## Make the Geom object
    geom = Geom( vdata )
    geom.addPrimitive( prim )
    
    ## Make the GeomNode
    node = GeomNode( name )
    node.addGeom( geom )
    
    return node

def create_Linestrip_with_GeomVertexData( strip, vdata, name = '' ):
    '''
    Returns a Panda3D GeomNode representing the line strip 'strip, a list of (0-indexed) indices
    into 'vdata', the GeomVertexData.
    (To create a GeomVertexData call create_GeomVertexData() or use the convenience
    function create_Linestrip_for_points().)
    '''
    
    ## Make the geometry primitive object and fill it with our triangles.
    prim = GeomLines( Geom.UHStatic )
    for index in strip: prim.addVertex( index )
    prim.closePrimitive()
    
    ## Make the Geom object
    geom = Geom( vdata )
    geom.addPrimitive( prim )
    
    ## Make the GeomNode
    node = GeomNode( name )
    node.addGeom( geom )
    
    return node

def create_TrianglesNode_for_TriMesh( mesh, name = '' ):
    '''
    Returns a Panda3D GeomNode given a TriMesh 'mesh', really anything with .vs and .faces fields where
    faces is a list of CCW 0-indexed triplets into .vs.
    '''
    
    vdata = create_GeomVertexData( mesh.vs, name )
    node = create_Triangles_with_GeomVertexData( mesh.faces, vdata, name )
    return node

def create_LinestripNode_for_points( points, name = '' ):
    '''
    Returns a Panda3D GeomNode given a list of 3d points.
    '''
    
    vdata = create_GeomVertexData( points, name )
    node = create_Linestrip_with_GeomVertexData( list(range(len(points))), vdata, name )
    return node

def merge_GeomNodes( node1, node2 ):
    '''
    Given two Panda3D GeomNodes, returns a new node merging them so that they share
    the same vertex data.
    This is convenient so you can create two Geoms and then merge them so that
    subsequent updates of the vertex data don't copy-on-write.
    '''
    
    ## I should assert this but I don't know the right call.
    # assert node1.getVertexData().number_of_rows == node2.getVertexData().number_of_rows
    
    assert node1.getNumGeoms() == 1
    assert node2.getNumGeoms() == 1
    
    geom1 = node1.getGeom(0)
    geom2 = node2.getGeom(0)
    
    ## Make the Geom object
    merged_geom = Geom( geom1.getVertexData() )
    for i in geom1.getNumPrimitives(): merged_geom.addPrimitive( geom1.getPrimitive(i) )
    for i in geom2.getNumPrimitives(): merged_geom.addPrimitive( geom2.getPrimitive(i) )
    
    ## Make the GeomNode
    node = GeomNode( node1.getName() )
    node.addGeom( merged_geom )
    
    return node

def update_GeomNode_vertices( node, vs ):
    '''
    Given a Panda3D GeomNode, for example created by create_Triangles_for_TriMesh(),
    update its vertices using the 'vs' parameter, obtainable as the field '.vs' of a TriMesh.
    '''
    update_GeomVertexData( node.modifyVertexData(), vs )

def update_GeomVertexData_vertices( vdata, vs, which = None ):
    '''
    Given a Panda3D GeomVertexData, for example created by create_GeomVertexData(),
    update its vertices using the 'vs' parameter, a list of (x,y,z) triplets with the same
    length as the GeomVertexData.
    If the optional 'which' parameter is not None, it is treated as a list of vertex indices to
    update.  The updated value for each vertex is obtained from the corresponding
    xyz triplet in 'vs'.  The length of 'which' and 'vs' must be the same.
    '''
    
    vertex = GeomVertexWriter( vdata, 'vertex' )
    if which is None:
        #while not vertex.isAtEnd():
        for v in vs:
            vertex.setData3f( *v )
        assert vertex.isAtEnd()
    else:
        for vi, v in zip( which, vs ):
            vertex.setRow( vi )
            vertex.setData3f( *v )
