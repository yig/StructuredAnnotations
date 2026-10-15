from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *

from OpenGL.GL import *
from OpenGL.GLU import *
## GLUT launches X11
# from OpenGL.GLUT import *

class Ellipsoid( object ):
    def __init__( self ):
        self.center = zerosf( 3 )
        ## These should be orthogonal.  The magnitudes determine the shape of the sphere.
        self.axes = arrayf( [ (1,0,0), (0,1,0), (0,0,1) ] )
    
    def RayIntersection( self, ray, ray_is_line = False ):
        '''
        Given a ray ( 3d point, 3d direction ), returns the closest intersection
        point with the ellipsoid along 'ray', or None if they don't intersect.
        If 'ray_is_line' is True, then the closest point may lie along the negative ray direction.
        '''
        
        from vecutil import closest_point_to_ellipse_on_ray
        closest_pt, intersected = closest_point_to_ellipse_on_ray(
            ( self.center, self.axes ),
            ray,
            ray_is_line
            )
        if intersected:
            return closest_pt
        else:
            return None
    
    def BoundingBox( self ):
        '''
        Returns a pair of points, the maximum x,y,z corner of the xyz-axis-aligned bounding box
        containing the ellipse, and the minimum x,y,z corner.
        '''
        max_abs_axes = abs( self.axes ).max(0)
        return (
            self.center + max_abs_axes,
            self.center - max_abs_axes
            )
    
    def glTransform( self ):
        '''
        Multiplies the matrix stack by the transform to take points in 
        '''
        glTranslated( *self.center )
        
        ## This is the sphere to ellipse rotation and scale matrix.
        rot_scale_mat = zeros( (4,4) )
        rot_scale_mat[ :3, :3 ] = self.axes.T
        rot_scale_mat[3,3] = 1.
        glMultMatrixd( rot_scale_mat.T.flatten() )
    
    def Draw( self ):
        
        glMatrixMode( GL_MODELVIEW )
        glPushMatrix()
        
        self.glTransform()
        
        # glutSolidSphere( 1., 40, 40 )
        from mesh_gui import draw_mesh_faces_smooth
        draw_mesh_faces_smooth( get_sphere_mesh() )
        
        glPopMatrix()
    
    def GenMesh( self ):
        '''
        Returns a TriMesh object sampling the ellipsoid.
        '''
        
        from copy import deepcopy
        result = deepcopy( get_sphere_mesh() )
        
        result.vs = dot( self.axes.T, asarray( result.vs ).T ).T + self.center
        
        ## Normals need to be transformed by the inverse transpose.
        inv_transpose = linalg.inv( self.axes.T ).T
        result.face_normals[:] = dot( inv_transpose, result.face_normals.T ).T
        result.vertex_normals[:] = dot( inv_transpose, result.vertex_normals.T ).T
        
        return result

_sphere_mesh = None
def get_sphere_mesh():
    global _sphere_mesh
    
    if _sphere_mesh is None:
        
        sphere_fname = 'sphere4.out.obj'
        
        from mesh_gui import TriMesh
        import os.path
        mesh = TriMesh.FromOBJ_FileName( os.path.join( '.', 'ext', sphere_fname ) )
        mesh.vs = asarrayf( mesh.vs )
        mesh.faces = asarray( mesh.faces, dtype = int )
        from math import atan2, atan
        ## TODO: Seams.
        mesh.uvs = [
            ( .5 + .5 * atan2( y, x ) / pi, .5 + .5 * atan( z / sqrt( x**2 + y**2 ) ) / (pi/2.) )
            for x,y,z in mesh.vs
            ]
        
        ## Call all the accessors to precalculate them.
        mesh.get_face_normals()
        mesh.get_face_areas()
        mesh.get_vertex_normals()
        mesh.get_vertex_areas()
        mesh.get_edges()
        mesh.get_halfedges()
        
        _sphere_mesh = mesh
    
    return _sphere_mesh
