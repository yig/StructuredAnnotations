from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *

class Environment( object ):
    def __init__( self ):
        ## This is the spring constant
        self.k = 1.
        ## This is the bending constant
        self.b = 1.
        ## This is the gravitational constant
        #self.g = .01
        self.g = .02
        #self.g = .021
        ## This is the area expansion constant
        self.a = 0.
        ## This is the pressure in the material
        self.p_m = 1.
        ## This is the external (atmospheric) pressure
        self.p_atm = 1.
        ## This is the curvature^2 energy
        self.c = 1.
        ## This is the laplacian editing energy
        self.l = 1.
        
        self.kGravity = False
        self.kStretching = True
        self.kBending = False
        self.kExpansion = True
        self.kPressure = False
        self.kCurvature = False
        self.kCurvatureUndeformed = False
        self.kLaplacianEditing = False
        self.kLaplacianEditingUndeformed = False
    
    def DisableAllForces( self ):
        self.kGravity = False
        self.kStretching = False
        self.kBending = False
        self.kExpansion = False
        self.kPressure = False
        self.kCurvature = False
        self.kCurvatureUndeformed = False
        self.kLaplacianEditing = False
        self.kLaplacianEditingUndeformed = False

## This is our threshold for whether a number is approximately zero.
ZERO_THRESHOLD = 1e-8

kGradientDescent = False
#kGradientDescent = True

from helpers import x_intercept

def triangle_area( p_i, p_j, p_k ):
    '''
    Returns the triangle area of CCW triangle < p_i, p_j, p_k >.
    '''
    
    A = .5 * ( ( p_j[0] - p_i[0] )*( p_k[1] - p_j[1] ) - ( p_k[0] - p_j[0] )*( p_j[1] - p_i[1] ) )

    if A < 0:
        #print A
        A = 0

    #assert A >= 0.
    return A

def triangle_vertex_normal( p_i, p_j, p_k ):
    '''
    Returns the outward facing normal for vertex p_i of CCW triangle < p_i, p_j, p_k >.
    '''
    
    p_ij = p_i - p_j
    len_p_ij = sqrt( sum( p_ij ** 2 ) )
    if abs( len_p_ij ) > ZERO_THRESHOLD:
        p_ij *= 1/len_p_ij
    
    p_ik = p_i - p_k
    len_p_ik = sqrt( sum( p_ik ** 2 ) )
    if abs( len_p_ik ) > ZERO_THRESHOLD:
        p_ik *= 1/len_p_ik
    
    p_ijik = p_ij + p_ik
    len_p_ijik = sqrt( sum( p_ijik ** 2 ) )
    if abs( len_p_ijik ) > ZERO_THRESHOLD:
        p_ijik *= 1/len_p_ijik
    
    return p_ijik

def triangle_vertex_perimeter( p_i, p_j, p_k ):
    '''
    Returns the triangle's perimeter owned by vertex p_i of CCW triangle < p_i, p_j, p_k >.
    '''
    
    p_ij = p_i - p_j
    len_p_ij = sqrt( sum( p_ij ** 2 ) )
    
    p_ik = p_i - p_k
    len_p_ik = sqrt( sum( p_ik ** 2 ) )
    
    result = .5 * ( len_p_ij + len_p_ik )
    #assert result >= 0.
    if result < 0: result = 0
    return result

def mesh_area( p, triangles ):
    '''
    Returns the total area of the mesh given by positions 'p' and a list of (i,j,k) 'triangles'.
    '''
    
    A = 0.
    for i,j,k in triangles: A += triangle_area( p[i], p[j], p[k] )
    assert A >= 0.
    return A

def Fpressure_ijk( p,u, i,j,k, env ):
    '''
    Returns the pressure force of the triangle p_i, p_j, p_k acting on the node p_i.
    To get the force acting on the node p_j (or p_k), cycle the input parameters so that
    p_j is the first parameter.
    '''
    
    A0 = triangle_area( u[i], u[j], u[k] )
    A = triangle_area( p[i], p[j], p[k] )
    
    nijk = triangle_vertex_normal( p[i], p[j], p[k] )
    sijk = triangle_vertex_perimeter( p[i], p[j], p[k] )
    
    return ( ( env.p_m * A0/A - env.p_atm ) * sijk ) * nijk

def dFpressure_ijk_dpl( p,u, i,j,k, l, env ):
    '''
    Returns the derivative of the pressure force of the triangle p_i, p_j, p_k,
    acting on the node p_i, with respect to node p_l.
    This is zero unless l == i, l == j, or l == k.
    '''
    
    if l not in (i,j,k): return zerosf( (2,2) )
    
    from springs_backwards_adt import dFpressure_ijk_dp_ijk
    l_index = [i,j,k].index(l)
    return dFpressure_ijk_dp_ijk(
        p[i], p[j], p[k],
        u[i], u[j], u[k],
        env.p_m, env.p_atm
        )[ :, 2*l_index : 2*(l_index+1) ]

def dFpressure_ijk_dul( p,u, i,j,k, l, env ):
    '''
    Returns the derivative of the pressure force of the triangle p_i, p_j, p_k,
    acting on the node p_i, with respect to node u_l.
    This is zero unless l == i, l == j, or l == k.
    '''
    
    if l not in (i,j,k): return zerosf( (2,2) )
    
    from springs_backwards_adt import dFpressure_ijk_du_ijk
    l_index = [i,j,k].index(l)
    return dFpressure_ijk_du_ijk(
        p[i], p[j], p[k],
        u[i], u[j], u[k],
        env.p_m, env.p_atm
        )[ :, 2*l_index : 2*(l_index+1) ]

expansion_is_density = False
if expansion_is_density:
    expansion_density_area_is_undeformed = False
    
    def FDexpansion_ijk( p_i, p_j, p_k, env ):
        '''
        Returns the expansion force density of the triangle p_i, p_j, p_k acting on the node p_i.
        To get the force acting on the node p_j (or p_k), cycle the input parameters so that
        p_j is the first parameter.
        '''
        
        return env.a * .5 * arrayf( (
            p_j[1] - p_k[1],
            p_k[0] - p_j[0]
            ) )
    
    def Fexpansion_ijk( p,u, i,j,k, env ):
        '''
        Returns the expansion force of the triangle p_i, p_j, p_k acting on the node p_i.
        To get the force acting on the node p_j (or p_k), cycle the input parameters so that
        p_j is the first parameter.
        '''
        
        if expansion_density_area_is_undeformed:
            A = triangle_area( u[i], u[j], u[k] )
        else:
            A = triangle_area( p[i], p[j], p[k] )
        
        return A/3. * FDexpansion_ijk( p[i], p[j], p[k], env )
    
    def dFDexpansion_ijk_dl( i, j, k, l, env ):
        '''
        Returns the derivative of the expansion force density of the triangle p_i, p_j, p_k,
        acting on the node p_i, with respect to node p_l.
        This is zero unless l == j or l == k.
        '''
        
        ## NOTE: This should be static.
        result = zerosf( (2,2) )
        result[ 0,1 ] = .5 * env.a
        result[ 1,0 ] = -.5 * env.a
        
        if l == j: return result
        elif l == k: return -result
        else: return zerosf( (2,2) )
    
    def dFexpansion_ijk_dpl( p,u, i,j,k, l, env ):
        '''
        Returns the derivative of the expansion force of the triangle p_i, p_j, p_k,
        acting on the node p_i, with respect to node p_l.
        This is zero unless l == j or l == k.
        '''
        
        tri_p = [ p[i], p[j], p[k], p[i], p[j] ]
        if l == i: loff = 0
        elif l == j: loff = 1
        elif l == k: loff = 2
        
        if expansion_density_area_is_undeformed:
            A = triangle_area( u[i], u[j], u[k] )
            result = (
                A/3. * dFDexpansion_ijk_dl( i,j,k, l, env )
                )
        else:
            A = triangle_area( p[i], p[j], p[k] )
            result = 1./3. * (
                A * dFDexpansion_ijk_dl( i,j,k, l, env )
                + outer( FDexpansion_ijk( *(tri_p[ loff : loff + 3 ] + [env]) ), FDexpansion_ijk( *(tri_p[ loff : loff + 3 ] + [env]) ) )
                )
        return result
    
    
    def dFexpansion_ijk_dul( p,u, i,j,k, l, env ):
        '''
        Returns the derivative of the expansion force of the triangle p_i, p_j, p_k,
        acting on the node p_i, with respect to node u_l.
        This is zero unless l == j or l == k.
        '''
        
        tri_p = [ p[i], p[j], p[k], p[i], p[j] ]
        tri_u = [ u[i], u[j], u[k], u[i], u[j] ]
        if l == i: loff = 0
        elif l == j: loff = 1
        elif l == k: loff = 2
        
        if expansion_density_area_is_undeformed:
            A = triangle_area( p[i], p[j], p[k] )
            result = 1./3. * outer(
                FDexpansion_ijk( *(tri_u[ loff : loff + 3 ] + [env]) ),
                FDexpansion_ijk( *(tri_p[ loff : loff + 3 ] + [env]) )
                )
        else:
            result = zerosf( (2,2) )
        return result
else:
    def Fexpansion_ijk( p,u, i,j,k, env ):
        '''
        Returns the expansion force of the triangle p_i, p_j, p_k acting on the node p_i.
        To get the force acting on the node p_j (or p_k), cycle the input parameters so that
        p_j is the first parameter.
        '''
        
        return env.a * .5 * arrayf( (
            p[j][1] - p[k][1],
            p[k][0] - p[j][0]
            ) )
    
    def dFexpansion_ijk_dpl( p,u, i,j,k, l, env ):
        '''
        Returns the derivative of the expansion force of the triangle p_i, p_j, p_k,
        acting on the node p_i, with respect to node p_l.
        This is zero unless l == j or l == k.
        '''
        
        ## NOTE: This should be static.
        result = zerosf( (2,2) )
        result[ 0,1 ] = .5 * env.a
        result[ 1,0 ] = -.5 * env.a
        
        if l == j: return result
        elif l == k: return -result
        else: return zerosf( (2,2) )
    
    def dFexpansion_ijk_dul( p,u, i,j,k, l, env ):
        '''
        Returns the derivative of the expansion force of the triangle p_i, p_j, p_k,
        acting on the node p_i, with respect to node u_l.
        This is always zero.
        '''
        
        return zerosf( (2,2) )

def Fstretching_ij( p_i, p_j, u_i, u_j, env ):
    '''
    Returns the force of the spring from 'p_i' to 'p_j' with rest length r_ij = |u_i - u_j|
    acting on 'p_i'.
    Note that Fstretching_ij( p_i, p_j, u_i, u_j ) equals -Fstretching_ij( p_j, p_i, u_j, u_i ).
    '''
    
    p_ij = p_i - p_j
    len_p_ij = sqrt( sum( p_ij ** 2 ) )
    r_ij = sqrt( sum( (u_i - u_j) ** 2 ) )
    
    if abs( len_p_ij ) < ZERO_THRESHOLD:
        result = 0. * p_ij
    else:
        result = -env.k * ( len_p_ij - r_ij ) / len_p_ij * p_ij
    
    return result

def F( p, edges, u, triangles, env ):
    '''
    Returns a vector containing the force at every point.
    Note that the input 'p' is assumed to be a number-of-points by dimension array, where
    dimension is 2 for our example.  The result is flattened into a single vector
    of size number-of-points times dimension.
    '''
    
    dim = p.shape[1]
    
    Fp = zerosf( prod( p.shape ) )
    
    if env.kStretching:
        for (i,j) in edges:
            assert i != j
            
            Fij = Fstretching_ij( p[i], p[j], u[i], u[j], env )
            
            Fp[ i*dim : (i+1) * dim ] += Fij
            ## 'edges' contains edges uniquely, so 'edges' will contain (i,j) but not (j,i).
            ## This means that we must add -Fij to row j as well as Fij to row i.
            Fp[ j*dim : (j+1) * dim ] += -Fij
    
    if env.kGravity:
        for i in range( len(p) ):
            Fp[ i*dim : (i+1) * dim ] += -env.g * arrayf( (0,1) )
    
    if env.kBending:
        for i in range( 1,len(p)-1 ):
            Fp[ i*dim : (i+1) * dim ] += Fbending_i( p[i-1], p[i], p[i+1], env )
    
    if env.kCurvature:
        import Eline_ad
        ## I'm going to assume the dimension is 2 and append 0. to every point.
        assert dim == 2
        for i in range( 1,len(p)-1 ):
            E, grad = Eline_ad.Ecurvature_and_derivatives( append( p[ i-1 ], 0. ), append( p[ i ], 0. ), append( p[ i+1 ], 0. ), compute_hessian = False )
            Fp[ (i-1)*dim : (i+2) * dim ] += env.c * grad[ [0,1,3,4,6,7] ]
    
    if env.kCurvatureUndeformed:
        import Eline_ad
        ## I'm going to assume the dimension is 2 and append 0. to every point.
        assert dim == 2
        for i in range( 1,len(p)-1 ):
            E, grad = Eline_ad.Ecurvature_and_derivatives_undeformed(
                append( p[ i-1 ], 0. ), append( p[ i ], 0. ), append( p[ i+1 ], 0. ),
                append( u[ i-1 ], 0. ), append( u[ i ], 0. ), append( u[ i+1 ], 0. ),
                compute_hessian = False
                )
            Fp[ (i-1)*dim : (i+2) * dim ] += env.c * grad[ [0,1,3,4,6,7] ]
    
    if env.kLaplacianEditing:
        import Eline_ad
        ## I'm going to assume the dimension is 2 and append 0. to every point.
        assert dim == 2
        for i in range( 1,len(p)-1 ):
            E, grad = Eline_ad.Elaplacian_editing_and_derivatives( append( p[ i-1 ], 0. ), append( p[ i ], 0. ), append( p[ i+1 ], 0. ), compute_hessian = False )
            Fp[ (i-1)*dim : (i+2) * dim ] += env.l * grad[ [0,1,3,4,6,7] ]
    
    if env.kLaplacianEditingUndeformed:
        import Eline_ad
        ## I'm going to assume the dimension is 2 and append 0. to every point.
        assert dim == 2
        for i in range( 1,len(p)-1 ):
            E, grad = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                append( p[ i-1 ], 0. ), append( p[ i ], 0. ), append( p[ i+1 ], 0. ),
                append( u[ i-1 ], 0. ), append( u[ i ], 0. ), append( u[ i+1 ], 0. ),
                compute_hessian = False
                )
            Fp[ (i-1)*dim : (i+2) * dim ] += env.l * grad[ [0,1,3,4,6,7] ]
    
    if env.kExpansion:
        for (i,j,k) in triangles:
            assert i != j and i != k and j != k
            
            #print (i,j,k), Fexpansion_ijk( p[i], p[j], p[k], env )
            Fp[ i*dim : (i+1) * dim ] += Fexpansion_ijk( p,u, i, j, k, env )
            Fp[ j*dim : (j+1) * dim ] += Fexpansion_ijk( p,u, j, k, i, env )
            Fp[ k*dim : (k+1) * dim ] += Fexpansion_ijk( p,u, k, i, j, env )
    
    if env.kPressure:
        for (i,j,k) in triangles:
            assert i != j and i != k and j != k
            
            Fp[ i*dim : (i+1) * dim ] += Fpressure_ijk( p,u, i, j, k, env )
            Fp[ j*dim : (j+1) * dim ] += Fpressure_ijk( p,u, j, k, i, env )
            Fp[ k*dim : (k+1) * dim ] += Fpressure_ijk( p,u, k, i, j, env )
    
    return Fp

def Fbending_i( p_i_m_1, p_i, p_i_p_1, env ):
    return env.b * ( p_i_m_1 - 2 * p_i + p_i_p_1 )

def dFbending_i_dp_j( p_i_m_1, p_i, p_i_p_1, i, j, env ):
    if j == i: return -2 * env.b * identityf( 2 )
    if j == i+1: return env.b * identityf( 2 )
    if j == i-1: return env.b * identityf( 2 )
    else: return zerosf( (2,2) )

def dFstretching_ij_dp_i( p_i, p_j, u_i, u_j, env ):
    '''
    Returns the derivative with respect to 'p_i' of the force of the spring
    from 'p_i' to 'p_j' with rest length r_ij = |u_i - u_j| acting on 'p_i'.
    Our dimension is 2, so this is a 2x2 quantity.
    Note that dFstretching_ij_dp_i( p_i, p_j, u_i, u_j ) equals dFstretching_ij_dp_i( p_j, p_i, u_j, u_i ).
    '''
    
    dim = p_i.shape[0]
    
    p_ij = p_i - p_j
    len_p_ij = sqrt( sum( p_ij ** 2 ) )
    r_ij = sqrt( sum( (u_i - u_j) ** 2 ) )
    
    if abs( len_p_ij ) < ZERO_THRESHOLD:
        result = -env.k * identity( dim )
    else:
        result = env.k * ( -identity( dim ) - r_ij / len_p_ij**3 * outer( p_ij, p_ij ) + r_ij / len_p_ij * identity( dim ) )
    return result

def dFstretching_ij_du_i( p_i, p_j, u_i, u_j, env ):
    '''
    Returns the derivative with respect to 'u_i' of the force of the spring
    from 'p_i' to 'p_j' with rest length r_ij = |u_i - u_j| acting on 'p_i'.
    Our dimension is 2, so this is a 2x2 quantity.
    Note that dFstretching_ij_du_i( p_i, p_j, u_i, u_j ) equals -dFstretching_ij_du_j( p_i, p_j, u_i, u_j ).
    Note that dFstretching_ij_du_i( p_i, p_j, u_i, u_j ) equals dFstretching_ij_du_i( p_j, p_i, u_j, u_i ).
    '''
    
    dim = p_i.shape[0]
    
    p_ij = p_i - p_j
    len_p_ij = sqrt( sum( p_ij ** 2 ) )
    
    u_ij = u_i - u_j
    r_ij = sqrt( sum( u_ij ** 2 ) )
    
    if abs( len_p_ij ) < ZERO_THRESHOLD:
        result = zerosf( (dim,dim) )
    else:
        result = env.k / (len_p_ij * r_ij) * outer( p_ij, u_ij )
    return result

def Jp( p, edges, u, triangles, env ):
    '''
    Returns a matrix containing the derivative of the force at every point with respect to each point.
    Note that the inputs 'p' and 'u' are assumed to be a number-of-points by dimension array, where
    dimension is 2 for our example.
    'p' contains the deformed configuration of the points, while 'u' contains the undeformed
    configuration.
    The result is flattened is a square matrix (of type numpy.array), of size NxN, where
    N = number-of-points times dimension.
    '''
    
    assert p.shape == u.shape
    
    dim = p.shape[1]
    
    if kGradientDescent: return -10*identityf( prod( p.shape ) )
    
    Jp = zerosf( ( prod( p.shape ), prod( p.shape ) ) )
    
    if env.kStretching:
        ## stretching
        for (i,j) in edges:
            assert i != j
            
            dF = dFstretching_ij_dp_i( p[i], p[j], u[i], u[j], env )
            assert ( ( Jp[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] - zeros( ( dim, dim ) ) ) ** 2 ).sum().sum() == 0.
            
            Jp[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] = -dF
            Jp[ i*dim : (i+1) * dim, i*dim : (i+1) * dim ] += dF
            ## 'edges' contains edges uniquely, so 'edges' will contain (i,j) but not (j,i).
            ## This means that we must add dF to the right places in column j as well.
            Jp[ j*dim : (j+1) * dim, i*dim : (i+1) * dim ] = -dF
            Jp[ j*dim : (j+1) * dim, j*dim : (j+1) * dim ] += dF    
    
    if env.kGravity:
        ## gravity
        
        ## Jacobian of gravity is zero
        pass
        
        #for i in range( len(p) ):
        #    Jp[ i*dim : (i+1) * dim, i*dim : (i+1) * dim ] += zerosf( (2,2) )
    
    if env.kBending:
        ## bending (naive)
        #for i in range( len(p) ):
        #    for j in range( len(p) ):
        #        Jp[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] += dFbending_i_dp_j( p[i-1], p[i], p[i+1], i, j, env )
        
        ## bending
        for i in range( 1, len(p)-1 ):
            for j in range( i-1, i+1+1 ):
                Jp[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] += dFbending_i_dp_j( p[i-1], p[i], p[i+1], i, j, env )
    
    if env.kCurvature:
        import Eline_ad
        ## I'm going to assume the dimension is 2 and append 0. to every point.
        assert dim == 2
        for i in range( 1,len(p)-1 ):
            E, grad, hess = Eline_ad.Ecurvature_and_derivatives( append( p[ i-1 ], 0. ), append( p[ i ], 0. ), append( p[ i+1 ], 0. ), compute_hessian = True )
            #print "grad(i-1,i,i+1):"
            #print grad.round(3)
            #print "hess(i-1,i,i+1):"
            #print hess.round(3)
            ## TODO: Can I replace this double for loop with a reshape into block 2x2 and block 3x3 matrices?
            for ci in range(dim):
                for cj in range(dim):
                    Jp[ (i-1)*dim + ci : (i+2) * dim : dim, (i-1)*dim + cj : (i+2)*dim : dim ] += env.c * hess[ci::3,cj::3]
    
    if env.kCurvatureUndeformed:
        import Eline_ad
        ## I'm going to assume the dimension is 2 and append 0. to every point.
        assert dim == 2
        for i in range( 1,len(p)-1 ):
            E, grad, hess = Eline_ad.Ecurvature_and_derivatives_undeformed(
                append( p[ i-1 ], 0. ), append( p[ i ], 0. ), append( p[ i+1 ], 0. ),
                append( u[ i-1 ], 0. ), append( u[ i ], 0. ), append( u[ i+1 ], 0. ),
                compute_hessian = True
                )
            #print "grad(i-1,i,i+1):"
            #print grad.round(3)
            #print "hess(i-1,i,i+1):"
            #print hess.round(3)
            ## TODO: Can I replace this double for loop with a reshape into block 2x2 and block 3x3 matrices?
            for ci in range(dim):
                for cj in range(dim):
                    Jp[ (i-1)*dim + ci : (i+2) * dim : dim, (i-1)*dim + cj : (i+2)*dim : dim ] += env.c * hess[ci::3,cj::3]
    
    if env.kLaplacianEditing:
        import Eline_ad
        ## I'm going to assume the dimension is 2 and append 0. to every point.
        assert dim == 2
        for i in range( 1,len(p)-1 ):
            E, grad, hess = Eline_ad.Elaplacian_editing_and_derivatives( append( p[ i-1 ], 0. ), append( p[ i ], 0. ), append( p[ i+1 ], 0. ), compute_hessian = True )
            #print "grad(i-1,i,i+1):"
            #print grad.round(3)
            #print "hess(i-1,i,i+1):"
            #print hess.round(3)
            ## TODO: Can I replace this double for loop with a reshape into block 2x2 and block 3x3 matrices?
            for ci in range(dim):
                for cj in range(dim):
                    Jp[ (i-1)*dim + ci : (i+2) * dim : dim, (i-1)*dim + cj : (i+2)*dim : dim ] += env.l * hess[ci::3,cj::3]
    
    if env.kLaplacianEditingUndeformed:
        import Eline_ad
        ## I'm going to assume the dimension is 2 and append 0. to every point.
        assert dim == 2
        for i in range( 1,len(p)-1 ):
            E, grad, hess = Eline_ad.Elaplacian_editing_and_derivatives_undeformed(
                append( p[ i-1 ], 0. ), append( p[ i ], 0. ), append( p[ i+1 ], 0. ),
                append( u[ i-1 ], 0. ), append( u[ i ], 0. ), append( u[ i+1 ], 0. ),
                compute_hessian = True
                )
            #print "grad(i-1,i,i+1):"
            #print grad.round(3)
            #print "hess(i-1,i,i+1):"
            #print hess.round(3)
            ## TODO: Can I replace this double for loop with a reshape into block 2x2 and block 3x3 matrices?
            for ci in range(dim):
                for cj in range(dim):
                    Jp[ (i-1)*dim + ci : (i+2) * dim : dim, (i-1)*dim + cj : (i+2)*dim : dim ] += env.l * hess[ci::3,cj::3]
    
    if env.kExpansion:
        for (i,j,k) in triangles:
            assert i != j and i != k and j != k
            
            if expansion_is_density:
                Jp[ i*dim : (i+1) * dim, i*dim : (i+1) * dim ] += dFexpansion_ijk_dpl( p,u, i,j,k, i, env )
                Jp[ j*dim : (j+1) * dim, j*dim : (j+1) * dim ] += dFexpansion_ijk_dpl( p,u, j,k,i, j, env )
                Jp[ k*dim : (k+1) * dim, k*dim : (k+1) * dim ] += dFexpansion_ijk_dpl( p,u, k,i,j, k, env )
            else:
                ## We don't compute the [i,i] block because we know it's zero.
                pass
            
            Jp[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] += dFexpansion_ijk_dpl( p,u, i,j,k, j, env )
            Jp[ i*dim : (i+1) * dim, k*dim : (k+1) * dim ] += dFexpansion_ijk_dpl( p,u, i,j,k, k, env )
            
            Jp[ j*dim : (j+1) * dim, i*dim : (i+1) * dim ] += dFexpansion_ijk_dpl( p,u, j,k,i, i, env )
            Jp[ j*dim : (j+1) * dim, k*dim : (k+1) * dim ] += dFexpansion_ijk_dpl( p,u, j,k,i, k, env )
            
            Jp[ k*dim : (k+1) * dim, i*dim : (i+1) * dim ] += dFexpansion_ijk_dpl( p,u, k,i,j, i, env )
            Jp[ k*dim : (k+1) * dim, j*dim : (j+1) * dim ] += dFexpansion_ijk_dpl( p,u, k,i,j, j, env )
    
    if env.kPressure:
        for (i,j,k) in triangles:
            assert i != j and i != k and j != k
            
            Jp[ i*dim : (i+1) * dim, i*dim : (i+1) * dim ] += dFpressure_ijk_dpl( p,u, i,j,k, i, env )
            Jp[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] += dFpressure_ijk_dpl( p,u, i,j,k, j, env )
            Jp[ i*dim : (i+1) * dim, k*dim : (k+1) * dim ] += dFpressure_ijk_dpl( p,u, i,j,k, k, env )
            
            Jp[ j*dim : (j+1) * dim, j*dim : (j+1) * dim ] += dFpressure_ijk_dpl( p,u, j,k,i, j, env )
            Jp[ j*dim : (j+1) * dim, i*dim : (i+1) * dim ] += dFpressure_ijk_dpl( p,u, j,k,i, i, env )
            Jp[ j*dim : (j+1) * dim, k*dim : (k+1) * dim ] += dFpressure_ijk_dpl( p,u, j,k,i, k, env )
            
            Jp[ k*dim : (k+1) * dim, k*dim : (k+1) * dim ] += dFpressure_ijk_dpl( p,u, k,i,j, k, env )
            Jp[ k*dim : (k+1) * dim, i*dim : (i+1) * dim ] += dFpressure_ijk_dpl( p,u, k,i,j, i, env )
            Jp[ k*dim : (k+1) * dim, j*dim : (j+1) * dim ] += dFpressure_ijk_dpl( p,u, k,i,j, j, env )
    
    #print '|Jp.sum(1)|', mag( Jp.sum(1).ravel() )
    return Jp

def Ju( p, edges, u, triangles, env ):
    '''
    Returns a matrix containing the derivative of the force at every point with respect to each
    undeformed point.
    Note that the input 'p' is assumed to be a number-of-points by dimension array, where
    dimension is 2 for our example.
    The result is flattened is a square matrix (of type numpy.array), of size NxN, where
    N = number-of-points times dimension.
    '''
    
    assert p.shape == u.shape
    
    dim = p.shape[1]
    
    Ju = zerosf( ( prod( p.shape ), prod( p.shape ) ) )
    
    if env.kStretching:
        for (i,j) in edges:
            assert i != j
            
            dF_ij_du_i = dFstretching_ij_du_i( p[i], p[j], u[i], u[j], env )
            ## NOTE:
            ##  dF_ij_du_j == -dF_ij_du_i
            ##  dF_ji_du_j = dF_ij_du_i
            ##  dF_ji_du_i = -dF_ij_du_i
            assert ( ( Ju[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] - zeros( ( dim, dim ) ) ) ** 2 ).sum().sum() == 0.
            
            Ju[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] = -dF_ij_du_i
            Ju[ i*dim : (i+1) * dim, i*dim : (i+1) * dim ] += dF_ij_du_i
            ## 'edges' contains edges uniquely, so 'edges' will contain (i,j) but not (j,i).
            ## This means that we must add dF to the right places in column j as well.
            Ju[ j*dim : (j+1) * dim, i*dim : (i+1) * dim ] = -dF_ij_du_i
            Ju[ j*dim : (j+1) * dim, j*dim : (j+1) * dim ] += dF_ij_du_i
    
    if env.kGravity:
        ## jacobian of gravity is zero
        pass
        ## gravity
        #for i in range( len(p) ):
        #    Jp[ i*dim : (i+1) * dim, i*dim : (i+1) * dim ] += zerosf( (2,2) )
    
    if env.kBending:
        pass
        ## bending (naive)
        #for i in range( len(p) ):
        #    for j in range( len(p) ):
        #        Jp[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] += dFbending_i_dp_j( p[i-1], p[i], p[i+1], i, j, env )
        
        ## bending
        #for i in range( 1, len(p)-1 ):
        #    for j in range( i-1, i+1+1 ):
        #        Jp[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] += dFbending_i_dp_j( p[i-1], p[i], p[i+1], i, j, env )
    
    if env.kCurvature:
        raise NotImplementedError
    
    if env.kCurvatureUndeformed:
        raise NotImplementedError
    
    if env.kLaplacianEditing:
        raise NotImplementedError
    
    if env.kLaplacianEditingUndeformed:
        raise NotImplementedError
    
    if env.kExpansion:
        for (i,j,k) in triangles:
            assert i != j and i != k and j != k
            
            if expansion_is_density:
                Ju[ i*dim : (i+1) * dim, i*dim : (i+1) * dim ] += dFexpansion_ijk_dul( p,u, i,j,k, i, env )
                Ju[ j*dim : (j+1) * dim, j*dim : (j+1) * dim ] += dFexpansion_ijk_dul( p,u, j,k,i, j, env )
                Ju[ k*dim : (k+1) * dim, k*dim : (k+1) * dim ] += dFexpansion_ijk_dul( p,u, k,i,j, k, env )
            else:
                ## We don't compute the [i,i] block because we know it's zero.
                pass
            
            Ju[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] += dFexpansion_ijk_dul( p,u, i,j,k, j, env )
            Ju[ i*dim : (i+1) * dim, k*dim : (k+1) * dim ] += dFexpansion_ijk_dul( p,u, i,j,k, k, env )
            
            Ju[ j*dim : (j+1) * dim, i*dim : (i+1) * dim ] += dFexpansion_ijk_dul( p,u, j,k,i, i, env )
            Ju[ j*dim : (j+1) * dim, k*dim : (k+1) * dim ] += dFexpansion_ijk_dul( p,u, j,k,i, k, env )
            
            Ju[ k*dim : (k+1) * dim, i*dim : (i+1) * dim ] += dFexpansion_ijk_dul( p,u, k,i,j, i, env )
            Ju[ k*dim : (k+1) * dim, j*dim : (j+1) * dim ] += dFexpansion_ijk_dul( p,u, k,i,j, j, env )
    
    if env.kPressure:
        for (i,j,k) in triangles:
            assert i != j and i != k and j != k
            
            Ju[ i*dim : (i+1) * dim, i*dim : (i+1) * dim ] += dFpressure_ijk_dul( p,u, i,j,k, i, env )
            Ju[ i*dim : (i+1) * dim, j*dim : (j+1) * dim ] += dFpressure_ijk_dul( p,u, i,j,k, j, env )
            Ju[ i*dim : (i+1) * dim, k*dim : (k+1) * dim ] += dFpressure_ijk_dul( p,u, i,j,k, k, env )
            
            Ju[ j*dim : (j+1) * dim, j*dim : (j+1) * dim ] += dFpressure_ijk_dul( p,u, j,k,i, j, env )
            Ju[ j*dim : (j+1) * dim, i*dim : (i+1) * dim ] += dFpressure_ijk_dul( p,u, j,k,i, i, env )
            Ju[ j*dim : (j+1) * dim, k*dim : (k+1) * dim ] += dFpressure_ijk_dul( p,u, j,k,i, k, env )
            
            Ju[ k*dim : (k+1) * dim, k*dim : (k+1) * dim ] += dFpressure_ijk_dul( p,u, k,i,j, k, env )
            Ju[ k*dim : (k+1) * dim, i*dim : (i+1) * dim ] += dFpressure_ijk_dul( p,u, k,i,j, i, env )
            Ju[ k*dim : (k+1) * dim, j*dim : (j+1) * dim ] += dFpressure_ijk_dul( p,u, k,i,j, j, env )
    
    #print '|Ju.sum(1)|', mag( Ju.sum(1).ravel() )
    return Ju

def constrain_system( A, rhs, rows ):
    '''
    This function modifies its input parameters, a system matrix 'A' and
    right-hand-side vector 'rhs', such that for every index i in 'rows',
    the row i of A is set to row i of the identity matrix and rhs[i] is set to zero.
    '''
    
    for i in rows:
        A[ i, : ] = zeros( A.shape[1] )
        ## We can also zero the column, which keeps the matrix symmetric, because
        ## we are zeroing the corresponding entries in the right-hand-side (x*0 = 0).
        A[ :, i ] = zeros( A.shape[0] )
        A[ i, i ] = 1
        rhs[i] = 0
    
    return A, rhs

def solve_in_steps( solve_func, result_converged_func, lerp_args_func, min_step = 1e-2, verbose = False ):
    '''
    Meta-routine to find a solution which can't be reached immediately, but only slowly, in steps.
    The parameter 'lerp_args_func' is a function that takes a value from 0 to 1,
    where 0 is guaranteed to converge and 1 may not.
    The return value of 'lerp_args_func' is passed to the function parameter 'solve_func'.
    The return value of 'solve_func' is passed to the function parameter 'result_converged_func',
    which returns True if the call to 'solve_func' converged, and False otherwise.
    The routine attempts to reach successful convergence with 'lerp_args_func(1)', but
    terminates if the lerp parameter step decreases below 'min_step'.
    The routine always returns the result of the last call to 'solve_func'; this return value
    can be tested for convergence.
    
    TODO Q: Should I return a tuple, ( last converged lerp value, corresponding result of 'solve_func' )?
    '''
    
    assert min_step > 0.
    
    converged_t = 0.
    if verbose: print('solve_in_steps(): Attempting lerp value ', converged_t)
    last_result = solve_func( lerp_args_func( converged_t ) )
    
    ## 0 must always succeed.
    # assert result_converged_func( last_result )
    if not result_converged_func( last_result ): return last_result
    if verbose: print('solve_in_steps(): Lerp value converged')
    
    attempt_t = 1.
    while attempt_t - converged_t > min_step:
        assert attempt_t > converged_t
        
        if verbose: print('solve_in_steps(): Attempting lerp value ', attempt_t)
        attempt_t_args = lerp_args_func( attempt_t )
        last_result = solve_func( attempt_t_args )
        if result_converged_func( last_result ):
            if verbose: print('solve_in_steps(): Lerp value converged')
            
            ## Go for the goal.
            next_attempt_t = 1.
            ## TODO Q: Should I do something like:
            # next_attempt_t = min( 1., attempt_t + 1.5 * ( attempt_t - converged_t ) )
            ## A1: No, because attempting to jump directly to 1 seems to work pretty well.
            
            converged_t = attempt_t
            attempt_t = next_attempt_t
        else:
            if verbose: print('solve_in_steps(): Lerp value diverged')
            
            ## Decrease the step by half
            attempt_t = .5 * ( converged_t + attempt_t )
    
    return last_result

def static_forwards_solution_in_steps( p, edges, u, triangles, constraints, env, **kwargs ):
    '''
    Takes the same arguments as static_forwards_solution(), with the addition of an
    'in_steps_min_step' parameter, which gets passed as 'min_step' to solve_in_steps().
    '''
    
    ## Extract the keyword arguments for solve_in_steps().
    solve_in_steps_kwargs = {}
    if 'in_steps_min_step' in kwargs:
        solve_in_steps_kwargs[ 'min_step' ] = kwargs[ 'in_steps_min_step' ]
        del kwargs[ 'in_steps_min_step' ]
    if 'verbose' in kwargs:
        solve_in_steps_kwargs[ 'verbose' ] = kwargs[ 'verbose' ]
        del kwargs['verbose']
    
    ## Prepare the arguments to solve_in_steps().
    deformed = p.copy()
    def solve_func( env ):
        result = static_forwards_solution( deformed, edges, u, triangles, constraints, env, **kwargs )
        if result[1]: deformed.ravel()[:] = result[0].ravel()
        return result
    def result_converged_func( solve_func_result ):
        return solve_func_result[1]
    from copy import copy
    def lerp_args_func( t ):
        ## The only argument which gives us trouble is env.p_m, so we'll lerp that.
        args = copy( env )
        def lerp( a, b, t ): return a + t * ( b - a )
        args.p_m = lerp( env.p_atm, env.p_m, t )
        return args
    
    ## Call solve_in_steps().
    return solve_in_steps( solve_func, result_converged_func, lerp_args_func, **solve_in_steps_kwargs )

def static_forwards_solution( p, edges, u, triangles, constraints, env, **kwargs ):
    '''
    Given a list of points 'p' as an n-by-2 array, a list of (i,j) pairs 'edges' denoting an edge
    between points p[i] and p[j], a list of undeformed points 'u',
    and a list of position constraints specifying that certain p[i] don't change,
    uses Newton's method to solve for the positions where the forces are all zero.
    
    Returns two items: the updated positions (even if it diverges) and a boolean which
    is True if we converged.
    
    NOTE: 'edges' must not have both (i,j) and (j,i)
    '''
    
    p = asarray( p )
    u = asarray( u )
    assert p.shape[1] in (2,3)
    assert False not in [ i >= 0 and i < p.shape[0] for i, pos in constraints ]
    assert len( constraints ) == len( set( [ i for i, pos in constraints ] ) )
    assert False not in [ i >= 0 and j >= 0 and i < p.shape[0] and j < p.shape[0] for i,j in edges ]
    assert len( edges ) == len( set( [ frozenset( edge ) for edge in edges ] ) ) 
    assert False not in [ False not in [ i >= 0 and i < p.shape[0] for i in triangle ] for triangle in triangles ]
    assert len( triangles ) == len( set( [ frozenset( triangle ) for triangle in triangles ] ) )
    assert p.shape == u.shape
    
    kwargs.setdefault( 'verbose', False )
    kwargs.setdefault( 'xthresh', 1e-5 )
    kwargs.setdefault( 'fthresh', 1e-8 )
    kwargs.setdefault( 'max_iter', 100 )
    
    verbose = kwargs['verbose']
    
    XSTEP_THRESHOLD = kwargs['xthresh']
    F_THRESHOLD = kwargs['fthresh']
    MAX_ITERATIONS = kwargs['max_iter']
    
    p_n = p.copy().flatten()
    dim = p.shape[1]
    
    constrain_rows = []
    for i, pos in constraints:
        for c in range(len(pos)):
            if pos[c] is not None:
                constrain_rows.append( i*dim + c )
                p_n[ i*dim + c ] = pos[c]
    
    iteration = 0
    diverged = False
    while True:
        if verbose: print('-- iteration', iteration, '--')
        iteration += 1
        
        Jp_n = Jp( p_n.reshape( p.shape ), edges, u, triangles, env )
        Fp_n = F( p_n.reshape( p.shape ), edges, u, triangles, env )
        mag2_Fp_n = sum( Fp_n ** 2 )
        if verbose: print('| F( p_n ) |^2:', mag2_Fp_n)
        if mag2_Fp_n < F_THRESHOLD ** 2: break
        
        #debugger()
        if verbose:
            print("J:")
            print(Jp_n.round(3))
            print("F:")
            print(Fp_n.round(3))
            print("constraints:", constraints)
            print("constraint rows:", constrain_rows)
            
        constrain_system( Jp_n, Fp_n, constrain_rows )
        
        if verbose:
            print("J (constrained):")
            print(Jp_n.round(3))
            print("F (constrained):")
            print(Fp_n.round(3))
        
        # p_n_p_1 = p_n - dot( linalg.inv( Jp_n ), Fp_n )
        ## <=> p_n_p_1 - p_n = -linalg.inv( Jp_n ) * Fp_n
        ## <=> p_n - p_n_p_1 = linalg.inv( Jp_n ) * Fp_n
        ## <=> Jp_n * ( p_n - p_n_p_1 ) = Fp_n
        try:
            p_negative_delta = linalg.solve( Jp_n, Fp_n )
            ## p_n - ( p_n - p_n_p_1 ) = p_n_p_1
            p_n_p_1 = p_n - p_negative_delta
        except linalg.LinAlgError:
            print('LinAlgError')
            p_n_p_1 = p_n.copy()
            diverged = True
            break
        
        diff2 = sum( ( p_n_p_1 - p_n ) ** 2 )
        if verbose: print('| p_n+1 - p_n |^2:', diff2)
        p_n = p_n_p_1
        if diff2 < XSTEP_THRESHOLD ** 2: break
        
        if iteration >= MAX_ITERATIONS:
            print('Diverged.')
            diverged = True
            break
    
    return p_n.reshape( p.shape ), not diverged

def static_backwards_solution( p, edges, u, triangles, constraints, env, **kwargs ):
    '''
    Given a list of points 'p' as an n-by-2 array, a list of (i,j) pairs 'edges' denoting an edge
    between points p[i] and p[j], a list of undeformed points 'u',
    and a list of position constraints specifying that certain p[i] don't change,
    uses Newton's method to solve for the positions where the forces are all zero.
    
    Returns two items: the updated positions (even if it diverges) and a boolean which
    is True if we converged.
    
    NOTE: 'edges' must not have both (i,j) and (j,i)
    '''
    
    p = asarray( p )
    u = asarray( u )
    assert p.shape[1] in (2,3)
    assert False not in [ i >= 0 and i < p.shape[0] for i, pos in constraints ]
    assert len( constraints ) == len( set( [ i for i, pos in constraints ] ) )
    assert False not in [ i >= 0 and j >= 0 and i < p.shape[0] and j < p.shape[0] for i,j in edges ]
    assert len( edges ) == len( set( [ frozenset( edge ) for edge in edges ] ) ) 
    assert False not in [ False not in [ i >= 0 and i < p.shape[0] for i in triangle ] for triangle in triangles ]
    assert len( triangles ) == len( set( [ frozenset( triangle ) for triangle in triangles ] ) )
    assert p.shape == u.shape
    
    kwargs.setdefault( 'verbose', False )
    kwargs.setdefault( 'xthresh', 1e-5 )
    kwargs.setdefault( 'fthresh', 1e-8 )
    kwargs.setdefault( 'max_iter', 100 )
    
    verbose = kwargs['verbose']
    
    XSTEP_THRESHOLD = kwargs['xthresh']
    F_THRESHOLD = kwargs['fthresh']
    MAX_ITERATIONS = kwargs['max_iter']
    
    u_n = u.copy().flatten()
    dim = u.shape[1]
    
    constrain_rows = []
    for i, pos in constraints:
        for c in range(len(pos)):
            if pos[c] is not None:
                constrain_rows.append( i*dim + c )
                p_n[ i*dim + c ] = pos[c]
    
    iteration = 0
    diverged = False
    while True:
        if verbose: print('-- iteration', iteration, '--')
        iteration += 1
        
        Ju_n = Ju( p, edges, u_n.reshape( u.shape ), triangles, env )
        Fu_n = F( p, edges, u_n.reshape( u.shape ), triangles, env )
        mag2_Fu_n = sum( Fu_n ** 2 )
        if verbose: print('| F( u_n ) |^2:', mag2_Fu_n)
        if mag2_Fu_n < F_THRESHOLD: break
        
        constrain_system( Ju_n, Fu_n, constrain_rows )
        
        # u_n_p_1 = u_n - dot( linalg.inv( Ju_n ), Fu_n )
        ## <=> u_n_p_1 - u_n = -linalg.inv( Ju_n ) * Fu_n
        ## <=> u_n - u_n_p_1 = linalg.inv( Ju_n ) * Fu_n
        ## <=> Ju_n * ( u_n - u_n_p_1 ) = Fu_n
        try:
            u_negative_delta = linalg.solve( Ju_n, Fu_n )
            ## u_n - ( u_n - u_n_p_1 ) = u_n_p_1
            u_n_p_1 = u_n - u_negative_delta
        except linalg.LinAlgError:
            print('LinAlgError')
            u_n_p_1 = u_n.copy()
            diverged = True
            break
        
        diff2 = sum( ( u_n_p_1 - u_n ) ** 2 )
        if verbose: print('| u_n+1 - u_n |^2:', diff2)
        u_n = u_n_p_1
        if diff2 < XSTEP_THRESHOLD: break
        
        if iteration >= MAX_ITERATIONS:
            print('Diverged.')
            diverged = True
            break
    
    return u_n.reshape( u.shape ), not diverged

def static_backwards_solution_with_given_volume_expansion( p, edges, u, triangles, constraints, env, volume_expansion_factor, **kwargs ):
    '''
    Returns ( undeformed points, a, converged ) such that static_forwards_solution() called on the
    returned undeformed points and with env.a = a results in points 'p',
    and volume_expansion_factor * area( undeformed points ) = area(p).
    If 'converged' is False, then the process diverged and the returned values are the last
    step taken.
    '''
    
    p = asarray( p )
    u = asarray( u )
    from copy import copy
    env = copy( env )
    assert p.shape[1] in (2,3)
    assert False not in [ i >= 0 and j >= 0 and i < p.shape[0] and j < p.shape[0] for i,j in edges ]
    assert len( edges ) == len( set( [ frozenset( edge ) for edge in edges ] ) ) 
    assert False not in [ False not in [ i >= 0 and i < p.shape[0] for i in triangle ] for triangle in triangles ]
    assert len( triangles ) == len( set( [ frozenset( triangle ) for triangle in triangles ] ) )
    assert p.shape == u.shape
    
    kwargs.setdefault( 'verbose', False )
    kwargs.setdefault( 'xthresh', 1e-5 )
    kwargs.setdefault( 'fthresh', 1e-8 )
    kwargs.setdefault( 'max_iter', 100 )
    kwargs.setdefault( 'pressure_not_expansion', False )
    
    verbose = kwargs['verbose']
    
    XSTEP_THRESHOLD = kwargs['xthresh']
    F_THRESHOLD = kwargs['fthresh']
    MAX_ITERATIONS = kwargs['max_iter']
    
    ## We want to find env.a such that volume_expansion_factor = area_p / area_u.
    area_p = mesh_area( p, triangles )
    
    a_step = .01
    area_delta = 0.
    if kwargs['pressure_not_expansion']:
        env.p_m = env.p_atm + a_step
    else:
        env.a = a_step
    
    def geta():
        if kwargs['pressure_not_expansion']:
            return env.p_m
        else:
            return env.a
    def seta( val ):
        if kwargs['pressure_not_expansion']:
            env.p_m = val
        else:
            env.a = val
    
    u_n = u.copy()
    
    iteration = 0
    diverged = False
    while True:
        if verbose: print('-- volume iteration %s, a: %s --' % (iteration, geta()))
        iteration += 1
        
        if verbose: print('=> static_backwards_solution()')
        u_n_p_1 = static_backwards_solution( p, edges, u_n, triangles, constraints, env )[0]#, **kwargs )
        if verbose: print('<= static_backwards_solution()')
        
        #print 'would static_forwards_solution converge?:', static_forwards_solution( u_n_p_1, edges, u_n_p_1, triangles, constraints, env )[1]
        
        diff2 = sum( ( u_n_p_1 - u_n ) ** 2 )
        if verbose: print('| u_n+1 - u_n |^2:', diff2)
        u_n = u_n_p_1
        
        area_u_n = mesh_area( u_n, triangles )
        area_delta = area_p / area_u_n - volume_expansion_factor
        if verbose: print('| area_delta |:', area_delta)
        
        if diff2 < XSTEP_THRESHOLD: break
        if area_delta**2 < F_THRESHOLD: break
        
        if iteration >= MAX_ITERATIONS:
            print('Diverged.')
            diverged = True
            break
        
        if iteration > 1:
            ## Linear estimate for correct 'a' based on previous iteration
            ## The two points of the line are
            ## (a_n_m_1, area_delta_m_1) and (env.a, area_delta).
            a_n_p_1 = x_intercept( (a_n_m_1, area_delta_m_1), (geta(), area_delta) )
        else:
            a_n_p_1 = geta()
            if area_delta > 0.:
                a_n_p_1 -= a_step
            if area_delta < 0.:
                a_n_p_1 += a_step
        
        a_n_m_1 = geta()
        area_delta_m_1 = area_delta
        seta( a_n_p_1 )
    
    return u_n, geta(), not diverged

def test1():
    print('===== test1() =====')
    
    #debugger()
    
    import springs_setup
    u, edges, triangles, constraints = springs_setup.make_trigrid( 4,4 )
    
    print('u.shape:', u.shape)
    print('u:', u)
    
    env = Environment()
    
    print('## This should immediately return')
    env.a = 0.
    p_check = static_forwards_solution( u, edges, u, triangles, constraints, env )[0]
    print('| p_check - u |', mag( p_check.flatten() - u.flatten() ))
    print('| F( p_check, u ) |^2', mag2( F( p_check, edges, u, triangles, env ).flatten() ))
    
    print('''## Most basic recovery test:
    ## 1 Simulate forwards with some expansion force
    ## 2 Plug the solution as the initial configuration for the backwards solver (should terminate immediately).''')
    
    ## 1
    env.a = 1.
    p_exp = static_forwards_solution( u, edges, u, triangles, constraints, env, verbose = True )[0]
    print('p_exp:', p_exp)
    print('| p_exp - u |', mag( p_exp.flatten() - u.flatten() ))
    print('| F( p_exp, u ) |^2', mag2( F( p_exp, edges, u, triangles, env ).flatten() ))
    
    ## 2
    u_check = static_backwards_solution( p_exp, edges, u, triangles, constraints, env, verbose = True )[0]
    print('u_check:', u_check)
    print('| u_check - u |', mag( u_check.flatten() - u.flatten() ))
    print('| F( p_exp, u_check ) |^2', mag2( F( p_exp, edges, u_check, triangles, env ).flatten() ))
    
    
    print('## Now try to recover the backwards solution without a perfect initial guess')
    u_try = static_backwards_solution( p_exp, edges, p_exp, triangles, constraints, env, verbose = True )[0]
    print('u_try:', u_try)
    print('| u_try - u |', mag( u_try.flatten() - u.flatten() ))
    print('| F( p_exp, u_try ) |^2', mag2( F( p_exp, edges, u_try, triangles, env ).flatten() ))
    print('| u_try - u_check |', mag( u_try.flatten() - u_check.flatten() ))
    
    print('## Now try to recover a backwards solution which expands 0% (should return immediately)')
    u_no_expansion_solution, a_no_exp, converged = static_backwards_solution_with_given_volume_expansion( p_exp, edges, p_exp, triangles, constraints, env, 1., verbose = True )
    print('u_no_expansion_solution:', u_no_expansion_solution)
    print('a_no_exp:', a_no_exp)
    print('| u_no_expansion_solution - p_exp |', mag( u_no_expansion_solution.flatten() - p_exp.flatten() ))
    print('| F( p_exp, u_no_expansion_solution ) |^2', mag2( F( p_exp, edges, u_no_expansion_solution, triangles, env ).flatten() ))
    
    print('## Now try to recover a backwards solution which expands 50% (should iterate the outer loop)')
    u_50_expansion_solution, a_50_exp, converged = static_backwards_solution_with_given_volume_expansion( p_exp, edges, p_exp, triangles, constraints, env, 1.5, verbose = True )
    print('u_50_expansion_solution:', u_50_expansion_solution)
    print('a_50_exp:', a_50_exp)
    print('| u_50_expansion_solution - p_exp |', mag( u_50_expansion_solution.flatten() - p_exp.flatten() ))
    print('| F( p_exp, u_50_expansion_solution ) |^2', mag2( F( p_exp, edges, u_50_expansion_solution, triangles, env ).flatten() ))
    
    print('## Now simulate forwards with the given undeformed positions and expansion coefficient')
    env.a = a_50_exp
    p_50_exp = static_forwards_solution( u_50_expansion_solution, edges, u_50_expansion_solution, triangles, constraints, env, verbose = True )[0]
    print('p_50_exp:', p_50_exp)
    print('| p_50_exp - p_exp |', mag( p_50_exp.flatten() - p_exp.flatten() ))
    print('| F( p_50_exp, u_50_expansion_solution ) |^2', mag2( F( p_50_exp, edges, u_50_expansion_solution, triangles, env ).flatten() ))
    area_p_50_exp = mesh_area( p_50_exp, triangles )
    area_u_50_expansion_solution = mesh_area( u_50_expansion_solution, triangles )
    print('area( p_50_exp ):', area_p_50_exp)
    print('area( u_50_expansion_solution ):', area_u_50_expansion_solution)
    print('area( p_50_exp ) / area( u_50_expansion_solution ):', area_p_50_exp / area_u_50_expansion_solution)

def test2():
    print('===== test2() =====')
    
    #debugger()
    
    import springs_setup
    u, edges, triangles, constraints = springs_setup.make_trigrid( 4,4 )
    #u, edges, triangles, constraints = springs_setup.make_trigrid( 3,3 )
    print('constraints:', constraints)
    max_iter = 25
    
    print('u.shape:', u.shape)
    print('u:', u)
    
    env = Environment()
    env.kExpansion = False
    env.kPressure = True
    
    print('## Equilibrium check')
    print('## This should immediately return')
    env.p_m = 1.
    p_check = static_forwards_solution_in_steps( u, edges, u, triangles, constraints, env, max_iter = max_iter )[0]
    print('| p_check - u |', mag( p_check.flatten() - u.flatten() ))
    print('| F( p_check, u ) |^2', mag2( F( p_check, edges, u, triangles, env ).flatten() ))
    
    print('''## Most basic recovery test:
    ## 1 Simulate forwards with some expansion force
    ## 2 Plug the solution as the initial configuration for the backwards solver (should terminate immediately).''')
    
    ## 1
    expansion_p_m = 1.1
    #expansion_p_m = 1.7
    env.p_m = expansion_p_m
    p_exp = static_forwards_solution_in_steps( u, edges, u, triangles, constraints, env, max_iter = max_iter, verbose = True )[0]
    print('p_exp:', p_exp.round(2))
    print('| p_exp - u |', mag( p_exp.flatten() - u.flatten() ))
    print('| F( p_exp, u ) |^2', mag2( F( p_exp, edges, u, triangles, env ).flatten() ))
    
    ## 2
    u_check = static_backwards_solution( p_exp, edges, u, triangles, constraints, env, max_iter = max_iter, verbose = True )[0]
    print('u_check:', u_check.round(2))
    print('| u_check - u |', mag( u_check.flatten() - u.flatten() ))
    print('| F( p_exp, u_check ) |^2', mag2( F( p_exp, edges, u_check, triangles, env ).flatten() ))
    
    
    print('## Now try to recover the backwards solution without a perfect initial guess')
    u_try = static_backwards_solution( p_exp, edges, p_exp, triangles, constraints, env, max_iter = max_iter, verbose = True )[0]
    print('u_try:', u_try.round(2))
    print('| u_try - u |', mag( u_try.flatten() - u.flatten() ))
    print('| F( p_exp, u_try ) |^2', mag2( F( p_exp, edges, u_try, triangles, env ).flatten() ))
    print('| u_try - u_check |', mag( u_try.flatten() - u_check.flatten() ))
    
    print('## Now try to recover a backwards solution which expands 0% (should return immediately)')
    u_no_expansion_solution, a_no_exp, converged = static_backwards_solution_with_given_volume_expansion( p_exp, edges, p_exp, triangles, constraints, env, 1., verbose = True, pressure_not_expansion = True )
    print('u_no_expansion_solution:', u_no_expansion_solution.round(2))
    print('p_m_no_exp:', a_no_exp)
    print('| u_no_expansion_solution - p_exp |', mag( u_no_expansion_solution.flatten() - p_exp.flatten() ))
    print('| F( p_exp, u_no_expansion_solution ) |^2', mag2( F( p_exp, edges, u_no_expansion_solution, triangles, env ).flatten() ))
    
    
    print('## Now try to recover a backwards solution which expands the same amount as our forwards solution (should iterate the outer loop)')
    env.p_m = expansion_p_m
    volume_expansion = mesh_area( p_exp, triangles ) / mesh_area( u_check, triangles )
    print('mesh_area( p_exp ) / mesh_area( u_check ):', volume_expansion)
    u_N_expansion_solution, a_N_exp, converged = static_backwards_solution_with_given_volume_expansion( p_exp, edges, p_exp, triangles, constraints, env, volume_expansion, verbose = True, pressure_not_expansion = True )
    print('u_N_expansion_solution:', u_N_expansion_solution.round(2))
    print('p_m_N_exp:', a_N_exp)
    print('| u_N_expansion_solution - p_exp |', mag( u_N_expansion_solution.flatten() - p_exp.flatten() ))
    print('| u_N_expansion_solution - u |', mag( u_N_expansion_solution.flatten() - u.flatten() ))
    print('| F( p_exp, u_N_expansion_solution ) |^2', mag2( F( p_exp, edges, u_N_expansion_solution, triangles, env ).flatten() ))
    
    print('## Now simulate forwards with the given undeformed positions and pressure coefficient')
    env.p_m = a_N_exp
    p_N_exp = static_forwards_solution_in_steps( u_N_expansion_solution, edges, u_N_expansion_solution, triangles, constraints, env, max_iter = max_iter, verbose = True )[0]
    print('p_N_exp:', p_N_exp.round(2))
    print('| p_N_exp - p_exp |', mag( p_N_exp.flatten() - p_exp.flatten() ))
    print('| F( p_N_exp, u_N_expansion_solution ) |^2', mag2( F( p_N_exp, edges, u_N_expansion_solution, triangles, env ).flatten() ))
    area_p_N_exp = mesh_area( p_N_exp, triangles )
    area_u_N_expansion_solution = mesh_area( u_N_expansion_solution, triangles )
    print('area( p_N_exp ):', area_p_N_exp)
    print('area( u_N_expansion_solution ):', area_u_N_expansion_solution)
    print('area( p_N_exp ) / area( u_N_expansion_solution ):', area_p_N_exp / area_u_N_expansion_solution)
    
    
    print('## Now try to recover a backwards solution which expands 50% (should iterate the outer loop)')
    env.p_m = expansion_p_m
    volume_expansion = 1.5
    u_50_expansion_solution, a_50_exp, converged = static_backwards_solution_with_given_volume_expansion( p_exp, edges, p_exp, triangles, constraints, env, volume_expansion, verbose = True, pressure_not_expansion = True )
    print('u_50_expansion_solution:', u_50_expansion_solution.round(2))
    print('p_m_50_exp:', a_50_exp)
    print('| u_50_expansion_solution - p_exp |', mag( u_50_expansion_solution.flatten() - p_exp.flatten() ))
    print('| F( p_exp, u_50_expansion_solution ) |^2', mag2( F( p_exp, edges, u_50_expansion_solution, triangles, env ).flatten() ))
    
    print('## Now simulate forwards with the given undeformed positions and pressure coefficient')
    env.p_m = a_50_exp
    p_50_exp = static_forwards_solution_in_steps( u_50_expansion_solution, edges, u_50_expansion_solution, triangles, constraints, env, max_iter = max_iter, verbose = True )[0]
    print('p_50_exp:', p_50_exp.round(2))
    print('| p_50_exp - p_exp |', mag( p_50_exp.flatten() - p_exp.flatten() ))
    print('| F( p_50_exp, u_50_expansion_solution ) |^2', mag2( F( p_50_exp, edges, u_50_expansion_solution, triangles, env ).flatten() ))
    area_p_50_exp = mesh_area( p_50_exp, triangles )
    area_u_50_expansion_solution = mesh_area( u_50_expansion_solution, triangles )
    print('area( p_50_exp ):', area_p_50_exp)
    print('area( u_50_expansion_solution ):', area_u_50_expansion_solution)
    print('area( p_50_exp ) / area( u_50_expansion_solution ):', area_p_50_exp / area_u_50_expansion_solution)

def main():
    #test1()
    test2()

if __name__ == '__main__': main()
