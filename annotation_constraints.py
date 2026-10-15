from math import *
from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.
from myarray import *



from helpers import Struct, solve_for_fixed_points
from vecutil import *

from annotations_ellipsoid import kEps

class AnnotationConstraintSolver( object ):
    def __init__( self, shape_annotations ):
        self.shape_annotations = shape_annotations
        
        ## Lists of the below classes.
        self.same_length_constraints_xy = []
        self.same_length_constraints_xyz = []
        #self.symmetry_constraints = []
        
        ## If an annotation is in this list, it won't change.
        self.constrained_annotations = []
        
        ## TODO: weights?
    
    def solve( self ):
        xforms = self.init_xforms()
        
        ## TODO: find xforms that minimize energy.
        
        ## Apply xform to each shape annotation
        for i in range(len( self.shape_annotations )):
            self.apply_xform_to_shape_annotation( xforms, i )
    
    def init_xforms( self ):
        '''
        'xforms' are the degrees-of-freedom we solve for.
        Each 'xform' is 6 numbers: scale for each axis, and axis/angle rotation represented as
        an axis vector whose magnitude is the rotation about the axis in radians.
        We initialize them to the current values for the shape annotations.
        '''
        
        xforms = zerosf( ( len( self.shape_annotations ), 6 ) )
        for i, ann in enumerate( self.shape_annotations ):
            ## TODO Q: Should I start these at 1, and scale the initial axes by this.
            ## A1: Yes.  That would make optimization more comparable.
            # xforms[i,:3] = sqrt( ( self.shape_annotations.ellipsoid.axes**2 ).sum(0) )
            xforms[i,:3] = 1.
            ## TODO Q: Decompose axes into an axis/angle rotation?
            ## A1: Axes are orthogonal, so the unit length axes are a rotation matrix, and all
            ##     rotations can be converted to a single axis angle rotation (proved by Euler).
            ## A2: Leave this as zeros, and apply it to the initial axes.
            xforms[i,3:] = 0.
        
        return xforms
    
    def apply_xform_to_shape_annotation( self, xforms, ann_index ):
        xform = xforms[ ann_index ]
        ann = self.shape_annotations[ ann_index ]
        
        ## Apply the scale.
        ann.ellipsoid.axes *= xform[:3]
        
        ## Apply the rotation.
        radians2 = mag2( xform[3:] )
        if radians2 > kEps:
            radians = sqrt( radians2 )
            axis = xform[3:] / radians
            R = get_rotation_matrix( axis, radians )
            ann.ellipsoid.axes[:,:] = R @ ann.ellipsoid.axes
    
    def xform_mag_xy( self, xforms, ann_index, axis_index ):
        orig_axis_xy = self.shape_annotations[ ann_index ].ellipsoid.axes[ axis_index, :2 ]
        return mag( xforms[ ann_index, axis_index ] * orig_axis_xy )
    
    def xform_mag_xyz( self, xforms, ann_index, axis_index ):
        orig_axis_xyz = self.shape_annotations[ ann_index ].ellipsoid.axes[ axis_index ]
        return mag( xforms[ ann_index, axis_index ] * orig_axis_xyz )
    
    def energy( self, xforms ):
        assert len( self.shape_annotations ) == len( xforms )
        
        E = 0.
        
        ## xy length constraints contribute [ xy_length - constraint0.xy_length ]^2 to the energy.
        if self.same_length_constraints_xy:
            c0 = self.same_length_constraints_xy[0]
            c0ann_index = self.shape_annotations.index( c0[0] )
            c0mag_xy = self.xform_mag_xy( xforms, c0ann_index, c0[1] )
            for c in self.same_length_constraints_xy[1:]:
                for ann, axis_index in c.annotations_and_axis_indices:
                    ann_index = self.shape_annotations.index( ann )
                    E += ( self.xform_mag_xy( xforms, ann_index, axis_index ) - c0mag_xy )**2
        
        ## xyz length constraints contribute [ xyz_length - constraint0.xyz_length ]^2 to the energy.
        if self.same_length_constraints_xyz:
            c0 = self.same_length_constraints_xyz[0]
            c0ann_index = self.shape_annotations.index( c0[0] )
            c0mag_xyz = self.xform_mag_xyz( xforms, c0ann_index, c0[1] )
            for c in self.same_length_constraints_xyz[1:]:
                for ann, axis_index in c.annotations_and_axis_indices:
                    ann_index = self.shape_annotations.index( ann )
                    E += ( self.xform_mag_xyz( xforms, ann_index, axis_index ) - c0mag_xyz )**2
        
        ## Regularization, assume the user's input is something we want to keep.
        for ann_index, ann in enumerate( self.shape_annotations ):
            E_ann = 0.
            '''
            ## Begin with the magnitude of the scale - magnitude of the original scale
            axes_mags = sqrt( ( ann.ellipsoid.axes**2 ).sum(1) )
            E_ann += ( ( axes_mags - ann_original_magnitude_scale )**2 ).sum()
            ## Also penalize the rotation.  Take the magnitude difference from the identity matrix.
            E_ann += ( ( ann.ellipsoid.axes / axes_mags - ann_original_axes_normalized )**2 ).sum()
            '''
            ## TODO Q: Why not just take the difference squared from the original axes.
            ## A1: This is a good idea.  I care about rotation and scale changes.
            #E_ann += ( ( ann.ellipsoid.axes - ann_original_axes )**2 ).sum()
            ## A2: Let's take the magnitude of the scale minus 1, and the magnitude of the rotation.
            E_ann += mag2( xforms[ ann_index, :3 ] - 1. )
            E_ann += mag2( xforms[ ann_index, 3: ] )
            
            w_ann = 1e7 if ann in self.constrained_annotations else 1.
            E += w_ann * E_ann
        
        return E
    
    def gradient( self ):
        ## TODO: This (autodiff)
        ## TODO Q: What are my dofs? for each ellipsoid: center, axis0, axis1 length + angle, axis2 length?
        ## A1: For each ellipsoid, dofs are scale for each axis, and a rotation represented as an axis (phi,theta) and angle.
        pass

class SameLengthConstraint( object ):
    def __init__( self, annotations_and_axis_indices ):
        self.annotations_and_axis_indices = annotations_and_axis_indices
        
        assert len( self.annotations_and_axis_indices ) > 0
        assert None not in [ ann for (ann, axis_index) in self.annotations_and_axis_indices ]
        assert False not in [ axis_index >= 0 and axis_index < 3 for (ann, axis_index) in self.annotations_and_axis_indices ]

class SameLengthConstraintXY( SameLengthConstraint ):
    def __init__( self, annotations_and_axis_indices ):
        SameLengthConstraint.__init__( annotations_and_axis_indices )

class SameLengthConstraintXYZ( SameLengthConstraint ):
    def __init__( self, annotations_and_axis_indices ):
        SameLengthConstraint.__init__( annotations_and_axis_indices )
