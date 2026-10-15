from numpy import *
from builtins import abs, bool, max, min, pow, round  ## NumPy 2 star-exports these; keep the builtins.

import annotations_sweep

class CuboidShape( annotations_sweep.SweepShape ):
    def __init__( self ):
        annotations_sweep.SweepShape.__init__( self )
        
        
        ## These are the defaults set by
        ## SweepShapeAnnotationCreateMode.SpineDrawer.doStrokeFinished():
        '''
        ## Add untilted circles at the first and last edges
        self.shape.tilt_radians[ 0. ] = 0.
        self.shape.tilt_radians[ 1. ] = 0.
        self.shape.radii[ 0. ] = GetSweepShapeCurrentRadius( self )
        self.shape.radii[ 1. ] = GetSweepShapeCurrentRadius( self )
        self.shape.symmetry_radians[ .5 ] = -.5*pi
        '''
        
        
        ## Make the default handles only one end;
        ## use t = 1. so they are close to the mouse.
        ## (Having default handles on both sides was always annoying
        ## with the cylinders, because I would often want to adjust
        ## the value along the entire cylinder, but I couldn't
        ## because there were two handles by default.)
        self.tilt_radians[ 1. ] = 0.
        ## Actually, I want symmetry handles on either end, since this lets
        ## me adjust the twist and the default is always correct
        ## (the edges of the cube are parallel to the view direction).
        self.symmetry_radians[ 0. ] = 0.
        self.symmetry_radians[ 1. ] = 0.
        
        ## We are creating a cuboid, so make the caps flat.
        ## TODO Q: Do I want the ends to be a little bit rounded,
        ##         like .1 * radius?
        ## A: Yes.  I do this in FromCurveAndRadius().
        self.cap_radii[0] = 0.
        self.cap_radii[1] = 0.
        
        ## Make the cross section a square.
        self.cross_sections[ 0. ] = annotations_sweep.normalize_arbitrary_cross_section(
                        default_cuboid_cross_section()
                        )
    
    def FromCurveAndRadius( curve, radius ):
        result = CuboidShape()
        
        result.radii[ 0. ] = radius
        #result.radii[ 1. ] = radius
        
        result.cap_radii[0] = .1 * radius
        result.cap_radii[1] = .1 * radius
        
        result.set_spine2d( curve )
        
        return result
    
    FromCurveAndRadius = staticmethod( FromCurveAndRadius )

def default_cuboid_cross_section():
    ## This doesn't work because cross sections must have
    ## annotations_sweep.kCrossSectionSamples samples.
    # return [ ( -1, -1 ), ( 1, -1 ), ( 1, 1 ), ( -1, 1 ) ]
    
    ## Project every point on the unit circle out onto the unit square [-1,1]x[-1,1].
    from annotations_sweep import unit_circle
    result = unit_circle / abs( unit_circle ).max( axis = 1 )[...,newaxis]
    
    ## Make it a little bit rounded.
    rounded_radius = .9 * sqrt(2)
    result_mags = sqrt( ( result**2 ).sum(1) )
    result_scale = ones( result_mags.size )
    result_scale[ result_mags > rounded_radius ] = rounded_radius / result_mags[ result_mags > rounded_radius ]
    result *= result_scale[...,newaxis]
    
    return result
