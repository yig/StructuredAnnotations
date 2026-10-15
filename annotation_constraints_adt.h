#ifndef __annotation_constraints_adt_h__
#define __annotation_constraints_adt_h__

#ifdef __cplusplus
extern "C" {
#endif

typedef double real_t;

#ifdef _WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

DLLEXPORT
void same_xy_length_energy_gradient_hessian(
    /// These parameters are the annotation whose energy we are evaluating
    real_t* xform, // array of length 6
    real_t* annotation_axes, // array of length 9: axis[0].x, axis[0].y, axis[0].z, axis[1].x, ...
    int which_axis, // axis index (0,1, or 2)
    
    /// These parameters are the annotation whose length we want to match
    real_t* xform0, // array of length 6
    real_t* annotation0_axes, // array of length 9: axis[0].x, axis[0].y, axis[0].z, axis[1].x, ...
    int which_axis, // axis index (0,1, or 2)
    
    /// Output parameters
    real_t* energy_out, // 1 value
    real_t* gradient_out, // array of length 6
    real_t* hessian_out // row-major array of length 36
    );

DLLEXPORT
void same_xyz_length_energy_gradient_hessian(
    /// These parameters are the annotation whose energy we are evaluating
    real_t* xform, // array of length 6
    real_t* annotation_axes, // array of length 9: axis[0].x, axis[0].y, axis[0].z, axis[1].x, ...
    int which_axis, // axis index (0,1, or 2)
    
    /// These parameters are the annotation whose length we want to match
    real_t* xform0, // array of length 6
    real_t* annotation0_axes, // array of length 9: axis[0].x, axis[0].y, axis[0].z, axis[1].x, ...
    int which_axis, // axis index (0,1, or 2)
    
    /// Output parameters
    real_t* energy_out, // 1 value
    real_t* gradient_out, // array of length 6
    real_t* hessian_out // row-major array of length 36
    );

DLLEXPORT
void xform_energy_gradient_hessian(
    /// These parameters are the annotation whose energy we are evaluating
    real_t* xform, // array of length 6
    real_t* annotation_axes, // array of length 9: axis[0].x, axis[0].y, axis[0].z, axis[1].x, ...
    
    /// Output parameters
    real_t* energy_out, // 1 value
    real_t* gradient_out, // array of length 6
    real_t* hessian_out // row-major array of length 36
    );

#ifdef __cplusplus
}
#endif

#endif /* __annotation_constraints_adt_h__ */
