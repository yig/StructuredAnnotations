#if 0
make -f Makefile.fibermesh

MESHOPT_DIR=/Users/yotam/Work/meshopt/git-devel/meshopt ;\
EXT_DIR=${MESHOPT_DIR}/ext ;\
mpic++ fibermesh.cpp \
    -dynamiclib -o libfibermesh.dylib \
    -g -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare \
    -fno-inline \
    -I${MESHOPT_DIR}/src -I${EXT_DIR} \
    -I${EXT_DIR}/petsc-2.3 \
    -I${EXT_DIR}/petsc-2.3/bmake/darwin8.8.1-cxx-opt \
    -I${EXT_DIR}/petsc-2.3/include \
    ${MESHOPT_DIR}/libhiorder.a \
    ${EXT_DIR}/OpenMesh/Core/Darwin_gcc_opt/libOpenMesh_Core.dylib \
    ${EXT_DIR}/UMFPACK/Lib/libumfpack.a ${EXT_DIR}/AMD/Lib/libamd.a \
    \
    -L${EXT_DIR}/petsc-2.3/lib/darwin8.8.1-cxx-opt \
    -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc \
    -L/usr/X11R6/lib -lX11 \
    /System/Library/Frameworks/vecLib.framework/vecLib
    
    -fno-inline \
    -O1 \
    -L${EXT_DIR}/petsc-2.3/externalpackages/mpich2-1.0.4p1/darwin8.8.1-cxx-opt/lib \
    -lmpich -lpmpich -lmpichcxx

mpic++ instead of g++
#endif

#include <cassert>
#include <algorithm>
#include <string>
#include <cstdio>
#include "optimizers/ConstrainedSystem.h"
#include "common/SparseMatrix_PETSc.h"
#include "forces/LaplacianEditing.h" // LaplacianEditing::cotangentWeights() for computing area
#include "common/bobbles_mesh.h" // get_one_ring*()
#include "common/stl.h" // tostr()
#include <petscmat.h>
#include "common/Timing.hpp"

//#define FN_ENTER do { std::cerr << "Entering function " << __PRETTY_FUNCTION__ << "\n"; } while(0);
#define FN_ENTER ;
//#define FN_ENTER Tick t( __PRETTY_FUNCTION__ );
//#define FN_ENTER Tick t( __func__ );
// #define FN_EXIT do { std::cerr << "Exiting function " << __PRETTY_FUNCTION__ << "\n"; } while(0);
#define FN_EXIT ;

typedef SparseMatrix::value_t real_t;
typedef SparseMatrix::index_t int_t;

// This is so ugly.  I get a link error if I call these directly.  I think it has to
// do with the MPI object parameter.
// It goes away if I use a proper Makefile that includes PETSc's Makefile bit.
// These wrappers are hidden in meshopt/src/optimizers/build_constrained_system.cpp
// PetscErrorCode MatCreateSeqDense_wrapper( PetscInt m, PetscInt n, PetscScalar *data, Mat *A );
// PetscErrorCode VecCreateSeqWithArray_wrapper( PetscInt n, const PetscScalar array[], Vec *V );

extern "C"
{
struct optimize_surface2_params_t
{
    // The number of iterations to perform, unless 'absolute_convergence_threshold' is reached.
    int_t max_iterations;
    // If the maximum difference between two vertices is less than 'absolute_convergence_threshold',
    // stop iterating.
    real_t absolute_convergence_threshold;
    // Whether to use graph laplacian weights or cotangent weights.
    int_t graph_not_cotangent_weights;
    // Whether to impose average edge lengths (necessary for graph_not_cotangent_weights).
    int_t average_edge_lengths;
    // Whether to restrict vertices to only move in the normal direction (prevents tangential drift)
    int_t vertex_normal_motion_only;
};
}


namespace fibermesh
{
bool kClipNegativeCurvatureConstraints = true;
bool kPrintConstraints = false;

// used
void compute_target_curvatures_mass_laplacian(
    const TriMesh& mesh, const SparseMatrix_PETSc& L, const SparseMatrix_PETSc& L2,
    const ConstrainedSystem::Indices& constraint_vertices, ConstrainedSystem::Values constraint_curvatures,
    const real_t constraint_weight,
    ConstrainedSystem::Values& LMs_out,
    ConstrainedSystem** system_handle = NULL
    )
{FN_ENTER
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    assert( N == L.GetNumColumns() );
    assert( L.GetNumRows() == L.GetNumColumns() );
    assert( L.GetNumRows() == L2.GetNumRows() );
    assert( L.GetNumColumns() == L2.GetNumColumns() );
    
    assert( constraint_weight >= 0. );
    
    assert( constraint_vertices.size() == constraint_curvatures.size() );
    
#ifndef NDEBUG
    for( int i = 0; i < constraint_vertices.size(); ++i )
        assert( constraint_vertices.at(i) >= 0 && constraint_vertices.at(i) < N );
#endif
    
    std::cout << "negative constraint_curvatures: ";
    for( int i = 0; i < constraint_curvatures.size(); ++i )
    {
        if( constraint_curvatures.at(i) < 0. )
        {
            std::cout << constraint_curvatures.at(i) << ", ";
            if( kClipNegativeCurvatureConstraints ) constraint_curvatures.at(i) = 0.;
        }
    }
    std::cout << '\n';
    
    LMs_out.resize( N );
    ConstrainedSystem::Values b( N, 0. );
    ConstrainedSystem* system(0);
    if( system_handle && *system_handle )
    {
        std::cout << "Re-using existing curvature system\n";
        system = *system_handle;
        system->UpdateSystem( L2, b );
    }
    else
    {
        std::cout << "Creating a new curvature system\n";
        system = new ConstrainedSystem( L2, b );
        if( system_handle ) *system_handle = system;
    }
    assert( system );
    // ConstrainedSystem system( L2, b );
    
    if( system->GetLSQConstraintsWeight() != constraint_weight )
    {
        std::cout << "Replacing curvature constraint weight\n";
        system->SetLSQConstraintsWeight( constraint_weight );
    }
    
    // Since our constraint equations are all very simple, a single vertex dof = some value,
    // we'll use the vertex dof as our equation IDs.
    // If the set of equation IDs differs from the set of constrained vertex dofs, we'll
    // simply clear them all and set them anew.  This isn't so inefficient since any change in the
    // equations (not the rhs's) necessitates rebuilding the constraint matrix anyhow.
    std::vector< int > all_eq_ids;
    system->all_eq_ids( all_eq_ids );
    if( std::set<int>( all_eq_ids.begin(), all_eq_ids.end() ) != std::set<int>( constraint_vertices.begin(), constraint_vertices.end() ) )
    {
        std::cout << "Replacing curvature constraints\n";
        for( int i = 0; i < all_eq_ids.size(); ++i ) system->remove_eq( all_eq_ids.at(i) );
        
        ConstrainedSystem::Indices eq_indices( 1, 0 );
        ConstrainedSystem::Values eq_vals( 1, 1. );
        for( int i = 0; i < constraint_vertices.size(); ++i )
        {
            eq_indices.at(0) = constraint_vertices.at(i);
            const int id = system->add_eq(
                eq_indices, eq_vals, constraint_curvatures.at(i),
                constraint_vertices.at(i)
                );
            assert( -1 != id );
        }
    }
    // But if the set of equation IDs is the same as the set of constrained vertex dofs, we still
    // need to update the right-hand-side values, since these change every time.
    else
    {
        std::cout << "Replacing curvature constraint values\n";
        for( int i = 0; i < constraint_vertices.size(); ++i )
        {
            system->update_eq( constraint_vertices.at(i), constraint_curvatures.at(i) );
        }
    }
    system->solve( &LMs_out[0] );
    
    std::cout << "negative LMs: ";
    for( int i = 0; i < LMs_out.size(); ++i )
    {
        if( LMs_out.at(i) < 0. )
        {
            std::cout << LMs_out.at(i) << ", ";
            LMs_out.at(i) = 0.;
        }
    }
    std::cout << '\n';
FN_EXIT}

namespace
{
template< typename forward_iter >
void assert_indices_inbounds( const forward_iter begin, const forward_iter end, const int lower, const int upper )
{
    for( forward_iter iter = begin; iter != end; ++iter )
    {
        const int val = *iter;
        assert( val >= lower && val <= upper );
    }
}
#if 0
template< typename T >
void assert_indices_inbounds( const std::vector< T >& vec, const int lower, const int upper )
{
    for( typename std::vector<T>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter )
    {
        const int val = *iter;
        assert( val >= lower && val <= upper );
    }
}
#endif

void print_hard_constraints(
    const ConstrainedSystem::Indices& hard_constraints_indices,
    const ConstrainedSystem::Values& hard_constraints_values
    )
{
    std::cout << " ([";
    for( int i = 0; i < hard_constraints_indices.size(); ++i )
    {
        std::cout
            << "("
            << hard_constraints_indices.at(i) << ", "
            << hard_constraints_values.at(i)
            << ")";
        if( i < hard_constraints_indices.size()-1 ) std::cout << ",\n   ";
    }
    std::cout << "])\n";
}

void print_linear_constraints(
    const ConstrainedSystem::Indices& linear_constraints_indices,
    const ConstrainedSystem::Values& linear_constraints_values,
    const real_t& linear_constraints_rhs
    )
{
    std::cout << "  ([";
    for( int i = 0; i < linear_constraints_indices.size(); ++i )
    {
        std::cout << linear_constraints_indices.at(i);
        if( i < linear_constraints_indices.size()-1 ) std::cout << ", ";
    }
    std::cout << "],\n";
    std::cout << "   [";
    for( int i = 0; i < linear_constraints_values.size(); ++i )
    {
        std::cout << linear_constraints_values.at(i);
        if( i < linear_constraints_values.size()-1 ) std::cout << ",\n    ";
    }
    std::cout << "],\n   ";
    std::cout << linear_constraints_rhs << ")\n";
}

void print_linear_constraints(
    const ConstrainedSystem::Indices& linear_constraints_indices,
    const ConstrainedSystem::Indices& linear_constraints_coords,
    const ConstrainedSystem::Values& linear_constraints_values,
    const real_t& linear_constraints_rhs
    )
{
    std::cout << "  ([";
    for( int i = 0; i < linear_constraints_indices.size(); ++i )
    {
        std::cout
            << "( "
            << linear_constraints_indices.at(i) << ", "
            << linear_constraints_coords.at(i)
            << " )";
        if( i < linear_constraints_indices.size()-1 ) std::cout << ", ";
    }
    std::cout << "],\n";
    std::cout << "   [";
    for( int i = 0; i < linear_constraints_values.size(); ++i )
    {
        std::cout << linear_constraints_values.at(i);
        if( i < linear_constraints_values.size()-1 ) std::cout << ",\n    ";
    }
    std::cout << "],\n   ";
    std::cout << linear_constraints_rhs << ")\n";
}

// used
void solve_uncoupled_system(
    SparseMatrix_PETSc& A,
    SparseMatrix::value_t* b,
    const SparseMatrix_PETSc* Afresh,
    SparseMatrix::value_t* solution_out,
    
    const ConstrainedSystem::Indices& hard_constraints_indices,
    const ConstrainedSystem::Values& hard_constraints_values,
    
    const std::vector< ConstrainedSystem::Indices >& linear_constraints_indices,
    const std::vector< ConstrainedSystem::Values >& linear_constraints_values,
    const ConstrainedSystem::Values& linear_constraints_rhss,
    
    real_t linear_constraint_weight,
    
    ConstrainedSystem** system_handle = NULL
    )
{FN_ENTER
    for( int i = 0; i < hard_constraints_indices.size(); ++i )
    {
        const int row = hard_constraints_indices.at(i);
        A.IdentityRow( row );
        b[ row ] = hard_constraints_values.at(i);
    }
    A.Finalize();
    
    // Without linear constraints, the overhead of using ConstrainedSystem is a straight copy
    // of the data.
    // At the same time, ConstrainedSystem is set up to call the fastest solver I have.
    ConstrainedSystem* system(0);
    bool new_system = true;
    if( system_handle && *system_handle )
    {
        std::cout << "Re-using existing position system\n";
        system = *system_handle;
        new_system = false;
        if( std::set<int>( hard_constraints_indices.begin(), hard_constraints_indices.end() ) == std::set<int>( system->GetRowMask().begin(), system->GetRowMask().end() ) )
        {
            system->UpdateSystem( A, b );
        }
        else
        {
            system->UpdateSystem( A, b, hard_constraints_indices );
        }
    }
    else
    {
        std::cout << "Creating a new position system\n";
        system = new ConstrainedSystem( A, b, hard_constraints_indices );
        if( system_handle ) *system_handle = system;
    }
    assert( system );
    // ConstrainedSystem system( A, b, hard_constraints_indices );
    // system.SetLSQConstraintsWeight( linear_constraint_weight );
    if( system->GetLSQConstraintsWeight() != linear_constraint_weight )
    {
        std::cout << "Replacing position constraint weight\n";
        system->SetLSQConstraintsWeight( linear_constraint_weight );
    }
    
    // We are assuming that neither constraint equations nor the right-hand-side values ever change.
    if( new_system )
    {
        for( int i = 0; i < linear_constraints_indices.size(); ++i )
        {
            system->add_eq(
                linear_constraints_indices.at(i),
                linear_constraints_values.at(i),
                linear_constraints_rhss.at(i)
                );
        }
    }
    else
    {
        std::cout << "Assuming position linear constraints don't change every iteration.\n";
    }
    
    system->solve( solution_out );
    
    if( Afresh )
    {
        assert( Afresh->GetNumRows() == A.GetNumRows() );
        assert( Afresh->GetNumColumns() == A.GetNumColumns() );
        
        ConstrainedSystem::Indices indices;
        ConstrainedSystem::Values values;
        for( int i = 0; i < hard_constraints_indices.size(); ++i )
        {
            const int row = hard_constraints_indices.at(i);
            Afresh->GetRowValues( row, indices, values );
            // The matrix needs to be finalized before calling ZeroRow().
            A.Finalize();
            A.ZeroRow( row );
            A.SetRowValues( row, indices, values );
        }
        A.Finalize();
    }
FN_EXIT}
}

// used
void optimize_positions_for_mass_laplacians_complex_constraints(
    TriMesh& mesh,
    const SparseMatrix_PETSc& L, const SparseMatrix_PETSc& L2,
    const ConstrainedSystem::Values& LMs,
    
    const ConstrainedSystem::Indices& x_hard_constraints_indices,
    const ConstrainedSystem::Values& x_hard_constraints_values,
    
    const std::vector< ConstrainedSystem::Indices >& x_linear_constraints_indices,
    const std::vector< ConstrainedSystem::Values >& x_linear_constraints_coeffs,
    const ConstrainedSystem::Values& x_linear_constraints_rhss,
    
    
    const ConstrainedSystem::Indices& y_hard_constraints_indices,
    const ConstrainedSystem::Values& y_hard_constraints_values,
    
    const std::vector< ConstrainedSystem::Indices >& y_linear_constraints_indices,
    const std::vector< ConstrainedSystem::Values >& y_linear_constraints_coeffs,
    const ConstrainedSystem::Values& y_linear_constraints_rhss,
    
    
    const ConstrainedSystem::Indices& z_hard_constraints_indices,
    const ConstrainedSystem::Values& z_hard_constraints_values,
    
    const std::vector< ConstrainedSystem::Indices >& z_linear_constraints_indices,
    const std::vector< ConstrainedSystem::Values >& z_linear_constraints_coeffs,
    const ConstrainedSystem::Values& z_linear_constraints_rhss,
    
    
    const std::vector< ConstrainedSystem::Indices >& general_linear_constraints_vertex_indices,
    const std::vector< ConstrainedSystem::Indices >& general_linear_constraints_coordinate_offsets,
    const std::vector< ConstrainedSystem::Values >& general_linear_constraints_coeffs,
    const ConstrainedSystem::Values& general_linear_constraints_rhss,
    
    
    const optimize_surface2_params_t& params,
    
    const real_t constraint_weight = 1.,
    
    ConstrainedSystem** system_handle = NULL
    )
{FN_ENTER
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    assert( N == L.GetNumColumns() );
    assert( L.GetNumRows() == L.GetNumColumns() );
    assert( L.GetNumRows() == L2.GetNumRows() );
    assert( L.GetNumColumns() == L2.GetNumColumns() );
    assert( N == LMs.size() );
    
    assert( constraint_weight >= 0. );
    
    assert( x_hard_constraints_indices.size() == x_hard_constraints_values.size() );
    assert( y_hard_constraints_indices.size() == y_hard_constraints_values.size() );
    assert( z_hard_constraints_indices.size() == z_hard_constraints_values.size() );
    
    assert( x_linear_constraints_indices.size() == x_linear_constraints_coeffs.size() );
    assert( x_linear_constraints_indices.size() == x_linear_constraints_rhss.size() );
    
    assert( y_linear_constraints_indices.size() == y_linear_constraints_coeffs.size() );
    assert( y_linear_constraints_indices.size() == y_linear_constraints_rhss.size() );
    
    assert( z_linear_constraints_indices.size() == z_linear_constraints_coeffs.size() );
    assert( z_linear_constraints_indices.size() == z_linear_constraints_rhss.size() );
    
    assert( general_linear_constraints_vertex_indices.size() == general_linear_constraints_coordinate_offsets.size() );
    assert( general_linear_constraints_vertex_indices.size() == general_linear_constraints_coeffs.size() );
    assert( general_linear_constraints_vertex_indices.size() == general_linear_constraints_rhss.size() );
    
#ifndef NDEBUG
    for( int i = 0; i < x_linear_constraints_indices.size(); ++i )
    {
        assert( x_linear_constraints_indices.at(i).size() == x_linear_constraints_coeffs.at(i).size() );
        assert_indices_inbounds( x_linear_constraints_indices.at(i).begin(), x_linear_constraints_indices.at(i).end(), 0, N-1 );
    }
    for( int i = 0; i < y_linear_constraints_indices.size(); ++i )
    {
        assert( y_linear_constraints_indices.at(i).size() == y_linear_constraints_coeffs.at(i).size() );
        assert_indices_inbounds( y_linear_constraints_indices.at(i).begin(), y_linear_constraints_indices.at(i).end(), 0, N-1 );
    }
    for( int i = 0; i < z_linear_constraints_indices.size(); ++i )
    {
        assert( z_linear_constraints_indices.at(i).size() == z_linear_constraints_coeffs.at(i).size() );
        assert_indices_inbounds( z_linear_constraints_indices.at(i).begin(), z_linear_constraints_indices.at(i).end(), 0, N-1 );
    }
    for( int i = 0; i < general_linear_constraints_vertex_indices.size(); ++i )
    {
        assert( general_linear_constraints_vertex_indices.at(i).size() == general_linear_constraints_coordinate_offsets.at(i).size() );
        assert( general_linear_constraints_vertex_indices.at(i).size() == general_linear_constraints_coeffs.at(i).size() );
        assert_indices_inbounds( general_linear_constraints_vertex_indices.at(i).begin(), general_linear_constraints_vertex_indices.at(i).end(), 0, N-1 );
        assert_indices_inbounds( general_linear_constraints_coordinate_offsets.at(i).begin(), general_linear_constraints_coordinate_offsets.at(i).end(), 0, 2 );
    }
#endif
    
    ConstrainedSystem::Values solution( N*3 );
    ConstrainedSystem::Values b( N*3 );
    {
        for( int i = 0; i < N; ++i )
        {
            const TriMesh::Normal& normal = mesh.normal( TriMesh::VertexHandle( i ) );
            b[ i + 0*N ] = normal[0];
            b[ i + 1*N ] = normal[1];
            b[ i + 2*N ] = normal[2];
        }
        
        int ierr(0);
#if 1 // For some reason I get a link error when I call MatCreateSeqDense().
        // Create a dense PETSc matrix containing the normals.
        ConstrainedSystem::Values Ln( N*3 );
        Mat normals;
        ierr = MatCreateSeqDense( PETSC_COMM_WORLD, N, 3, &b[0], &normals ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        Mat Lnmat;
        ierr = MatCreateSeqDense( PETSC_COMM_WORLD, N, 3, &Ln[0], &Lnmat ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        
        // PETSc aborts if this hasn't been called.
        ierr = MatAssemblyBegin( normals, MAT_FINAL_ASSEMBLY ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        ierr = MatAssemblyEnd( normals, MAT_FINAL_ASSEMBLY ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        ierr = MatAssemblyBegin( Lnmat, MAT_FINAL_ASSEMBLY ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        ierr = MatAssemblyEnd( Lnmat, MAT_FINAL_ASSEMBLY ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        
        ierr = MatMatMult( L.PETScMat(), normals, MAT_REUSE_MATRIX, 1., &Lnmat ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        // Destroy the matrix.
        MatDestroy( normals );
        MatDestroy( Lnmat );
        // Replace b with Ln
        std::copy( Ln.begin(), Ln.end(), b.begin() );
#elif 0 // For some reason I get a link error when I call VecCreateSeqWithArray().
        Vec normalcoord;
        // x coordinate
        
        this is broken, because MatMult() cannot work in-place. in other words, we cannot do x = A*x
        // ConstrainedSystem::Values normalcoorddata( N );
        // std::copy( b.begin() + 0*N, b.begin() + 1*N, normalcoorddata.begin() );
        
        ierr = VecCreateSeqWithArray( PETSC_COMM_WORLD, N, &b[0*N], &normalcoord );
        ierr = MatMult( L.PETScMat(), normalcoord, normalcoord ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        VecDestroy( normalcoord );
        // y coordinate
        // ierr = VecCreateSeqWithArray_wrapper( N, &b[1*N], &normalcoord );
        ierr = VecCreateSeqWithArray( PETSC_COMM_WORLD, N, &b[1*N], &normalcoord );
        ierr = MatMult( L.PETScMat(), normalcoord, normalcoord ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        VecDestroy( normalcoord );
        // z coordinate
        // ierr = VecCreateSeqWithArray_wrapper( N, &b[2*N], &normalcoord );
        ierr = VecCreateSeqWithArray( PETSC_COMM_WORLD, N, &b[2*N], &normalcoord );
        ierr = MatMult( L.PETScMat(), normalcoord, normalcoord ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        VecDestroy( normalcoord );
#elif 1 // What?  Not even VecCreate() works?
        Vec Con_rhs_vec;
        ierr = VecCreate( PETSC_COMM_WORLD, &Con_rhs_vec ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
#endif
        // Scale each normal by LM.
        for( int i = 0; i < N; ++i )
        {
            b[ i + 0*N ] *= LMs[i];
            b[ i + 1*N ] *= LMs[i];
            b[ i + 2*N ] *= LMs[i];
        }
    }
    
    if( kPrintConstraints )
    {
        print_hard_constraints( x_hard_constraints_indices, x_hard_constraints_values );
        for( int i = 0; i < x_linear_constraints_indices.size(); ++i )
            print_linear_constraints(
                x_linear_constraints_indices.at(i),
                x_linear_constraints_coeffs.at(i),
                x_linear_constraints_rhss.at(i)
                );
        print_hard_constraints( y_hard_constraints_indices, y_hard_constraints_values );
        for( int i = 0; i < y_linear_constraints_indices.size(); ++i )
            print_linear_constraints(
                y_linear_constraints_indices.at(i),
                y_linear_constraints_coeffs.at(i),
                y_linear_constraints_rhss.at(i)
                );
        print_hard_constraints( z_hard_constraints_indices, z_hard_constraints_values );
        for( int i = 0; i < z_linear_constraints_indices.size(); ++i )
            print_linear_constraints(
                z_linear_constraints_indices.at(i),
                z_linear_constraints_coeffs.at(i),
                z_linear_constraints_rhss.at(i)
                );
        
        for( int i = 0; i < general_linear_constraints_vertex_indices.size(); ++i )
            print_linear_constraints(
                general_linear_constraints_vertex_indices.at(i),
                general_linear_constraints_coordinate_offsets.at(i),
                general_linear_constraints_coeffs.at(i),
                general_linear_constraints_rhss.at(i)
                );
    }
    
    // NOTE: both code paths have been tested
    if( general_linear_constraints_vertex_indices.size() == 0 )
    {
        SparseMatrix_PETSc A( L2 );
        
        // x
        solve_uncoupled_system(
            A, &b[ N*0 ], &L2, &solution[ N*0 ],
            x_hard_constraints_indices,
            x_hard_constraints_values,
            x_linear_constraints_indices,
            x_linear_constraints_coeffs,
            x_linear_constraints_rhss,
            constraint_weight
            );
        // y
        solve_uncoupled_system(
            A, &b[ N*1 ], &L2, &solution[ N*1 ],
            y_hard_constraints_indices,
            y_hard_constraints_values,
            y_linear_constraints_indices,
            y_linear_constraints_coeffs,
            y_linear_constraints_rhss,
            constraint_weight
            );
        // z
        solve_uncoupled_system(
            A, &b[ N*2 ], 0, &solution[ N*2 ],
            z_hard_constraints_indices,
            z_hard_constraints_values,
            z_linear_constraints_indices,
            z_linear_constraints_coeffs,
            z_linear_constraints_rhss,
            constraint_weight
            );
    }
    else
    {
        // Build a block diagonal matrix of size 3*N with 3 copies of L2
        std::vector< int > nnzs( 3*N, 0 );
        for( int i = 0; i < N; ++i )
        {
            const int nz = L2.GetRowNNZ( i );
            nnzs.at( i + 0*N ) = nz;
            nnzs.at( i + 1*N ) = nz;
            nnzs.at( i + 2*N ) = nz;
        }
        SparseMatrix_PETSc A( 3*N, 3*N, nnzs );
        {
            ConstrainedSystem::Indices indices;
            ConstrainedSystem::Values vals;
            ConstrainedSystem::Indices::iterator it;
            for( int i = 0; i < N; ++i )
            {
                L2.GetRowValues( i, indices, vals );
                // The first copy goes in the usual place.
                A.SetRowValues( i + 0*N, indices, vals );
                // The second copy has indices offset by N.
                for( it = indices.begin(); it != indices.end(); ++it ) *it += N;
                A.SetRowValues( i + 1*N, indices, vals );
                // The third copy has indices offset by N again.
                for( it = indices.begin(); it != indices.end(); ++it ) *it += N;
                A.SetRowValues( i + 2*N, indices, vals );
            }
        }
        A.Finalize();
        
        
        // Now build the list of all constraints.
        ConstrainedSystem::Indices hard_constraints_indices;
        ConstrainedSystem::Values hard_constraints_values;
        std::vector< ConstrainedSystem::Indices > linear_constraints_indices;
        std::vector< ConstrainedSystem::Values > linear_constraints_coeffs;
        ConstrainedSystem::Values linear_constraints_rhss;
        
        // Concatenate all hard constraints.
        const int num_hard_constraints =
            x_hard_constraints_indices.size() +
            y_hard_constraints_indices.size() +
            z_hard_constraints_indices.size()
            ;
        hard_constraints_indices.reserve( num_hard_constraints );
        hard_constraints_values.reserve( num_hard_constraints );
        for( int i = 0; i < x_hard_constraints_indices.size(); ++i )
        {
            hard_constraints_indices.push_back( 0*N + x_hard_constraints_indices.at(i) );
            hard_constraints_values.push_back( x_hard_constraints_values.at(i) );
        }
        for( int i = 0; i < y_hard_constraints_indices.size(); ++i )
        {
            hard_constraints_indices.push_back( 1*N + y_hard_constraints_indices.at(i) );
            hard_constraints_values.push_back( y_hard_constraints_values.at(i) );
        }
        for( int i = 0; i < z_hard_constraints_indices.size(); ++i )
        {
            hard_constraints_indices.push_back( 2*N + z_hard_constraints_indices.at(i) );
            hard_constraints_values.push_back( z_hard_constraints_values.at(i) );
        }
        assert( hard_constraints_indices.size() == hard_constraints_values.size() );
        assert( hard_constraints_indices.size() == num_hard_constraints );
        
        // Concatenate all linear constraints.
        const int num_linear_constraints =
            x_linear_constraints_indices.size() +
            y_linear_constraints_indices.size() +
            z_linear_constraints_indices.size() +
            general_linear_constraints_vertex_indices.size()
            ;
        linear_constraints_indices.reserve( num_linear_constraints );
        linear_constraints_coeffs.reserve( num_linear_constraints );
        linear_constraints_rhss.reserve( num_linear_constraints );
        ConstrainedSystem::Indices::iterator it;
        for( int i = 0; i < x_linear_constraints_indices.size(); ++i )
        {
            linear_constraints_indices.push_back( x_linear_constraints_indices.at(i) );
            linear_constraints_coeffs.push_back( x_linear_constraints_coeffs.at(i) );
            linear_constraints_rhss.push_back( x_linear_constraints_rhss.at(i) );
        }
        for( int i = 0; i < y_linear_constraints_indices.size(); ++i )
        {
            linear_constraints_indices.push_back( y_linear_constraints_indices.at(i) );
            for( it = linear_constraints_indices.back().begin(); it != linear_constraints_indices.back().end(); ++it ) *it += N;
            linear_constraints_coeffs.push_back( y_linear_constraints_coeffs.at(i) );
            linear_constraints_rhss.push_back( y_linear_constraints_rhss.at(i) );
        }
        for( int i = 0; i < z_linear_constraints_indices.size(); ++i )
        {
            linear_constraints_indices.push_back( z_linear_constraints_indices.at(i) );
            for( it = linear_constraints_indices.back().begin(); it != linear_constraints_indices.back().end(); ++it ) *it += 2*N;
            linear_constraints_coeffs.push_back( z_linear_constraints_coeffs.at(i) );
            linear_constraints_rhss.push_back( z_linear_constraints_rhss.at(i) );
        }
        for( int i = 0; i < general_linear_constraints_vertex_indices.size(); ++i )
        {
            linear_constraints_indices.push_back( general_linear_constraints_vertex_indices.at(i) );
            for( int j = 0; j < linear_constraints_indices.back().size(); ++j )
            {
                linear_constraints_indices.back().at( j ) += N * general_linear_constraints_coordinate_offsets.at(i).at(j);
            }
            linear_constraints_coeffs.push_back( general_linear_constraints_coeffs.at(i) );
            linear_constraints_rhss.push_back( general_linear_constraints_rhss.at(i) );
        }
        assert( linear_constraints_indices.size() == num_linear_constraints );
        assert( linear_constraints_coeffs.size() == num_linear_constraints );
        assert( linear_constraints_rhss.size() == num_linear_constraints );
        
        solve_uncoupled_system(
            A, &b[ N*0 ], 0, &solution[ N*0 ],
            hard_constraints_indices,
            hard_constraints_values,
            linear_constraints_indices,
            linear_constraints_coeffs,
            linear_constraints_rhss,
            constraint_weight,
            system_handle
            );
    }
    
    // HACK 1: computing the sphere via graph laplacian, then switching to normal-motion-only 
    //         cotangent weights works well.
    //      2: computing the sphere via normal-motion-only cotangent weights, then switching to
    //         tangent-motion-only graph laplacian, then switching back to normal-motion-only
    //         cotangent weights works well, too.
    //      In fact, the two spheres are nearly identical.  They look better than the
    //      simply cotangent normal-motion-only solution and the simple graph laplacian solution.
    const bool kTangentMotionOnly = false;
    
    // Update the mesh.
    if( kTangentMotionOnly || params.vertex_normal_motion_only )
    {
        for( int i = 0; i < N; ++i )
        {
            TriMesh::Point& pos = mesh.point( TriMesh::VertexHandle( i ) );
            
            const TriMesh::Normal& n = mesh.normal( TriMesh::VertexHandle( i ) );
            assert( fabs( n.sqrnorm() - 1. ) < 1e-5 );
            
            const TriMesh::Point new_pos(
                solution.at( i + 0*N ),
                solution.at( i + 1*N ),
                solution.at( i + 2*N )
                );
            
            if( kTangentMotionOnly )
                pos = ( new_pos - pos ) - dot( ( new_pos - pos ), n ) * n + pos;
            else
                pos = dot( ( new_pos - pos ), n ) * n + pos;
        }
    }
    else
    {
        for( int i = 0; i < N; ++i )
        {
            TriMesh::Point& pos = mesh.point( TriMesh::VertexHandle( i ) );
            pos[0] = solution.at( i + 0*N );
            pos[1] = solution.at( i + 1*N );
            pos[2] = solution.at( i + 2*N );
        }
    }
    
#ifndef NDEBUG
    for( int i = 0; i < x_hard_constraints_indices.size(); ++i )
    {
        const int row = x_hard_constraints_indices.at(i);
        assert( fabs( solution.at( row + 0*N ) - x_hard_constraints_values.at(i) ) < 1e-5 );
    }
    for( int i = 0; i < y_hard_constraints_indices.size(); ++i )
    {
        const int row = y_hard_constraints_indices.at(i);
        assert( fabs( solution.at( row + 1*N ) - y_hard_constraints_values.at(i) ) < 1e-5 );
    }
    for( int i = 0; i < z_hard_constraints_indices.size(); ++i )
    {
        const int row = z_hard_constraints_indices.at(i);
        assert( fabs( solution.at( row + 2*N ) - z_hard_constraints_values.at(i) ) < 1e-5 );
    }
#endif
FN_EXIT}

// used
void get_current_curvatures(
    const TriMesh& mesh, const SparseMatrix_PETSc& L,
    ConstrainedSystem::Values& LMs_out
    )
{FN_ENTER
    const int N = mesh.n_vertices();
    
    assert( N > 0 );
    assert( N == L.GetNumColumns() );
    assert( L.GetNumRows() == L.GetNumColumns() );
    
    LMs_out.resize( N );
    
    ConstrainedSystem::Values laps( N*3 );
    {
        ConstrainedSystem::Values vs( N*3 );
        
        for( int i = 0; i < N; ++i )
        {
            const TriMesh::Point& point = mesh.point( TriMesh::VertexHandle( i ) );
            vs[ i + 0*N ] = point[0];
            vs[ i + 1*N ] = point[1];
            vs[ i + 2*N ] = point[2];
        }
        
        int ierr(0);
        // Create a dense PETSc matrix containing the normals.
        Mat vs_mat;
        ierr = MatCreateSeqDense( PETSC_COMM_WORLD, N, 3, &vs[0], &vs_mat ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        Mat laps_mat;
        ierr = MatCreateSeqDense( PETSC_COMM_WORLD, N, 3, &laps[0], &laps_mat ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        
        // PETSc aborts if this hasn't been called.
        ierr = MatAssemblyBegin( vs_mat, MAT_FINAL_ASSEMBLY ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        ierr = MatAssemblyEnd( vs_mat, MAT_FINAL_ASSEMBLY ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        ierr = MatAssemblyBegin( laps_mat, MAT_FINAL_ASSEMBLY ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        ierr = MatAssemblyEnd( laps_mat, MAT_FINAL_ASSEMBLY ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        
        ierr = MatMatMult( L.PETScMat(), vs_mat, MAT_REUSE_MATRIX, 1., &laps_mat ); CHKERRABORT( PETSC_COMM_WORLD, ierr );
        // Destroy the matrix.
        MatDestroy( vs_mat );
        MatDestroy( laps_mat );
    }
    
    for( int i = 0; i < N; ++i )
    {
        const TriMesh::Normal& n = mesh.normal( TriMesh::VertexHandle( i ) );
        LMs_out.at(i) = laps[ i + 0*N ]*n[0] + laps[ i + 1*N ]*n[1] + laps[ i + 2*N ]*n[2];
    }
FN_EXIT}

// used
void get_current_average_edge_lengths(
    const TriMesh& mesh,
    ConstrainedSystem::Values& avg_es_out
    )
{FN_ENTER
    const int N = mesh.n_vertices();
    
    assert( N > 0 );
    
    avg_es_out.resize( N );
    
    ConstrainedSystem::Indices neighs;
    for( int i = 0; i < N; ++i )
    {
        const TriMesh::Point& pi = mesh.point( TriMesh::VertexHandle(i) );
        
        get_one_ring( mesh, i, neighs );
        real_t avg_e = 0.;
        for( int j = 0; j < neighs.size(); ++j )
        {
            avg_e += ( pi - mesh.point( TriMesh::VertexHandle( neighs.at(j) ) ) ).length();
        }
        avg_e /= real_t( neighs.size() );
        
        avg_es_out.at(i) = avg_e;
    }
FN_EXIT}

// used
void calc_vertex_areas( const TriMesh& mesh, ConstrainedSystem::Values& areas_out )
{FN_ENTER
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    areas_out.resize( N );
    
    Vector3d vertex_pos;
    std::vector< int > one_ring_ids;
    std::vector< Vector3d > one_ring_pos;
    std::vector< Real > cot_alpha_out, cot_beta_out;
    for( int i = 0; i < N; ++i )
    {
        const bool has_boundary = get_one_ring_info(
            mesh, i,
            vertex_pos, one_ring_ids, one_ring_pos
            );
        
        cot_alpha_out.resize( one_ring_pos.size() );
        cot_beta_out.resize( one_ring_pos.size() );
        LaplacianEditing::cotangentWeights(
            vertex_pos, &one_ring_pos[0], one_ring_pos.size(),
            has_boundary,
            &cot_alpha_out[0], &cot_beta_out[0],
            &areas_out.at(i)
            );
    }
FN_EXIT}

// used
void gen_weighted_graph_laplacian_matrix( const TriMesh& mesh, SparseMatrix_PETSc** L, SparseMatrix_PETSc** L2 )
{FN_ENTER
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    assert( L && L2 );
    
    if( !*L )
    {
        std::vector< int > nnzs( N, 0 );
        for( int i = 0; i < N; ++i ) nnzs.at(i) = mesh.valence( TriMesh::VertexHandle(i) ) + 1;
        *L = new SparseMatrix_PETSc( N, N, nnzs );
    }
    else
    {
        (*L)->ZeroMatrix();
    }
    
    ConstrainedSystem::Values areas;
    calc_vertex_areas( mesh, areas );
    
    ConstrainedSystem::Indices indices;
    ConstrainedSystem::Values coeffs;
    for( int i = 0; i < N; ++i )
    {
        const real_t w = 1./areas.at(i);
        get_one_ring( mesh, i, indices );
        
        const real_t neigh_val = -w/indices.size();
        
        coeffs.clear();
        coeffs.resize( indices.size(), neigh_val );
        
        indices.push_back( i );
        coeffs.push_back( w );
        
        (*L)->SetRowValues( i, indices, coeffs );
    }
    (*L)->Finalize();
    
    
    // L*L seems to have ~2.7 times the number of non-zero entries.
    if( !*L2 )
    {
        Mat L2pet;
        int ierr = MatMatMult( (*L)->PETScMat(), (*L)->PETScMat(), MAT_INITIAL_MATRIX, 2.8 / 2., &L2pet );
        CHKERRABORT( PETSC_COMM_WORLD, ierr );
        *L2 = new SparseMatrix_PETSc( L2pet );
        MatDestroy( L2pet );
    }
    else
    {
        int ierr = MatMatMult( (*L)->PETScMat(), (*L)->PETScMat(), MAT_REUSE_MATRIX, 2.8 / 2., &(*L2)->PETScMat() );
        CHKERRABORT( PETSC_COMM_WORLD, ierr );
    }
    (*L2)->Finalize();
FN_EXIT}

// used
void gen_weighted_cotangent_laplacian_matrix( const TriMesh& mesh, SparseMatrix_PETSc** L, SparseMatrix_PETSc** L2 )
{FN_ENTER
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    assert( L && L2 );
    
    if( !*L )
    {
        std::vector< int > nnzs( N, 0 );
        for( int i = 0; i < N; ++i ) nnzs.at(i) = mesh.valence( TriMesh::VertexHandle(i) ) + 1;
        *L = new SparseMatrix_PETSc( N, N, nnzs );
    }
    else
    {
        (*L)->ZeroMatrix();
    }
    
    Vector3d vertex_pos;
    std::vector< int > one_ring_ids;
    std::vector< Vector3d > one_ring_pos;
    std::vector< Real > cot_alpha, cot_beta;
    
    ConstrainedSystem::Indices indices;
    ConstrainedSystem::Values coeffs;
    for( int i = 0; i < N; ++i )
    {
        const bool has_boundary = get_one_ring_info(
            mesh, i,
            vertex_pos, one_ring_ids, one_ring_pos
            );
        
        Real area;
        cot_alpha.resize( one_ring_pos.size() );
        cot_beta.resize( one_ring_pos.size() );
        LaplacianEditing::cotangentWeights(
            vertex_pos, &one_ring_pos[0], one_ring_pos.size(),
            has_boundary,
            &cot_alpha[0], &cot_beta[0],
            &area
            );
        
        indices.resize( one_ring_ids.size() + 1 );
        coeffs.resize( indices.size() );
        LaplacianEditing::Hn_coefficients(
            i,
            &one_ring_ids[0], one_ring_ids.size(),
            has_boundary,
            &cot_alpha[0], &cot_beta[0], area,
            &indices[0], &coeffs[0]
            );
        
        (*L)->SetRowValues( i, indices, coeffs );
    }
    (*L)->Finalize();
    
    
    // L*L seems to have ~2.7 times the number of non-zero entries.
    if( !*L2 )
    {
        Mat L2pet;
        int ierr = MatMatMult( (*L)->PETScMat(), (*L)->PETScMat(), MAT_INITIAL_MATRIX, 2.8 / 2., &L2pet );
        CHKERRABORT( PETSC_COMM_WORLD, ierr );
        *L2 = new SparseMatrix_PETSc( L2pet );
        MatDestroy( L2pet );
    }
    else
    {
        int ierr = MatMatMult( (*L)->PETScMat(), (*L)->PETScMat(), MAT_REUSE_MATRIX, 2.8 / 2., &(*L2)->PETScMat() );
        CHKERRABORT( PETSC_COMM_WORLD, ierr );
    }
    (*L2)->Finalize();
FN_EXIT}
void return_weighted_laplacian_matrix( SparseMatrix_PETSc* L, SparseMatrix_PETSc* L2 )
{FN_ENTER
    delete L;
    delete L2;
FN_EXIT}

bool petsc_initialized = false;

// This is a C++ port of fibermesh.py's optimize_surface2().
// It's (deliberately) not well-commented; please see fibermesh.py's optimize_surface2(), as this
// file is intended to simply reflect that.
//
// used
real_t optimize_surface2(
    TriMesh& mesh,
    
    
    const ConstrainedSystem::Indices& x_hard_constraints_indices,
    const ConstrainedSystem::Values& x_hard_constraints_values,
    
    const std::vector< ConstrainedSystem::Indices >& x_linear_constraints_indices,
    const std::vector< ConstrainedSystem::Values >& x_linear_constraints_coeffs,
    const ConstrainedSystem::Values& x_linear_constraints_rhss,
    
    
    const ConstrainedSystem::Indices& y_hard_constraints_indices,
    const ConstrainedSystem::Values& y_hard_constraints_values,
    
    const std::vector< ConstrainedSystem::Indices >& y_linear_constraints_indices,
    const std::vector< ConstrainedSystem::Values >& y_linear_constraints_coeffs,
    const ConstrainedSystem::Values& y_linear_constraints_rhss,
    
    
    const ConstrainedSystem::Indices& z_hard_constraints_indices,
    const ConstrainedSystem::Values& z_hard_constraints_values,
    
    const std::vector< ConstrainedSystem::Indices >& z_linear_constraints_indices,
    const std::vector< ConstrainedSystem::Values >& z_linear_constraints_coeffs,
    const ConstrainedSystem::Values& z_linear_constraints_rhss,
    
    
    const std::vector< ConstrainedSystem::Indices >& general_linear_constraints_vertex_indices_in,
    const std::vector< ConstrainedSystem::Indices >& general_linear_constraints_coordinate_offsets_in,
    const std::vector< ConstrainedSystem::Values >& general_linear_constraints_coeffs_in,
    const ConstrainedSystem::Values& general_linear_constraints_rhss_in,
    
    
    const optimize_surface2_params_t& params
    )
{FN_ENTER
    std::vector< ConstrainedSystem::Indices > general_linear_constraints_vertex_indices( general_linear_constraints_vertex_indices_in );
    std::vector< ConstrainedSystem::Indices > general_linear_constraints_coordinate_offsets( general_linear_constraints_coordinate_offsets_in );
    std::vector< ConstrainedSystem::Values > general_linear_constraints_coeffs( general_linear_constraints_coeffs_in );
    ConstrainedSystem::Values general_linear_constraints_rhss( general_linear_constraints_rhss_in );
    
    if( !petsc_initialized )
    {
        int argc = 0;
        char** argv = 0;
        PetscInitialize( &argc, &argv, 0, 0 );
        petsc_initialized = true;
    }
    
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    assert( params.max_iterations >= 0 );
    
    const real_t w_curvature = .1;
    //const real_t w_position = 10000.;
    const real_t w_position = 1e7;
    
    ConstrainedSystem::Indices range_len_N( N );
    for( int i = 0; i < N; ++i ) range_len_N.at(i) = i;
    
    // Initialize last_vs.
    std::vector< TriMesh::Point > last_vs( N );
    for( int i = 0; i < N; ++i ) last_vs.at(i) = mesh.point( TriMesh::VertexHandle(i) );
    // real_t min_max_diff = 1e20;
    real_t max_diff = params.absolute_convergence_threshold + 1.;
    
    // The list of general linear equations for constrained edges.
    // Specifically, a value in 'general_linear_constraint_avg_edge_eqs'
    // corresponds to the offset of the x linear constraint equation for the constrained edge
    // (vertex indices can be found by looking at
    // general_linear_constraints_vertex_indices.at(value)).
    // The y and z linear constraint equations directly follow the x equation.
    std::vector< int > general_linear_constraint_avg_edge_eqs;
    if( params.average_edge_lengths )
    for( int i = 0; i < x_hard_constraints_indices.size(); ++i )
    {
        const int vi = x_hard_constraints_indices.at(i);
        
        ConstrainedSystem::Indices two_int( 2 );
        ConstrainedSystem::Values two_real( 2 );
        
        ConstrainedSystem::Indices one_ring_ids;
        get_one_ring( mesh, vi, one_ring_ids );
        for( int j = 0; j < one_ring_ids.size(); ++j )
        {
            const int vj = one_ring_ids.at(j);
            
            // Skip vj that are also boundary constraints (and so the edge is
            // already totally constrained).
            if( std::find(
                    x_hard_constraints_indices.begin(), x_hard_constraints_indices.end(),
                    vj
                    ) != x_hard_constraints_indices.end()
                ) continue;
            
            general_linear_constraint_avg_edge_eqs.push_back( general_linear_constraints_vertex_indices.size() );
            
            two_int.at(0) = vi;
            two_int.at(1) = vj;
            general_linear_constraints_vertex_indices.push_back( two_int );
            general_linear_constraints_vertex_indices.push_back( two_int );
            general_linear_constraints_vertex_indices.push_back( two_int );
            
            two_int.at(0) = 0;
            two_int.at(1) = 0;
            general_linear_constraints_coordinate_offsets.push_back( two_int );
            two_int.at(0) = 1;
            two_int.at(1) = 1;
            general_linear_constraints_coordinate_offsets.push_back( two_int );
            two_int.at(0) = 2;
            two_int.at(1) = 2;
            general_linear_constraints_coordinate_offsets.push_back( two_int );
            
            two_real.at(0) = 1.;
            two_real.at(1) = -1.;
            general_linear_constraints_coeffs.push_back( two_real );
            general_linear_constraints_coeffs.push_back( two_real );
            general_linear_constraints_coeffs.push_back( two_real );
            
            general_linear_constraints_rhss.push_back( -31337. );
            general_linear_constraints_rhss.push_back( -31337. );
            general_linear_constraints_rhss.push_back( -31337. );
        }
    }
    
    SparseMatrix_PETSc* L(0);
    SparseMatrix_PETSc* L2(0);
    ConstrainedSystem* curvature_system(0);
    ConstrainedSystem* position_system(0);
    for(
        int iter = 0;
        iter < params.max_iterations && max_diff >= params.absolute_convergence_threshold;
        ++iter
        )
    {
        /*
        // This is wrong.  I don't want to change the (before weighting) matrix every iteration.
        // What I want to do is, begin with the graph laplacian matrix, then when we converge,
        // update the graph laplacian to a cotangent-weight laplacian matrix based on the converged
        // mesh, then when we converge with that, update the cotangent-weights again, etc, until
        // the whole thing converges.
        if( min_max_diff < 1e-4 )
            gen_weighted_cotangent_laplacian_matrix( mesh, &L, &L2 );
        else
        */
        if( params.graph_not_cotangent_weights )
            gen_weighted_graph_laplacian_matrix( mesh, &L, &L2 );
        else
            gen_weighted_cotangent_laplacian_matrix( mesh, &L, &L2 );
        
        ConstrainedSystem::Values LMs_cur;
        get_current_curvatures( mesh, *L, LMs_cur );
        
        ConstrainedSystem::Values LMs;
        compute_target_curvatures_mass_laplacian(
            mesh, *L, *L2, range_len_N, LMs_cur, w_curvature, LMs,
            &curvature_system
            );
        
        if( params.average_edge_lengths )
        {
            ConstrainedSystem::Values es_cur;
            get_current_average_edge_lengths( mesh, es_cur );
            
            ConstrainedSystem::Values es;
            compute_target_curvatures_mass_laplacian(
                mesh, *L, *L2, range_len_N, es_cur, w_curvature, es,
                // TODO: It's the same system, but a different constraint right-hand-side.
                &curvature_system
                );
            
            for( int i = 0; i < general_linear_constraint_avg_edge_eqs.size(); ++i )
            {
                const int glci = general_linear_constraint_avg_edge_eqs.at(i);
                
                assert( general_linear_constraints_vertex_indices.at(glci).size() == 2 );
                const int vi = general_linear_constraints_vertex_indices.at(glci).at(0);
                const int vj = general_linear_constraints_vertex_indices.at(glci).at(1);
                const TriMesh::Point nij =
                    (
                        mesh.point( TriMesh::VertexHandle(vi) ) -
                        mesh.point( TriMesh::VertexHandle(vj) )
                    ).normalize()
                    * .5 * ( es.at( vi ) + es.at( vj ) )
                    ;
                general_linear_constraints_rhss.at( glci + 0 ) = nij[0];
                general_linear_constraints_rhss.at( glci + 1 ) = nij[1];
                general_linear_constraints_rhss.at( glci + 2 ) = nij[2];
            }
        }
        
        optimize_positions_for_mass_laplacians_complex_constraints(
            mesh, *L, *L2, LMs,
            
            x_hard_constraints_indices,
            x_hard_constraints_values,
            
            x_linear_constraints_indices,
            x_linear_constraints_coeffs,
            x_linear_constraints_rhss,
            
            y_hard_constraints_indices,
            y_hard_constraints_values,
            
            y_linear_constraints_indices,
            y_linear_constraints_coeffs,
            y_linear_constraints_rhss,
            
            z_hard_constraints_indices,
            z_hard_constraints_values,
            
            z_linear_constraints_indices,
            z_linear_constraints_coeffs,
            z_linear_constraints_rhss,
            
            general_linear_constraints_vertex_indices,
            general_linear_constraints_coordinate_offsets,
            general_linear_constraints_coeffs,
            general_linear_constraints_rhss,
            
            params,
            
            w_position,
            
            // TODO: It's the same system, but a different constraint right-hand-side.
            ( params.average_edge_lengths ? 0 : &position_system )
            );
        mesh.update_normals( TriMesh::ANGLE_WEIGHTED );
        OpenMesh::IO::write_mesh( mesh, std::string("optimized") + stl::tostr( iter ) + std::string(".obj") );
        
        max_diff = 0.;
        for( int i = 0; i < mesh.n_vertices(); ++i )
        {
            const TriMesh::Point& pos_i = mesh.point( TriMesh::VertexHandle(i) );
            max_diff = std::max( max_diff, real_t( ( last_vs[i] - pos_i ).length() ) );
            last_vs[i] = pos_i;
        }
        std::cout << "max( | mesh.vs - last_vs | ): " << max_diff << '\n';
        // min_max_diff = std::min( min_max_diff, max_diff );
    }
    return_weighted_laplacian_matrix( L, L2 );
    L = 0;
    L2 = 0;
    delete curvature_system;
    delete position_system;
    
    return max_diff;
FN_EXIT}

// This is a C++ port of fibermesh.py's optimize_surface_linear().
// It's (deliberately) not well-commented; please see fibermesh.py's optimize_surface_linear(),
// as this file is intended to simply reflect that.
//
// used
real_t optimize_surface_linear(
    TriMesh& mesh,
    
    
    const ConstrainedSystem::Indices& x_hard_constraints_indices,
    const ConstrainedSystem::Values& x_hard_constraints_values,
    
    const std::vector< ConstrainedSystem::Indices >& x_linear_constraints_indices,
    const std::vector< ConstrainedSystem::Values >& x_linear_constraints_coeffs,
    const ConstrainedSystem::Values& x_linear_constraints_rhss,
    
    
    const ConstrainedSystem::Indices& y_hard_constraints_indices,
    const ConstrainedSystem::Values& y_hard_constraints_values,
    
    const std::vector< ConstrainedSystem::Indices >& y_linear_constraints_indices,
    const std::vector< ConstrainedSystem::Values >& y_linear_constraints_coeffs,
    const ConstrainedSystem::Values& y_linear_constraints_rhss,
    
    
    const ConstrainedSystem::Indices& z_hard_constraints_indices,
    const ConstrainedSystem::Values& z_hard_constraints_values,
    
    const std::vector< ConstrainedSystem::Indices >& z_linear_constraints_indices,
    const std::vector< ConstrainedSystem::Values >& z_linear_constraints_coeffs,
    const ConstrainedSystem::Values& z_linear_constraints_rhss,
    
    
    const std::vector< ConstrainedSystem::Indices >& general_linear_constraints_vertex_indices_in,
    const std::vector< ConstrainedSystem::Indices >& general_linear_constraints_coordinate_offsets_in,
    const std::vector< ConstrainedSystem::Values >& general_linear_constraints_coeffs_in,
    const ConstrainedSystem::Values& general_linear_constraints_rhss_in,
    
    
    const optimize_surface2_params_t& params
    )
{FN_ENTER
    std::vector< ConstrainedSystem::Indices > general_linear_constraints_vertex_indices( general_linear_constraints_vertex_indices_in );
    std::vector< ConstrainedSystem::Indices > general_linear_constraints_coordinate_offsets( general_linear_constraints_coordinate_offsets_in );
    std::vector< ConstrainedSystem::Values > general_linear_constraints_coeffs( general_linear_constraints_coeffs_in );
    ConstrainedSystem::Values general_linear_constraints_rhss( general_linear_constraints_rhss_in );
    
    if( !petsc_initialized )
    {
        int argc = 0;
        char** argv = 0;
        PetscInitialize( &argc, &argv, 0, 0 );
        petsc_initialized = true;
    }
    
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    //const real_t w_position = 10000.;
    const real_t w_position = 1e7;
    
    SparseMatrix_PETSc* L(0);
    SparseMatrix_PETSc* L2(0);
    
    if( params.graph_not_cotangent_weights )
        gen_weighted_graph_laplacian_matrix( mesh, &L, &L2 );
    else
        gen_weighted_cotangent_laplacian_matrix( mesh, &L, &L2 );
    
    assert( w_position >= 0. );
    
    assert( x_hard_constraints_indices.size() == x_hard_constraints_values.size() );
    assert( y_hard_constraints_indices.size() == y_hard_constraints_values.size() );
    assert( z_hard_constraints_indices.size() == z_hard_constraints_values.size() );
    
    assert( x_linear_constraints_indices.size() == x_linear_constraints_coeffs.size() );
    assert( x_linear_constraints_indices.size() == x_linear_constraints_rhss.size() );
    
    assert( y_linear_constraints_indices.size() == y_linear_constraints_coeffs.size() );
    assert( y_linear_constraints_indices.size() == y_linear_constraints_rhss.size() );
    
    assert( z_linear_constraints_indices.size() == z_linear_constraints_coeffs.size() );
    assert( z_linear_constraints_indices.size() == z_linear_constraints_rhss.size() );
    
    assert( general_linear_constraints_vertex_indices.size() == general_linear_constraints_coordinate_offsets.size() );
    assert( general_linear_constraints_vertex_indices.size() == general_linear_constraints_coeffs.size() );
    assert( general_linear_constraints_vertex_indices.size() == general_linear_constraints_rhss.size() );
    
#ifndef NDEBUG
    for( int i = 0; i < x_linear_constraints_indices.size(); ++i )
    {
        assert( x_linear_constraints_indices.at(i).size() == x_linear_constraints_coeffs.at(i).size() );
        assert_indices_inbounds( x_linear_constraints_indices.at(i).begin(), x_linear_constraints_indices.at(i).end(), 0, N-1 );
    }
    for( int i = 0; i < y_linear_constraints_indices.size(); ++i )
    {
        assert( y_linear_constraints_indices.at(i).size() == y_linear_constraints_coeffs.at(i).size() );
        assert_indices_inbounds( y_linear_constraints_indices.at(i).begin(), y_linear_constraints_indices.at(i).end(), 0, N-1 );
    }
    for( int i = 0; i < z_linear_constraints_indices.size(); ++i )
    {
        assert( z_linear_constraints_indices.at(i).size() == z_linear_constraints_coeffs.at(i).size() );
        assert_indices_inbounds( z_linear_constraints_indices.at(i).begin(), z_linear_constraints_indices.at(i).end(), 0, N-1 );
    }
    for( int i = 0; i < general_linear_constraints_vertex_indices.size(); ++i )
    {
        assert( general_linear_constraints_vertex_indices.at(i).size() == general_linear_constraints_coordinate_offsets.at(i).size() );
        assert( general_linear_constraints_vertex_indices.at(i).size() == general_linear_constraints_coeffs.at(i).size() );
        assert_indices_inbounds( general_linear_constraints_vertex_indices.at(i).begin(), general_linear_constraints_vertex_indices.at(i).end(), 0, N-1 );
        assert_indices_inbounds( general_linear_constraints_coordinate_offsets.at(i).begin(), general_linear_constraints_coordinate_offsets.at(i).end(), 0, 2 );
    }
#endif
    
    if( kPrintConstraints )
    {
        print_hard_constraints( x_hard_constraints_indices, x_hard_constraints_values );
        for( int i = 0; i < x_linear_constraints_indices.size(); ++i )
            print_linear_constraints(
                x_linear_constraints_indices.at(i),
                x_linear_constraints_coeffs.at(i),
                x_linear_constraints_rhss.at(i)
                );
        print_hard_constraints( y_hard_constraints_indices, y_hard_constraints_values );
        for( int i = 0; i < y_linear_constraints_indices.size(); ++i )
            print_linear_constraints(
                y_linear_constraints_indices.at(i),
                y_linear_constraints_coeffs.at(i),
                y_linear_constraints_rhss.at(i)
                );
        print_hard_constraints( z_hard_constraints_indices, z_hard_constraints_values );
        for( int i = 0; i < z_linear_constraints_indices.size(); ++i )
            print_linear_constraints(
                z_linear_constraints_indices.at(i),
                z_linear_constraints_coeffs.at(i),
                z_linear_constraints_rhss.at(i)
                );
        
        for( int i = 0; i < general_linear_constraints_vertex_indices.size(); ++i )
            print_linear_constraints(
                general_linear_constraints_vertex_indices.at(i),
                general_linear_constraints_coordinate_offsets.at(i),
                general_linear_constraints_coeffs.at(i),
                general_linear_constraints_rhss.at(i)
                );
    }
    
    {
        // Build a block diagonal matrix of size 3*N with 3 copies of L2
        std::vector< int > nnzs( 3*N, 0 );
        for( int i = 0; i < N; ++i )
        {
            const int nz = L2->GetRowNNZ( i );
            nnzs.at( i + 0*N ) = nz;
            nnzs.at( i + 1*N ) = nz;
            nnzs.at( i + 2*N ) = nz;
        }
        SparseMatrix_PETSc A( 3*N, 3*N, nnzs );
        {
            ConstrainedSystem::Indices indices;
            ConstrainedSystem::Values vals;
            ConstrainedSystem::Indices::iterator it;
            for( int i = 0; i < N; ++i )
            {
                L2->GetRowValues( i, indices, vals );
                // The first copy goes in the usual place.
                A.SetRowValues( i + 0*N, indices, vals );
                // The second copy has indices offset by N.
                for( it = indices.begin(); it != indices.end(); ++it ) *it += N;
                A.SetRowValues( i + 1*N, indices, vals );
                // The third copy has indices offset by N again.
                for( it = indices.begin(); it != indices.end(); ++it ) *it += N;
                A.SetRowValues( i + 2*N, indices, vals );
            }
        }
        A.Finalize();
        
        
        // Now build the list of all constraints.
        ConstrainedSystem::Indices hard_constraints_indices;
        ConstrainedSystem::Values hard_constraints_values;
        std::vector< ConstrainedSystem::Indices > linear_constraints_indices;
        std::vector< ConstrainedSystem::Values > linear_constraints_coeffs;
        ConstrainedSystem::Values linear_constraints_rhss;
        
        // Concatenate all hard constraints.
        const int num_hard_constraints =
            x_hard_constraints_indices.size() +
            y_hard_constraints_indices.size() +
            z_hard_constraints_indices.size()
            ;
        hard_constraints_indices.reserve( num_hard_constraints );
        hard_constraints_values.reserve( num_hard_constraints );
        for( int i = 0; i < x_hard_constraints_indices.size(); ++i )
        {
            hard_constraints_indices.push_back( 0*N + x_hard_constraints_indices.at(i) );
            hard_constraints_values.push_back( x_hard_constraints_values.at(i) );
        }
        for( int i = 0; i < y_hard_constraints_indices.size(); ++i )
        {
            hard_constraints_indices.push_back( 1*N + y_hard_constraints_indices.at(i) );
            hard_constraints_values.push_back( y_hard_constraints_values.at(i) );
        }
        for( int i = 0; i < z_hard_constraints_indices.size(); ++i )
        {
            hard_constraints_indices.push_back( 2*N + z_hard_constraints_indices.at(i) );
            hard_constraints_values.push_back( z_hard_constraints_values.at(i) );
        }
        assert( hard_constraints_indices.size() == hard_constraints_values.size() );
        assert( hard_constraints_indices.size() == num_hard_constraints );
        
        // Concatenate all linear constraints.
        const int num_linear_constraints =
            x_linear_constraints_indices.size() +
            y_linear_constraints_indices.size() +
            z_linear_constraints_indices.size() +
            general_linear_constraints_vertex_indices.size()
            ;
        linear_constraints_indices.reserve( num_linear_constraints );
        linear_constraints_coeffs.reserve( num_linear_constraints );
        linear_constraints_rhss.reserve( num_linear_constraints );
        ConstrainedSystem::Indices::iterator it;
        for( int i = 0; i < x_linear_constraints_indices.size(); ++i )
        {
            linear_constraints_indices.push_back( x_linear_constraints_indices.at(i) );
            linear_constraints_coeffs.push_back( x_linear_constraints_coeffs.at(i) );
            linear_constraints_rhss.push_back( x_linear_constraints_rhss.at(i) );
        }
        for( int i = 0; i < y_linear_constraints_indices.size(); ++i )
        {
            linear_constraints_indices.push_back( y_linear_constraints_indices.at(i) );
            for( it = linear_constraints_indices.back().begin(); it != linear_constraints_indices.back().end(); ++it ) *it += N;
            linear_constraints_coeffs.push_back( y_linear_constraints_coeffs.at(i) );
            linear_constraints_rhss.push_back( y_linear_constraints_rhss.at(i) );
        }
        for( int i = 0; i < z_linear_constraints_indices.size(); ++i )
        {
            linear_constraints_indices.push_back( z_linear_constraints_indices.at(i) );
            for( it = linear_constraints_indices.back().begin(); it != linear_constraints_indices.back().end(); ++it ) *it += 2*N;
            linear_constraints_coeffs.push_back( z_linear_constraints_coeffs.at(i) );
            linear_constraints_rhss.push_back( z_linear_constraints_rhss.at(i) );
        }
        for( int i = 0; i < general_linear_constraints_vertex_indices.size(); ++i )
        {
            linear_constraints_indices.push_back( general_linear_constraints_vertex_indices.at(i) );
            for( int j = 0; j < linear_constraints_indices.back().size(); ++j )
            {
                linear_constraints_indices.back().at( j ) += N * general_linear_constraints_coordinate_offsets.at(i).at(j);
            }
            linear_constraints_coeffs.push_back( general_linear_constraints_coeffs.at(i) );
            linear_constraints_rhss.push_back( general_linear_constraints_rhss.at(i) );
        }
        assert( linear_constraints_indices.size() == num_linear_constraints );
        assert( linear_constraints_coeffs.size() == num_linear_constraints );
        assert( linear_constraints_rhss.size() == num_linear_constraints );
        
        ConstrainedSystem::Values solution( N*3 );
        ConstrainedSystem::Values b( N*3 );
        solve_uncoupled_system(
            A, &b[ N*0 ], 0, &solution[ N*0 ],
            hard_constraints_indices,
            hard_constraints_values,
            linear_constraints_indices,
            linear_constraints_coeffs,
            linear_constraints_rhss,
            w_position
            );
        
        // Update the mesh.
        if( params.vertex_normal_motion_only )
        {
            for( int i = 0; i < N; ++i )
            {
                TriMesh::Point& pos = mesh.point( TriMesh::VertexHandle( i ) );
                
                const TriMesh::Normal& n = mesh.normal( TriMesh::VertexHandle( i ) );
                assert( fabs( n.sqrnorm() - 1. ) < 1e-5 );
                
                const TriMesh::Point new_pos(
                    solution.at( i + 0*N ),
                    solution.at( i + 1*N ),
                    solution.at( i + 2*N )
                    );
                
                pos = dot( ( new_pos - pos ), n ) * n + pos;
            }
        }
        else
        {
            for( int i = 0; i < N; ++i )
            {
                TriMesh::Point& pos = mesh.point( TriMesh::VertexHandle( i ) );
                pos[0] = solution.at( i + 0*N );
                pos[1] = solution.at( i + 1*N );
                pos[2] = solution.at( i + 2*N );
            }
        }
        
#ifndef NDEBUG
        for( int i = 0; i < x_hard_constraints_indices.size(); ++i )
        {
            const int row = x_hard_constraints_indices.at(i);
            assert( fabs( solution.at( row + 0*N ) - x_hard_constraints_values.at(i) ) < 1e-5 );
        }
        for( int i = 0; i < y_hard_constraints_indices.size(); ++i )
        {
            const int row = y_hard_constraints_indices.at(i);
            assert( fabs( solution.at( row + 1*N ) - y_hard_constraints_values.at(i) ) < 1e-5 );
        }
        for( int i = 0; i < z_hard_constraints_indices.size(); ++i )
        {
            const int row = z_hard_constraints_indices.at(i);
            assert( fabs( solution.at( row + 2*N ) - z_hard_constraints_values.at(i) ) < 1e-5 );
        }
#endif
    }
    
    mesh.update_normals( TriMesh::ANGLE_WEIGHTED );
    OpenMesh::IO::write_mesh( mesh, std::string("optimized-linear.obj") );
    
    return_weighted_laplacian_matrix( L, L2 );
    L = 0;
    L2 = 0;
    
    return 0;
FN_EXIT}
} // ~namespace fibermesh

namespace
{
// used
void copy_linear_parameters(
    // source
    const int_t linear_constraints_num,
    const int_t* linear_constraints_eq_num,
    const int_t* linear_constraints_indices,
    // may be NULL
    const int_t* linear_constraints_coordinate_offsets,
    const real_t* linear_constraints_coeffs,
    const real_t* linear_constraints_rhss,
    
    // destination
    std::vector< ConstrainedSystem::Indices >& linear_constraints_indices_in,
    // will be untouched if linear_constraints_coordinate_offsets is NULL, so pass in anything, even 'linear_constraints_indices_in' again.
    std::vector< ConstrainedSystem::Indices >& linear_constraints_coordinate_offsets_in,
    std::vector< ConstrainedSystem::Values >& linear_constraints_coeffs_in,
    ConstrainedSystem::Values& linear_constraints_rhss_in
    )
{FN_ENTER
    linear_constraints_indices_in.resize( linear_constraints_num );
    linear_constraints_coeffs_in.resize( linear_constraints_num );
    if( linear_constraints_coordinate_offsets ) linear_constraints_coordinate_offsets_in.resize( linear_constraints_num );
    linear_constraints_rhss_in.resize( linear_constraints_num );
    int off = 0;
    for( int i = 0; i < linear_constraints_num; ++i )
    {
        const int eq_num = linear_constraints_eq_num[i];
        linear_constraints_indices_in.at(i).resize( eq_num );
        linear_constraints_coeffs_in.at(i).resize( eq_num );
        std::copy( linear_constraints_indices + off, linear_constraints_indices + off + eq_num, linear_constraints_indices_in.at(i).begin() );
        std::copy( linear_constraints_coeffs + off, linear_constraints_coeffs + off + eq_num, linear_constraints_coeffs_in.at(i).begin() );
        
        if( linear_constraints_coordinate_offsets )
        {
            linear_constraints_coordinate_offsets_in.at(i).resize( eq_num );
            std::copy( linear_constraints_coordinate_offsets + off, linear_constraints_coordinate_offsets + off + eq_num, linear_constraints_coordinate_offsets_in.at(i).begin() );
        }
        
        off += eq_num;
    }
    std::copy( linear_constraints_rhss, linear_constraints_rhss + linear_constraints_num, linear_constraints_rhss_in.begin() );
FN_EXIT}

// used
void copy_parameters(
    // source
    const int_t hard_constraints_num,
    const int_t* hard_constraints_indices,
    const real_t* hard_constraints_values,
    
    const int_t linear_constraints_num,
    const int_t* linear_constraints_eq_num,
    const int_t* linear_constraints_indices,
    const real_t* linear_constraints_coeffs,
    const real_t* linear_constraints_rhss,
    
    // destination
    ConstrainedSystem::Indices& hard_constraints_indices_in,
    ConstrainedSystem::Values& hard_constraints_values_in,
    std::vector< ConstrainedSystem::Indices >& linear_constraints_indices_in,
    std::vector< ConstrainedSystem::Values >& linear_constraints_coeffs_in,
    ConstrainedSystem::Values& linear_constraints_rhss_in
    )
{FN_ENTER
    hard_constraints_indices_in.resize( hard_constraints_num );
    std::copy( hard_constraints_indices, hard_constraints_indices + hard_constraints_num, hard_constraints_indices_in.begin() );
    hard_constraints_values_in.resize( hard_constraints_num );
    std::copy( hard_constraints_values, hard_constraints_values + hard_constraints_num, hard_constraints_values_in.begin() );
    
    copy_linear_parameters(
        linear_constraints_num,
        linear_constraints_eq_num,
        linear_constraints_indices,
        0,
        linear_constraints_coeffs,
        linear_constraints_rhss,
        
        linear_constraints_indices_in,
        linear_constraints_indices_in,
        linear_constraints_coeffs_in,
        linear_constraints_rhss_in
        );
FN_EXIT}
}

extern "C"
{
// untested
void initialize_optimize_surface2_params( struct optimize_surface2_params_t* params )
{
    assert( params );
    params->max_iterations = 4;
    params->absolute_convergence_threshold = 0.;
    params->graph_not_cotangent_weights = true;
    params->average_edge_lengths = false;
    params->vertex_normal_motion_only = false;
}

// used, and constraint transmission has been tested.
real_t optimize_surface2(
    void* mesh_handle,
    
    
    // These are a pair of equal-length 1d arrays.  The num parameter is the length.
    int_t x_hard_constraints_num,
    int_t* x_hard_constraints_indices,
    real_t* x_hard_constraints_values,
    
    // These are three arrays.  The indices and coeffs parameters represent an array of sub-arrays,
    // but have been flattened end to end into two 1d arrays.
    // There are x_linear_constraints_num arrays and the length of each sub-array is stored in
    // x_linear_constraints_eq_num.
    // The rhss parameter is a 1d array with x_linear_constraints_num values.
    int_t x_linear_constraints_num,
    int_t* x_linear_constraints_eq_num,
    int_t* x_linear_constraints_indices,
    real_t* x_linear_constraints_coeffs,
    real_t* x_linear_constraints_rhss,
    
    
    int_t y_hard_constraints_num,
    int_t* y_hard_constraints_indices,
    real_t* y_hard_constraints_values,
    
    int_t y_linear_constraints_num,
    int_t* y_linear_constraints_eq_num,
    int_t* y_linear_constraints_indices,
    real_t* y_linear_constraints_coeffs,
    real_t* y_linear_constraints_rhss,
    
    
    int_t z_hard_constraints_num,
    int_t* z_hard_constraints_indices,
    real_t* z_hard_constraints_values,
    
    int_t z_linear_constraints_num,
    int_t* z_linear_constraints_eq_num,
    int_t* z_linear_constraints_indices,
    real_t* z_linear_constraints_coeffs,
    real_t* z_linear_constraints_rhss,
    
    
    // These are four arrays.  The vertex_indices and coordinate_offsets and coeffs parameters
    // represent an array of sub-arrays, but have been flattened end to end into three 1d arrays.
    // There are general_linear_constraints_num arrays and the length of each sub-array is stored in
    // general_linear_constraints_eq_num.
    // The rhss parameter is a 1d array with general_linear_constraints_num values.
    int_t general_linear_constraints_num,
    int_t* general_linear_constraints_eq_num,
    int_t* general_linear_constraints_vertex_indices,
    int_t* general_linear_constraints_coordinate_offsets,
    real_t* general_linear_constraints_coeffs,
    real_t* general_linear_constraints_rhss,
    
    
    struct optimize_surface2_params_t* params
    )
{FN_ENTER
    assert( mesh_handle );
    assert( params );
    
    TriMesh& mesh = *static_cast< TriMesh* >( mesh_handle );
    
    ConstrainedSystem::Indices x_hard_constraints_indices_in;
    ConstrainedSystem::Indices y_hard_constraints_indices_in;
    ConstrainedSystem::Indices z_hard_constraints_indices_in;
    
    ConstrainedSystem::Values x_hard_constraints_values_in;
    ConstrainedSystem::Values y_hard_constraints_values_in;
    ConstrainedSystem::Values z_hard_constraints_values_in;
    
    std::vector< ConstrainedSystem::Indices > x_linear_constraints_indices_in;
    std::vector< ConstrainedSystem::Indices > y_linear_constraints_indices_in;
    std::vector< ConstrainedSystem::Indices > z_linear_constraints_indices_in;
    
    std::vector< ConstrainedSystem::Values > x_linear_constraints_coeffs_in;
    std::vector< ConstrainedSystem::Values > y_linear_constraints_coeffs_in;
    std::vector< ConstrainedSystem::Values > z_linear_constraints_coeffs_in;
    
    ConstrainedSystem::Values x_linear_constraints_rhss_in;
    ConstrainedSystem::Values y_linear_constraints_rhss_in;
    ConstrainedSystem::Values z_linear_constraints_rhss_in;
    
    copy_parameters(
        x_hard_constraints_num,
        x_hard_constraints_indices,
        x_hard_constraints_values,
        
        x_linear_constraints_num,
        x_linear_constraints_eq_num,
        x_linear_constraints_indices,
        x_linear_constraints_coeffs,
        x_linear_constraints_rhss,
        
        x_hard_constraints_indices_in,
        x_hard_constraints_values_in,
        x_linear_constraints_indices_in,
        x_linear_constraints_coeffs_in,
        x_linear_constraints_rhss_in
        );
    
    copy_parameters(
        y_hard_constraints_num,
        y_hard_constraints_indices,
        y_hard_constraints_values,
        
        y_linear_constraints_num,
        y_linear_constraints_eq_num,
        y_linear_constraints_indices,
        y_linear_constraints_coeffs,
        y_linear_constraints_rhss,
        
        y_hard_constraints_indices_in,
        y_hard_constraints_values_in,
        y_linear_constraints_indices_in,
        y_linear_constraints_coeffs_in,
        y_linear_constraints_rhss_in
        );
    
    copy_parameters(
        z_hard_constraints_num,
        z_hard_constraints_indices,
        z_hard_constraints_values,
        
        z_linear_constraints_num,
        z_linear_constraints_eq_num,
        z_linear_constraints_indices,
        z_linear_constraints_coeffs,
        z_linear_constraints_rhss,
        
        z_hard_constraints_indices_in,
        z_hard_constraints_values_in,
        z_linear_constraints_indices_in,
        z_linear_constraints_coeffs_in,
        z_linear_constraints_rhss_in
        );
    
    
    std::vector< ConstrainedSystem::Indices > general_linear_constraints_vertex_indices_in;
    std::vector< ConstrainedSystem::Indices > general_linear_constraints_coordinate_offsets_in;
    std::vector< ConstrainedSystem::Values > general_linear_constraints_coeffs_in;
    ConstrainedSystem::Values general_linear_constraints_rhss_in;
    
    copy_linear_parameters(
        general_linear_constraints_num,
        general_linear_constraints_eq_num,
        general_linear_constraints_vertex_indices,
        general_linear_constraints_coordinate_offsets,
        general_linear_constraints_coeffs,
        general_linear_constraints_rhss,
        
        general_linear_constraints_vertex_indices_in,
        general_linear_constraints_coordinate_offsets_in,
        general_linear_constraints_coeffs_in,
        general_linear_constraints_rhss_in
        );
    
    return fibermesh::optimize_surface2(
        mesh,
        
        
        x_hard_constraints_indices_in,
        x_hard_constraints_values_in,
        
        x_linear_constraints_indices_in,
        x_linear_constraints_coeffs_in,
        x_linear_constraints_rhss_in,
        
        
        y_hard_constraints_indices_in,
        y_hard_constraints_values_in,
        
        y_linear_constraints_indices_in,
        y_linear_constraints_coeffs_in,
        y_linear_constraints_rhss_in,
        
        
        z_hard_constraints_indices_in,
        z_hard_constraints_values_in,
        
        z_linear_constraints_indices_in,
        z_linear_constraints_coeffs_in,
        z_linear_constraints_rhss_in,
        
        
        general_linear_constraints_vertex_indices_in,
        general_linear_constraints_coordinate_offsets_in,
        general_linear_constraints_coeffs_in,
        general_linear_constraints_rhss_in,
        
        
        *params
        );
FN_EXIT}

// used, and constraint transmission has been tested.
real_t optimize_surface_linear(
    void* mesh_handle,
    
    
    // These are a pair of equal-length 1d arrays.  The num parameter is the length.
    int_t x_hard_constraints_num,
    int_t* x_hard_constraints_indices,
    real_t* x_hard_constraints_values,
    
    // These are three arrays.  The indices and coeffs parameters represent an array of sub-arrays,
    // but have been flattened end to end into two 1d arrays.
    // There are x_linear_constraints_num arrays and the length of each sub-array is stored in
    // x_linear_constraints_eq_num.
    // The rhss parameter is a 1d array with x_linear_constraints_num values.
    int_t x_linear_constraints_num,
    int_t* x_linear_constraints_eq_num,
    int_t* x_linear_constraints_indices,
    real_t* x_linear_constraints_coeffs,
    real_t* x_linear_constraints_rhss,
    
    
    int_t y_hard_constraints_num,
    int_t* y_hard_constraints_indices,
    real_t* y_hard_constraints_values,
    
    int_t y_linear_constraints_num,
    int_t* y_linear_constraints_eq_num,
    int_t* y_linear_constraints_indices,
    real_t* y_linear_constraints_coeffs,
    real_t* y_linear_constraints_rhss,
    
    
    int_t z_hard_constraints_num,
    int_t* z_hard_constraints_indices,
    real_t* z_hard_constraints_values,
    
    int_t z_linear_constraints_num,
    int_t* z_linear_constraints_eq_num,
    int_t* z_linear_constraints_indices,
    real_t* z_linear_constraints_coeffs,
    real_t* z_linear_constraints_rhss,
    
    
    // These are four arrays.  The vertex_indices and coordinate_offsets and coeffs parameters
    // represent an array of sub-arrays, but have been flattened end to end into three 1d arrays.
    // There are general_linear_constraints_num arrays and the length of each sub-array is stored in
    // general_linear_constraints_eq_num.
    // The rhss parameter is a 1d array with general_linear_constraints_num values.
    int_t general_linear_constraints_num,
    int_t* general_linear_constraints_eq_num,
    int_t* general_linear_constraints_vertex_indices,
    int_t* general_linear_constraints_coordinate_offsets,
    real_t* general_linear_constraints_coeffs,
    real_t* general_linear_constraints_rhss,
    
    
    struct optimize_surface2_params_t* params
    )
{FN_ENTER
    assert( mesh_handle );
    assert( params );
    
    TriMesh& mesh = *static_cast< TriMesh* >( mesh_handle );
    
    ConstrainedSystem::Indices x_hard_constraints_indices_in;
    ConstrainedSystem::Indices y_hard_constraints_indices_in;
    ConstrainedSystem::Indices z_hard_constraints_indices_in;
    
    ConstrainedSystem::Values x_hard_constraints_values_in;
    ConstrainedSystem::Values y_hard_constraints_values_in;
    ConstrainedSystem::Values z_hard_constraints_values_in;
    
    std::vector< ConstrainedSystem::Indices > x_linear_constraints_indices_in;
    std::vector< ConstrainedSystem::Indices > y_linear_constraints_indices_in;
    std::vector< ConstrainedSystem::Indices > z_linear_constraints_indices_in;
    
    std::vector< ConstrainedSystem::Values > x_linear_constraints_coeffs_in;
    std::vector< ConstrainedSystem::Values > y_linear_constraints_coeffs_in;
    std::vector< ConstrainedSystem::Values > z_linear_constraints_coeffs_in;
    
    ConstrainedSystem::Values x_linear_constraints_rhss_in;
    ConstrainedSystem::Values y_linear_constraints_rhss_in;
    ConstrainedSystem::Values z_linear_constraints_rhss_in;
    
    copy_parameters(
        x_hard_constraints_num,
        x_hard_constraints_indices,
        x_hard_constraints_values,
        
        x_linear_constraints_num,
        x_linear_constraints_eq_num,
        x_linear_constraints_indices,
        x_linear_constraints_coeffs,
        x_linear_constraints_rhss,
        
        x_hard_constraints_indices_in,
        x_hard_constraints_values_in,
        x_linear_constraints_indices_in,
        x_linear_constraints_coeffs_in,
        x_linear_constraints_rhss_in
        );
    
    copy_parameters(
        y_hard_constraints_num,
        y_hard_constraints_indices,
        y_hard_constraints_values,
        
        y_linear_constraints_num,
        y_linear_constraints_eq_num,
        y_linear_constraints_indices,
        y_linear_constraints_coeffs,
        y_linear_constraints_rhss,
        
        y_hard_constraints_indices_in,
        y_hard_constraints_values_in,
        y_linear_constraints_indices_in,
        y_linear_constraints_coeffs_in,
        y_linear_constraints_rhss_in
        );
    
    copy_parameters(
        z_hard_constraints_num,
        z_hard_constraints_indices,
        z_hard_constraints_values,
        
        z_linear_constraints_num,
        z_linear_constraints_eq_num,
        z_linear_constraints_indices,
        z_linear_constraints_coeffs,
        z_linear_constraints_rhss,
        
        z_hard_constraints_indices_in,
        z_hard_constraints_values_in,
        z_linear_constraints_indices_in,
        z_linear_constraints_coeffs_in,
        z_linear_constraints_rhss_in
        );
    
    
    std::vector< ConstrainedSystem::Indices > general_linear_constraints_vertex_indices_in;
    std::vector< ConstrainedSystem::Indices > general_linear_constraints_coordinate_offsets_in;
    std::vector< ConstrainedSystem::Values > general_linear_constraints_coeffs_in;
    ConstrainedSystem::Values general_linear_constraints_rhss_in;
    
    copy_linear_parameters(
        general_linear_constraints_num,
        general_linear_constraints_eq_num,
        general_linear_constraints_vertex_indices,
        general_linear_constraints_coordinate_offsets,
        general_linear_constraints_coeffs,
        general_linear_constraints_rhss,
        
        general_linear_constraints_vertex_indices_in,
        general_linear_constraints_coordinate_offsets_in,
        general_linear_constraints_coeffs_in,
        general_linear_constraints_rhss_in
        );
    
    return fibermesh::optimize_surface_linear(
        mesh,
        
        
        x_hard_constraints_indices_in,
        x_hard_constraints_values_in,
        
        x_linear_constraints_indices_in,
        x_linear_constraints_coeffs_in,
        x_linear_constraints_rhss_in,
        
        
        y_hard_constraints_indices_in,
        y_hard_constraints_values_in,
        
        y_linear_constraints_indices_in,
        y_linear_constraints_coeffs_in,
        y_linear_constraints_rhss_in,
        
        
        z_hard_constraints_indices_in,
        z_hard_constraints_values_in,
        
        z_linear_constraints_indices_in,
        z_linear_constraints_coeffs_in,
        z_linear_constraints_rhss_in,
        
        
        general_linear_constraints_vertex_indices_in,
        general_linear_constraints_coordinate_offsets_in,
        general_linear_constraints_coeffs_in,
        general_linear_constraints_rhss_in,
        
        
        *params
        );
FN_EXIT}

// Given a flattened list of vertices (succession of 'n_vertices' xyz coordinates)
// and a flattened list of faces (succession of 'n_faces' vertex indices),
// returns a new OpenMesh TriMesh handle.
//
// tested
void* new_openmesh_from_data( int_t n_vertices, real_t* vertices, int_t n_faces, int_t* faces )
{FN_ENTER
    assert( n_vertices > 0 );
    assert( n_faces > 0 );
    
    assert( vertices );
    assert( faces );
    
    TriMesh* mesh = new TriMesh();
    
    std::vector< TriMesh::VertexHandle > vertex_handles( n_vertices );
    for( int i = 0; i < n_vertices; ++i )
    {
        vertex_handles[i] = mesh->add_vertex( TriMesh::Point(
            vertices[ 3*i + 0 ],
            vertices[ 3*i + 1 ],
            vertices[ 3*i + 2 ]
            ) );
    }
    
    std::vector< TriMesh::VertexHandle > face( 3 );
    for( int i = 0; i < n_faces; ++i )
    {
        face[0] = vertex_handles.at( faces[ 3*i + 0 ] );
        face[1] = vertex_handles.at( faces[ 3*i + 1 ] );
        face[2] = vertex_handles.at( faces[ 3*i + 2 ] );
        
        mesh->add_face( face );
    }
    
    // According to the tutorial "OpenMesh/Doc/Tutorial/01-build_cube/build_cube.cc", we're done.
    // (That is, the tutorial doesn't call any sort of update-data-structures function.)
    
    mesh->update_normals( TriMesh::ANGLE_WEIGHTED );
    return mesh;
FN_EXIT}

// Frees the memory associated with the OpenMesh TriMesh.
//
// tested
void delete_openmesh( void* mesh_handle )
{FN_ENTER
    TriMesh* mesh = static_cast< TriMesh* >( mesh_handle );
    delete mesh;
FN_EXIT}

// Updates the OpenMesh TriMesh 'mesh_handle' with the data in 'vertices_in', a succession of
// mesh->n_vertices() xyz coordinates.
//
// tested
void update_openmesh_vertices( void* mesh_handle, real_t* vertices_in )
{FN_ENTER
    assert( mesh_handle );
    TriMesh& mesh = *static_cast< TriMesh* >( mesh_handle );
    
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    assert( vertices_in );
    
    for( int i = 0; i < N; ++i )
    {
        TriMesh::Point& pos_i = mesh.point( TriMesh::VertexHandle(i) );
        pos_i[0] = vertices_in[ 3*i + 0 ];
        pos_i[1] = vertices_in[ 3*i + 1 ];
        pos_i[2] = vertices_in[ 3*i + 2 ];
    }
    
    mesh.update_normals( TriMesh::ANGLE_WEIGHTED );
FN_EXIT}

// Copies the vertices of the OpenMesh TriMesh 'mesh_handle' into 'vertices_out' as a succession of
// mesh->n_vertices() xyz coordinates.
//
// tested
void extract_openmesh_vertices( void* mesh_handle, real_t* vertices_out )
{FN_ENTER
    assert( mesh_handle );
    TriMesh& mesh = *static_cast< TriMesh* >( mesh_handle );
    
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    assert( vertices_out );
    
    for( int i = 0; i < N; ++i )
    {
        const TriMesh::Point& pos_i = mesh.point( TriMesh::VertexHandle(i) );
        vertices_out[ 3*i + 0 ] = pos_i[0];
        vertices_out[ 3*i + 1 ] = pos_i[1];
        vertices_out[ 3*i + 2 ] = pos_i[2];
    }
FN_EXIT}

// Copies the vertex normals of the OpenMesh TriMesh 'mesh_handle' into 'vertex_normals_out'
// as a succession of mesh->n_vertices() xyz coordinates.
//
// tested
void extract_openmesh_vertex_normals( void* mesh_handle, real_t* vertex_normals_out )
{FN_ENTER
    assert( mesh_handle );
    TriMesh& mesh = *static_cast< TriMesh* >( mesh_handle );
    
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    assert( vertex_normals_out );
    
    // NOTE: We don't need to call update_vertex_normals() (optionally with TriMesh::ANGLE_WEIGHTED),
    //       because update_openmesh_vertices() does it for us.
    
    for( int i = 0; i < N; ++i )
    {
        const TriMesh::Normal& normal_i = mesh.normal( TriMesh::VertexHandle(i) );
        vertex_normals_out[ 3*i + 0 ] = normal_i[0];
        vertex_normals_out[ 3*i + 1 ] = normal_i[1];
        vertex_normals_out[ 3*i + 2 ] = normal_i[2];
    }
FN_EXIT}

// Copies the face normals of the OpenMesh TriMesh 'mesh_handle' into 'face_normals_out'
// as a succession of mesh->n_faces() xyz coordinates.
//
// tested
void extract_openmesh_face_normals( void* mesh_handle, real_t* face_normals_out )
{FN_ENTER
    assert( mesh_handle );
    TriMesh& mesh = *static_cast< TriMesh* >( mesh_handle );
    
    const int N = mesh.n_faces();
    assert( N > 0 );
    
    assert( face_normals_out );
    
    // NOTE: We don't need to call update_face_normals() because update_openmesh_vertices()
    //       does it for us.
    
    for( int i = 0; i < N; ++i )
    {
        const TriMesh::Normal& normal_i = mesh.normal( TriMesh::FaceHandle(i) );
        face_normals_out[ 3*i + 0 ] = normal_i[0];
        face_normals_out[ 3*i + 1 ] = normal_i[1];
        face_normals_out[ 3*i + 2 ] = normal_i[2];
    }
FN_EXIT}

// Copies the face areas of the OpenMesh TriMesh 'mesh_handle' into 'face_areas_out'
// as a succession of mesh->n_faces() scalars.
//
// tested
void extract_openmesh_face_areas( void* mesh_handle, real_t* face_areas_out )
{FN_ENTER
    assert( mesh_handle );
    TriMesh& mesh = *static_cast< TriMesh* >( mesh_handle );
    
    const int N = mesh.n_faces();
    assert( N > 0 );
    
    assert( face_areas_out );
    
    TriMesh::Point face_points[3];
    for( int i = 0; i < N; ++i )
    {
        int v = 0;
        for( TriMesh::ConstFaceVertexIter fvit = mesh.cfv_iter( TriMesh::FaceHandle( i ) ); fvit; ++fvit, ++v )
        {
            face_points[v] = mesh.point( fvit.handle() );
        }
        
        face_areas_out[ i ] = .5 * cross(
            face_points[1] - face_points[0],
            face_points[2] - face_points[0]
            ).length();
    }
FN_EXIT}

// Returns the number of vertices in the OpenMesh TriMesh 'mesh_handle'.
//
// tested
int_t openmesh_num_vertices( void* mesh_handle )
{FN_ENTER
    assert( mesh_handle );
    TriMesh& mesh = *static_cast< TriMesh* >( mesh_handle );
    
    const int N = mesh.n_vertices();
    assert( N > 0 );
    
    return N;
FN_EXIT}

// Returns the number of faces in the OpenMesh TriMesh 'mesh_handle'.
//
// tested
int_t openmesh_num_faces( void* mesh_handle )
{FN_ENTER
    assert( mesh_handle );
    TriMesh& mesh = *static_cast< TriMesh* >( mesh_handle );
    
    const int N = mesh.n_faces();
    assert( N > 0 );
    
    return N;
FN_EXIT}
}
