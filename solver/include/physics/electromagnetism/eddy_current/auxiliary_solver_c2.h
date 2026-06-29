#pragma once

#include "math/fem/block_rack.h"
#include "math/fem/bc_dirichlet.h"


namespace simu::T_Omega_AMS_c {

typedef struct AMS_Context 
{
    Mat P, G, I, LpM;
    Mat H1, H2, J;
    Mat Qc, Qi;
    //Mat X;

    KSP inner_LpM_ksp;
    KSP inner_Q1_ksp, inner_Q2_ksp;   
    Vec tmp_1, tmp_2, tmp_3, tmp_4;
    Vec rho,   gamma;                 // size 3*nV (L-space RHS / sol)
    Vec rho_1,   gamma_1;                 // size 3*nV (L-space RHS / sol)
    Vec rho_2,   gamma_2;                 // size 3*nV (L-space RHS / sol)
    Vec rho_3,   gamma_3;                 // size 3*nV (L-space RHS / sol)
    Vec rho_4,   gamma_4;                 // size 3*nV (L-space RHS / sol)

    Vec zeta,  kappa;                 // size   nV (Q-space RHS / sol)
    Block_Rack* br_system;

    Vec zeta_1,  kappa_1;                 // size   nV (Q-space RHS / sol)
    Vec zeta_2,  kappa_2;                 // size   nV (Q-space RHS / sol)
    Vec zeta_3,  kappa_3;                 // size   nV (Q-space RHS / sol)
    Vec zeta_4,  kappa_4;                 // size   nV (Q-space RHS / sol)



    std::vector<Vec> block_x{4};
    std::vector<Vec> block_r{4};



    KSP inner_Qc_ksp, inner_Qi_ksp;     // built once, reused per PCApply
    //KSP inner_X_ksp; 

    Dirichlet_BC* bc_T_1_H1_s;
    Dirichlet_BC* bc_T_1_H1_v;
    Dirichlet_BC* bc_O_o;
    Dirichlet_BC* bc_O_i;
    Dirichlet_BC* global_H1;


} AMS_Context;

typedef struct AMS_Info 
{
    PetscInt  n_iteration = 0;
    PetscReal condition_number = 0.;
} AMS_Info;




PetscErrorCode AMS_apply_decoupled(PC pc, Vec r, Vec x);
PetscErrorCode AMS_apply_coupled(PC pc, Vec r, Vec x);
PetscErrorCode AMS_apply_global(PC pc, Vec r, Vec x);
PetscErrorCode AMS_apply_fully_coupled(PC pc, Vec r, Vec x);



PetscErrorCode AMS_destroy(PC pc);




PetscErrorCode solve_AMS(
    AMS_Info& AMS_info,
    int algorithm_id,
    Mat P, Mat G, Mat I, Mat LpM, 
    Mat Qc, Mat Qi,
    Mat H1, Mat H2,
    Mat J,
    //Mat X,
    Block_Rack* br_system,
    Dirichlet_BC* bc_T_1_H1_s,
    Dirichlet_BC* bc_T_1_H1_v,
    Dirichlet_BC* bc_O_o,
    Dirichlet_BC* bc_O_i,
    Dirichlet_BC* global_H1,
    PetscReal rtol = 1e-10, PetscInt max_iters = PETSC_DEFAULT, bool enable_monitor = false);


}