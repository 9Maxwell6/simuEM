#pragma once

#include "math/fem/block_rack.h"
#include "math/fem/bc_dirichlet.h"


namespace simu::T_Omega_AMS_c {

typedef struct AMS_Context 
{
    Mat P, G, I, LpM, H, J;
    Mat Qc, Qi;
    //Mat X;

    KSP inner_LpM_ksp, inner_Q_ksp;     // built once, reused per PCApply
    Vec tmp_1, tmp_2;
    Vec rho,   gamma;                 // size 3*nV (L-space RHS / sol)
    Vec zeta,  kappa;                 // size   nV (Q-space RHS / sol)
    Block_Rack* br_system;

    Vec zeta_1,  kappa_1;                 // size   nV (Q-space RHS / sol)
    Vec zeta_2,  kappa_2;                 // size   nV (Q-space RHS / sol)

    Vec rho_1, gamma_1;
    Vec rho_2, gamma_2;


    KSP inner_Qc_ksp, inner_Qi_ksp;     // built once, reused per PCApply
    //KSP inner_X_ksp; 

    Dirichlet_BC* bc_T_1_H1_s;
    Dirichlet_BC* bc_T_1_H1_v;
    Dirichlet_BC* bc_O_o;


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
    Mat Qc, Mat Qi, Mat H, Mat J,
    //Mat X,
    Block_Rack* br_system,
    Dirichlet_BC* bc_T_1_H1_s,
    Dirichlet_BC* bc_T_1_H1_v,
    Dirichlet_BC* bc_O_o,
    PetscReal rtol = 1e-10, PetscInt max_iters = PETSC_DEFAULT, bool enable_monitor = false);


}