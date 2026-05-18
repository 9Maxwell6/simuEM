#pragma once

#include "math/fem/block_rack.h"
#include "math/fem/bc_dirichlet.h"


namespace simu::T_Omega_AMS {

typedef struct AMS_Context 
{
    Mat P, G, I, L, Q;
    KSP inner_L_ksp, inner_Q_ksp;     // built once, reused per PCApply
    Vec tmp_1, tmp_2;
    Vec rho,   gamma;                 // size 3*nV (L-space RHS / sol)
    Vec zeta,  kappa;                 // size   nV (Q-space RHS / sol)
    Block_Rack* br_system;
    Dirichlet_BC *bc_v;                   // for L (vector / 3D nodal space)
    Dirichlet_BC *bc_s;                   // for Q (scalar nodal space)
} AMS_Context;

typedef struct AMS_Info 
{
    PetscInt n_iteration = 0;
} AMS_Info;



PetscErrorCode AMS_apply(PC pc, Vec r, Vec x);



PetscErrorCode AMS_destroy(PC pc);




PetscErrorCode solve_AMS(
    AMS_Info& AMS_info,
    Mat P, Mat G, Mat I, Mat L, Mat Q,
    Block_Rack* br_system,
    Dirichlet_BC *bc_v, Dirichlet_BC *bc_s,
    PetscReal rtol = 1e-10, PetscInt max_iters = PETSC_DEFAULT, bool enable_monitor = false);


}