#include "physics/electromagnetism/eddy_current/auxiliary_solver_c.h"


using namespace simu;


static inline void solve_2x2(const PetscScalar a, const PetscScalar b,
                             const PetscScalar c, const PetscScalar d,
                             const PetscScalar sp, const PetscScalar sq,
                             PetscScalar &xp, PetscScalar &xq)
{
    const PetscScalar det = a * d - b * c;
    if (PetscAbsScalar(det) > PETSC_SMALL * (PetscAbsScalar(a * d) + PetscAbsScalar(b * c)) ||
        PetscAbsScalar(det) > 1e-300)
    {
        const PetscScalar inv = 1.0 / det;
        xp = ( d * sp - b * sq) * inv;
        xq = (-c * sp + a * sq) * inv;
    }
    else
    {
        xp = (a != 0.0) ? sp / a : 0.0;
        xq = (d != 0.0) ? sq / d : 0.0;
    }
}
 
// One or more symmetric (forward + backward) 2x2-block Gauss-Seidel sweeps
// on A x = b, x updated in place.
//
//   nT : size of one T field (Re(T) and Im(T) each have nT dofs)
//   nO : size of one O field
//
// Block k pairs the dofs
//   k < nT :  ( k,            nT + k )            = ( ReT_i, ImT_i )
//   k >= nT:  ( 2nT + j,      2nT + nO + j )      = ( ReO_j, ImO_j ),  j = k - nT
inline PetscErrorCode block_ssgs_sweep(Mat A, Vec b, Vec x,
                                       PetscInt nT, PetscInt nO,
                                       PetscInt nsweeps = 1)
{
    PetscFunctionBeginUser;
 
    PetscInt        n;
    const PetscInt *ia, *ja;
    PetscBool       done;
    PetscCall(MatGetRowIJ(A, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done));
    PetscCheck(done, PETSC_COMM_SELF, PETSC_ERR_SUP,
               "block_ssgs_sweep: MatGetRowIJ not supported; A must be MATSEQAIJ");
    PetscCheck(n == 2 * (nT + nO), PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
               "block_ssgs_sweep: matrix size %" PetscInt_FMT
               " != 2*(nT+nO) = %" PetscInt_FMT, n, 2 * (nT + nO));
 
    const PetscScalar *aa;
    PetscCall(MatSeqAIJGetArrayRead(A, &aa));
 
    const PetscScalar *bb;
    PetscScalar       *xx;
    PetscCall(VecGetArrayRead(b, &bb));
    PetscCall(VecGetArray(x, &xx));
 
    const PetscInt nb = nT + nO;
 
    // Relax one 2x2 block (p, q): Gauss-Seidel, uses current xx for all
    // off-block columns, then solves the block exactly.
    const auto relax = [&](const PetscInt p, const PetscInt q)
    {
        PetscScalar sp = bb[p], sq = bb[q];
        PetscScalar app = 0.0, apq = 0.0, aqp = 0.0, aqq = 0.0;
 
        for (PetscInt t = ia[p]; t < ia[p + 1]; ++t) {
            const PetscInt j = ja[t];
            if      (j == p) app = aa[t];
            else if (j == q) apq = aa[t];
            else             sp -= aa[t] * xx[j];
        }
        for (PetscInt t = ia[q]; t < ia[q + 1]; ++t) {
            const PetscInt j = ja[t];
            if      (j == q) aqq = aa[t];
            else if (j == p) aqp = aa[t];
            else             sq -= aa[t] * xx[j];
        }
        solve_2x2(app, apq, aqp, aqq, sp, sq, xx[p], xx[q]);
    };
 
    const auto pair = [&](const PetscInt k, PetscInt &p, PetscInt &q)
    {
        if (k < nT) { p = k;                  q = nT + k; }
        else        { const PetscInt j = k - nT;
                      p = 2 * nT + j;         q = 2 * nT + nO + j; }
    };
 
    for (PetscInt s = 0; s < nsweeps; ++s) {
        PetscInt p, q;
        for (PetscInt k = 0; k < nb; ++k)       { pair(k, p, q); relax(p, q); }  // forward
        for (PetscInt k = nb - 1; k >= 0; --k)  { pair(k, p, q); relax(p, q); }  // backward
    }
 
    PetscCall(VecRestoreArray(x, &xx));
    PetscCall(VecRestoreArrayRead(b, &bb));
    PetscCall(MatSeqAIJRestoreArrayRead(A, &aa));
    PetscCall(MatRestoreRowIJ(A, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done));
    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode T_Omega_AMS_c::AMS_apply_decoupled(PC pc, Vec r, Vec x)
{
    PetscFunctionBeginUser;
    
    //PetscCall(VecCopy(r, x));
    //PetscFunctionReturn(0);

    //*
    AMS_Context *ctx;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));

    const Mat A = ctx->br_system->get_lhs();


    // all the matrices are already encode the sign
    const Mat Kc_r  = ctx->br_system->get_block_lhs(0,0); 
    const Mat Mc_r  = ctx->br_system->get_block_lhs(0,1);
    const Mat Xc_r  = ctx->br_system->get_block_lhs(0,3);
    const Mat Mc_c  = ctx->br_system->get_block_lhs(1,0);
    const Mat Kc_c  = ctx->br_system->get_block_lhs(1,1); 
    const Mat Xc_c  = ctx->br_system->get_block_lhs(1,2);
    const Mat Xi_r  = ctx->br_system->get_block_lhs(2,0);
    const Mat Ki_r  = ctx->br_system->get_block_lhs(2,2);
    const Mat Xi_c  = ctx->br_system->get_block_lhs(3,1);
    const Mat Ki_c  = ctx->br_system->get_block_lhs(3,3);




 
    PetscCall(VecZeroEntries(x));

    // pre-smoother.  single Gauss-Seidel sweep
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));
    

    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), x, ctx->block_x);
    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), r, ctx->block_r);


    const Vec x1 = ctx->block_x[0];
    const Vec x2 = ctx->block_x[1];
    const Vec x3 = ctx->block_x[2];
    const Vec x4 = ctx->block_x[3];

    const Vec r1 = ctx->block_r[0];
    const Vec r2 = ctx->block_r[1];
    const Vec r3 = ctx->block_r[2];
    const Vec r4 = ctx->block_r[3];

    Vec rho_c, gamma_c;

    std::vector<Vec> rho_c_l = {ctx->rho_1, ctx->rho_2};
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, nullptr, rho_c_l.data(), &rho_c));

    std::vector<Vec> gamma_c_l = {ctx->gamma_1, ctx->gamma_2};
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, nullptr, gamma_c_l.data(), &gamma_c));

    
    for(int n=0; n<1; ++n){
        //  rho_1 = Pt * (r1 - Kc_r*x1 - Mc_r*x2 - Xc_r*x4)
        PetscCall(MatMult(Kc_r, x1, ctx->tmp_1));                   
        PetscCall(MatMultAdd(Mc_r, x2, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(MatMultAdd(Xc_r, x4, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(VecAYPX(ctx->tmp_1, -1., r1));                  
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_1, ctx->rho_1));

        ctx->bc_T_1_H1_v->apply_to_vec(ctx->rho_1);  
        PetscCall(VecSet(ctx->gamma_1, 0.));
        PetscCall(KSPSolve(ctx->inner_LpM_ksp, ctx->rho_1,  ctx->gamma_1));
        PetscCall(MatMultAdd(ctx->P, ctx->gamma_1, x1, x1));


        //  rho_2 = Pt * (r2 - Mc_c*x1 - Kc_c*x2 - Xc_c*x3)
        PetscCall(MatMult(Mc_c, x1, ctx->tmp_2));                   
        PetscCall(MatMultAdd(Kc_c, x2, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(MatMultAdd(Xc_c, x3, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(VecAYPX(ctx->tmp_2, -1., r2));                  
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_2, ctx->rho_2));

        //ctx->bc_T_1_H1_v->apply_to_vec(ctx->rho_1);  
        ctx->bc_T_1_H1_v->apply_to_vec(ctx->rho_2);  

        //PetscCall(VecSet(ctx->gamma_1, 0.));
        PetscCall(VecSet(ctx->gamma_2, 0.));
        //PetscCall(KSPSolve(ctx->inner_LpM_ksp, ctx->rho_1,  ctx->gamma_1));
        PetscCall(KSPSolve(ctx->inner_LpM_ksp, ctx->rho_2,  ctx->gamma_2));
        //PetscCall(KSPSolve(ctx->inner_LpM_ksp, rho_c, gamma_c));


        //    xt = xt + P*gamma_1
        //PetscCall(MatMultAdd(ctx->P, ctx->gamma_1, x1, x1));
        PetscCall(MatMultAdd(ctx->P, ctx->gamma_2, x2, x2));

        //  zeta_1 = Gt * (r1 - Kc_r*x1 - Mc_r*x2 - Xc_r*x4)
        PetscCall(MatMult(Kc_r, x1, ctx->tmp_1));                   
        PetscCall(MatMultAdd(Mc_r, x2, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(MatMultAdd(Xc_r, x4, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(VecAYPX(ctx->tmp_1, -1., r1));                  
        PetscCall(MatMultTranspose(ctx->G, ctx->tmp_1, ctx->zeta_1));

        PetscCall(VecSet(ctx->kappa_1, 0.));
        ctx->bc_T_1_H1_s->apply_to_vec(ctx->zeta_1);
        PetscCall(KSPSolve(ctx->inner_Qc_ksp, ctx->zeta_1, ctx->kappa_1));
        PetscCall(MatMultAdd(ctx->G, ctx->kappa_1, x1, x1));

        //  zeta_2 = Gt * (r2 - Mc_c*x1 - Kc_c*x2 - Xc_c*x3)
        PetscCall(MatMult(Mc_c, x1, ctx->tmp_2));                   
        PetscCall(MatMultAdd(Kc_c, x2, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(MatMultAdd(Xc_c, x3, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(VecAYPX(ctx->tmp_2, -1., r2));                  
        PetscCall(MatMultTranspose(ctx->G, ctx->tmp_2, ctx->zeta_2));

        PetscCall(VecSet(ctx->kappa_2, 0.));
        ctx->bc_T_1_H1_s->apply_to_vec(ctx->zeta_2);
        PetscCall(KSPSolve(ctx->inner_Qc_ksp, ctx->zeta_2, ctx->kappa_2));
        PetscCall(MatMultAdd(ctx->G, ctx->kappa_2, x2, x2));


        //  zeta_3 = It * (r3 - Xi_r*x1 - Ki_r*x3)
        PetscCall(MatMult(Xi_r, x1, ctx->tmp_3));                   
        PetscCall(MatMultAdd(Ki_r, x3, ctx->tmp_3, ctx->tmp_3));   
        PetscCall(VecAYPX(ctx->tmp_3, -1., r3));                  
        PetscCall(MatMultTranspose(ctx->I, ctx->tmp_3, ctx->zeta_3));

        PetscCall(VecSet(ctx->kappa_3, 0.));
        ctx->bc_O_o->apply_to_vec(ctx->zeta_3);
        ctx->bc_O_i->apply_to_vec(ctx->zeta_3);
        PetscCall(KSPSolve(ctx->inner_Qi_ksp, ctx->zeta_3, ctx->kappa_3));
        PetscCall(MatMultAdd(ctx->I, ctx->kappa_3, x3, x3));

        //  zeta_4 = It * (r3 - Xi_c*x2 - Ki_c*x4)
        PetscCall(MatMult(Xi_c, x2, ctx->tmp_4));                   
        PetscCall(MatMultAdd(Ki_c, x4, ctx->tmp_4, ctx->tmp_4));   
        PetscCall(VecAYPX(ctx->tmp_4, -1., r4));                  
        PetscCall(MatMultTranspose(ctx->I, ctx->tmp_4, ctx->zeta_4));


        //PetscCall(VecSet(ctx->kappa_1, 0.));
        //PetscCall(VecSet(ctx->kappa_2, 0.));
        //PetscCall(VecSet(ctx->kappa_3, 0.));
        PetscCall(VecSet(ctx->kappa_4, 0.));


        //ctx->bc_T_1_H1_s->apply_to_vec(ctx->zeta_1);
        //ctx->bc_T_1_H1_s->apply_to_vec(ctx->zeta_2);

        //ctx->bc_O_o->apply_to_vec(ctx->zeta_3);
        //ctx->bc_O_i->apply_to_vec(ctx->zeta_3);
        ctx->bc_O_o->apply_to_vec(ctx->zeta_4);
        ctx->bc_O_i->apply_to_vec(ctx->zeta_4);
        
        //PetscCall(KSPSolve(ctx->inner_Qc_ksp, ctx->zeta_1, ctx->kappa_1));
        //PetscCall(KSPSolve(ctx->inner_Qc_ksp, ctx->zeta_2, ctx->kappa_2));

        //PetscCall(KSPSolve(ctx->inner_Qi_ksp, ctx->zeta_3, ctx->kappa_3));
        PetscCall(KSPSolve(ctx->inner_Qi_ksp, ctx->zeta_4, ctx->kappa_4));

        //PetscCall(MatMultAdd(ctx->G, ctx->kappa_1, x1, x1));
        //PetscCall(MatMultAdd(ctx->G, ctx->kappa_2, x2, x2));
        //PetscCall(MatMultAdd(ctx->I, ctx->kappa_3, x3, x3));
        PetscCall(MatMultAdd(ctx->I, ctx->kappa_4, x4, x4));
    }

    la_kernel::write_to_vec(ctx->block_x, x);
    // post-smoother (continues from updated z)
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    for (auto& b : ctx->block_x) PetscCall(VecDestroy(&b));
    for (auto& b : ctx->block_r) PetscCall(VecDestroy(&b));
    //PetscCall(VecDestroy(&rho_c));
    //PetscCall(VecDestroy(&gamma_c));
    //*/
    

    PetscFunctionReturn(0);
}


PetscErrorCode T_Omega_AMS_c::AMS_apply_global(PC pc, Vec r, Vec x)
{
    PetscFunctionBeginUser;
    AMS_Context *ctx;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));

    const Mat A = ctx->br_system->get_lhs();


    // all the matrices are already encode the sign
    const Mat Kc_r  = ctx->br_system->get_block_lhs(0,0); 
    const Mat Mc_r  = ctx->br_system->get_block_lhs(0,1);
    const Mat Xc_r  = ctx->br_system->get_block_lhs(0,3);
    const Mat Mc_c  = ctx->br_system->get_block_lhs(1,0);
    const Mat Kc_c  = ctx->br_system->get_block_lhs(1,1); 
    const Mat Xc_c  = ctx->br_system->get_block_lhs(1,2);
    const Mat Xi_r  = ctx->br_system->get_block_lhs(2,0);
    const Mat Ki_r  = ctx->br_system->get_block_lhs(2,2);
    const Mat Xi_c  = ctx->br_system->get_block_lhs(3,1);
    const Mat Ki_c  = ctx->br_system->get_block_lhs(3,3);




 
    PetscCall(VecZeroEntries(x));

    // pre-smoother.  single Gauss-Seidel sweep
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));
    

    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), x, ctx->block_x);
    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), r, ctx->block_r);


    const Vec x1 = ctx->block_x[0];
    const Vec x2 = ctx->block_x[1];
    const Vec x3 = ctx->block_x[2];
    const Vec x4 = ctx->block_x[3];

    const Vec r1 = ctx->block_r[0];
    const Vec r2 = ctx->block_r[1];
    const Vec r3 = ctx->block_r[2];
    const Vec r4 = ctx->block_r[3];

    Vec rho_c, gamma_c;

    //std::vector<Vec> rho_c_l = {ctx->rho_1, ctx->rho_2};
    //PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, nullptr, rho_c_l.data(), &rho_c));

    //std::vector<Vec> gamma_c_l = {ctx->gamma_1, ctx->gamma_2};
    //PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, nullptr, gamma_c_l.data(), &gamma_c));

    
    for(int n=0; n<1; ++n){
        //  rho_1 = Pt * (r1 - Kc_r*x1 - Mc_r*x2 - Xc_r*x4)
        PetscCall(MatMult(Kc_r, x1, ctx->tmp_1));                   
        PetscCall(MatMultAdd(Mc_r, x2, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(MatMultAdd(Xc_r, x4, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(VecAYPX(ctx->tmp_1, -1., r1));                  
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_1, ctx->rho_1));

        //  rho_2 = Pt * (r2 - Mc_c*x1 - Kc_c*x2 - Xc_c*x3)
        PetscCall(MatMult(Mc_c, x1, ctx->tmp_2));                   
        PetscCall(MatMultAdd(Kc_c, x2, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(MatMultAdd(Xc_c, x3, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(VecAYPX(ctx->tmp_2, -1., r2));                  
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_2, ctx->rho_2));

        ctx->bc_T_1_H1_v->apply_to_vec(ctx->rho_1);  
        ctx->bc_T_1_H1_v->apply_to_vec(ctx->rho_2);  

        PetscCall(VecSet(ctx->gamma_1, 0.));
        PetscCall(VecSet(ctx->gamma_2, 0.));
        PetscCall(KSPSolve(ctx->inner_LpM_ksp, ctx->rho_1,  ctx->gamma_1));
        PetscCall(KSPSolve(ctx->inner_LpM_ksp, ctx->rho_2,  ctx->gamma_2));
        //PetscCall(KSPSolve(ctx->inner_LpM_ksp, rho_c, gamma_c));


        //    xt = xt + P*gamma_1
        PetscCall(MatMultAdd(ctx->P, ctx->gamma_1, x1, x1));
        PetscCall(MatMultAdd(ctx->P, ctx->gamma_2, x2, x2));

        //  zeta_1 = Gt * (r1 - Kc_r*x1 - Mc_r*x2 - Xc_r*x4)
        PetscCall(MatMult(Kc_r, x1, ctx->tmp_1));                   
        PetscCall(MatMultAdd(Mc_r, x2, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(MatMultAdd(Xc_r, x4, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(VecAYPX(ctx->tmp_1, -1., r1));                  
        PetscCall(MatMultTranspose(ctx->G, ctx->tmp_1, ctx->zeta));

        PetscCall(VecScale(ctx->zeta, -1.));

        //  zeta_4 = It * (r4 - Xi_c*x2 - Ki_c*x4)
        PetscCall(MatMult(Xi_c, x2, ctx->tmp_4));                   
        PetscCall(MatMultAdd(Ki_c, x4, ctx->tmp_4, ctx->tmp_4));   
        PetscCall(VecAYPX(ctx->tmp_4, -1., r4));                  
        PetscCall(MatMultTransposeAdd(ctx->I, ctx->tmp_4, ctx->zeta, ctx->zeta));



        PetscCall(VecSet(ctx->kappa, 0.));
        ctx->bc_O_o->apply_to_vec(ctx->zeta);


        PetscCall(KSPSolve(ctx->inner_Q1_ksp, ctx->zeta, ctx->kappa));
        PetscCall(MatMultAdd(ctx->G, ctx->kappa, x2, x2));
        PetscCall(MatMultAdd(ctx->I, ctx->kappa, x4, x4));


        //  zeta_2 = Gt * (r2 - Mc_c*x1 - Kc_c*x2 - Xc_c*x3)
        PetscCall(MatMult(Mc_c, x1, ctx->tmp_2));                   
        PetscCall(MatMultAdd(Kc_c, x2, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(MatMultAdd(Xc_c, x3, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(VecAYPX(ctx->tmp_2, -1., r2));                  
        PetscCall(MatMultTranspose(ctx->G, ctx->tmp_2, ctx->zeta));

        //  zeta_3 = It * (r3 - Xi_r*x1 - Ki_r*x3)
        PetscCall(MatMult(Xi_r, x1, ctx->tmp_3));                   
        PetscCall(MatMultAdd(Ki_r, x3, ctx->tmp_3, ctx->tmp_3));   
        PetscCall(VecAYPX(ctx->tmp_3, -1., r3));                  
        PetscCall(MatMultTransposeAdd(ctx->I, ctx->tmp_3, ctx->zeta, ctx->zeta));

        PetscCall(VecSet(ctx->kappa, 0.));
        ctx->bc_O_o->apply_to_vec(ctx->zeta);

        PetscCall(KSPSolve(ctx->inner_Q1_ksp, ctx->zeta, ctx->kappa));
        
        PetscCall(MatMultAdd(ctx->G, ctx->kappa, x1, x1));
        PetscCall(MatMultAdd(ctx->I, ctx->kappa, x3, x3));


    }

    la_kernel::write_to_vec(ctx->block_x, x);
    // post-smoother (continues from updated z)
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    for (auto& b : ctx->block_x) PetscCall(VecDestroy(&b));
    for (auto& b : ctx->block_r) PetscCall(VecDestroy(&b));
    PetscFunctionReturn(0);
}


PetscErrorCode T_Omega_AMS_c::AMS_apply_coupled(PC pc, Vec r, Vec x)
{
    PetscFunctionBeginUser;

    AMS_Context *ctx;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));

    const Mat A = ctx->br_system->get_lhs();


    // all the matrices are already encode the sign
    const Mat Kc_r  = ctx->br_system->get_block_lhs(0,0); 
    const Mat Mc_r  = ctx->br_system->get_block_lhs(0,1);
    const Mat Xc_r  = ctx->br_system->get_block_lhs(0,3);
    const Mat Mc_c  = ctx->br_system->get_block_lhs(1,0);
    const Mat Kc_c  = ctx->br_system->get_block_lhs(1,1); 
    const Mat Xc_c  = ctx->br_system->get_block_lhs(1,2);
    const Mat Xi_r  = ctx->br_system->get_block_lhs(2,0);
    const Mat Ki_r  = ctx->br_system->get_block_lhs(2,2);
    const Mat Xi_c  = ctx->br_system->get_block_lhs(3,1);
    const Mat Ki_c  = ctx->br_system->get_block_lhs(3,3);


    const auto&    row_sizes = ctx->br_system->get_local_row_size();
    const PetscInt nT = row_sizes[0];   // = row_sizes[1]
    const PetscInt nO = row_sizes[2];   // = row_sizes[3]

 
    PetscCall(VecZeroEntries(x));

    // pre-smoother.  single Gauss-Seidel sweep
    //PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));
    PetscCall(block_ssgs_sweep(A, r, x, nT, nO, 1));
    
    

    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), x, ctx->block_x);
    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), r, ctx->block_r);


    const Vec x1 = ctx->block_x[0];
    const Vec x2 = ctx->block_x[1];
    const Vec x3 = ctx->block_x[2];
    const Vec x4 = ctx->block_x[3];

    const Vec r1 = ctx->block_r[0];
    const Vec r2 = ctx->block_r[1];
    const Vec r3 = ctx->block_r[2];
    const Vec r4 = ctx->block_r[3];

    Vec zeta_a, kappa_a;
    Vec zeta_b, kappa_b;


    Vec zeta_la[2] = {ctx->zeta_1, ctx->zeta_4};
    Vec kappa_la[2] = {ctx->kappa_2, ctx->kappa_4};

    Vec zeta_lb[2] = {ctx->zeta_2, ctx->zeta_3};
    Vec kappa_lb[2] = {ctx->kappa_1, ctx->kappa_3};

    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, zeta_la,  &zeta_a));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, kappa_la, &kappa_a));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, zeta_lb,  &zeta_b));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, kappa_lb, &kappa_b));

    //std::vector<Vec> rho_c_l = {ctx->rho_1, ctx->rho_2};
    //PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, nullptr, rho_c_l.data(), &rho_c));

    //std::vector<Vec> gamma_c_l = {ctx->gamma_1, ctx->gamma_2};
    //PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, nullptr, gamma_c_l.data(), &gamma_c));

    
    for(int n=0; n<1; ++n){
        //  rho_1 = Pt * (r1 - Kc_r*x1 - Mc_r*x2 - Xc_r*x4)
        PetscCall(MatMult(Kc_r, x1, ctx->tmp_1));                   
        PetscCall(MatMultAdd(Mc_r, x2, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(MatMultAdd(Xc_r, x4, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(VecAYPX(ctx->tmp_1, -1., r1));                  
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_1, ctx->rho_1));

        //ctx->bc_T_1_H1_v->apply_to_vec(ctx->rho_1);  
        //PetscCall(VecSet(ctx->gamma_1, 0.));
        //PetscCall(KSPSolve(ctx->inner_LpM_ksp, ctx->rho_1,  ctx->gamma_1));
        //    xt = xt + P*gamma_1
        //PetscCall(MatMultAdd(ctx->P, ctx->gamma_1, x1, x1));

        //  rho_2 = Pt * (r2 - Mc_c*x1 - Kc_c*x2 - Xc_c*x3)
        PetscCall(MatMult(Mc_c, x1, ctx->tmp_2));                   
        PetscCall(MatMultAdd(Kc_c, x2, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(MatMultAdd(Xc_c, x3, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(VecAYPX(ctx->tmp_2, -1., r2));                  
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_2, ctx->rho_2));

        ctx->bc_T_1_H1_v->apply_to_vec(ctx->rho_1);  
        ctx->bc_T_1_H1_v->apply_to_vec(ctx->rho_2);  

        PetscCall(VecSet(ctx->gamma_1, 0.));
        PetscCall(VecSet(ctx->gamma_2, 0.));
        PetscCall(KSPSolve(ctx->inner_LpM_ksp, ctx->rho_1,  ctx->gamma_1));
        PetscCall(KSPSolve(ctx->inner_LpM_ksp, ctx->rho_2,  ctx->gamma_2));
        //PetscCall(KSPSolve(ctx->inner_LpM_ksp, rho_c, gamma_c));


        //    xt = xt + P*gamma_2
        PetscCall(MatMultAdd(ctx->P, ctx->gamma_1, x1, x1));
        PetscCall(MatMultAdd(ctx->P, ctx->gamma_2, x2, x2));

        //  zeta_1 = Gt * (r1 - Kc_r*x1 - Mc_r*x2 - Xc_r*x4)
        PetscCall(MatMult(Kc_r, x1, ctx->tmp_1));                   
        PetscCall(MatMultAdd(Mc_r, x2, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(MatMultAdd(Xc_r, x4, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(VecAYPX(ctx->tmp_1, -1., r1));                  
        PetscCall(MatMultTranspose(ctx->G, ctx->tmp_1, ctx->zeta_1));

        // first row of 2x2 block system is applied with minus sign
        PetscCall(VecScale(ctx->zeta_1, -1.));


        //  zeta_4 = It * (r4 - Xi_c*x2 - Ki_c*x4)
        PetscCall(MatMult(Xi_c, x2, ctx->tmp_4));                   
        PetscCall(MatMultAdd(Ki_c, x4, ctx->tmp_4, ctx->tmp_4));   
        PetscCall(VecAYPX(ctx->tmp_4, -1., r4));                  
        PetscCall(MatMultTranspose(ctx->I, ctx->tmp_4, ctx->zeta_4));



        PetscCall(VecSet(kappa_a, 0.));
        ctx->bc_T_1_H1_s->apply_to_vec(ctx->zeta_1);
        ctx->bc_O_o->apply_to_vec(ctx->zeta_4);
        ctx->bc_O_i->apply_to_vec(ctx->zeta_4);


        //PetscCall(KSPSolve(ctx->inner_Q1_ksp, zeta_a, kappa_a));
        PetscCall(KSPSolve(ctx->inner_Q2_ksp, zeta_a, kappa_a));
        PetscCall(MatMultAdd(ctx->G, ctx->kappa_2, x2, x2));
        PetscCall(MatMultAdd(ctx->I, ctx->kappa_4, x4, x4));


        //  zeta_2 = Gt * (r2 - Mc_c*x1 - Kc_c*x2 - Xc_c*x3)
        PetscCall(MatMult(Mc_c, x1, ctx->tmp_2));                   
        PetscCall(MatMultAdd(Kc_c, x2, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(MatMultAdd(Xc_c, x3, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(VecAYPX(ctx->tmp_2, -1., r2));                  
        PetscCall(MatMultTranspose(ctx->G, ctx->tmp_2, ctx->zeta_2));

        //  zeta_3 = It * (r3 - Xi_r*x1 - Ki_r*x3)
        PetscCall(MatMult(Xi_r, x1, ctx->tmp_3));                   
        PetscCall(MatMultAdd(Ki_r, x3, ctx->tmp_3, ctx->tmp_3));   
        PetscCall(VecAYPX(ctx->tmp_3, -1., r3));                  
        PetscCall(MatMultTranspose(ctx->I, ctx->tmp_3, ctx->zeta_3));

        PetscCall(VecSet(kappa_b, 0.));
        ctx->bc_T_1_H1_s->apply_to_vec(ctx->zeta_2);
        ctx->bc_O_o->apply_to_vec(ctx->zeta_3);
        ctx->bc_O_i->apply_to_vec(ctx->zeta_3);

        PetscCall(KSPSolve(ctx->inner_Q2_ksp, zeta_b, kappa_b));
        
        PetscCall(MatMultAdd(ctx->G, ctx->kappa_1, x1, x1));
        PetscCall(MatMultAdd(ctx->I, ctx->kappa_3, x3, x3));


    }

    la_kernel::write_to_vec(ctx->block_x, x);
    // post-smoother (continues from updated z)
    //PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));
    PetscCall(block_ssgs_sweep(A, r, x, nT, nO, /*nsweeps=*/1));

    

    for (auto& b : ctx->block_x) PetscCall(VecDestroy(&b));
    for (auto& b : ctx->block_r) PetscCall(VecDestroy(&b));
    PetscFunctionReturn(0);
}



PetscErrorCode T_Omega_AMS_c::AMS_apply_fully_coupled(PC pc, Vec r, Vec x)
{

    PetscFunctionReturn(0);
}


PetscErrorCode T_Omega_AMS_c::AMS_destroy(PC pc)
{
    AMS_Context *ctx;
    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));
    PetscCall(KSPDestroy(&ctx->inner_LpM_ksp));
    PetscCall(KSPDestroy(&ctx->inner_Q1_ksp));
    PetscCall(KSPDestroy(&ctx->inner_Q2_ksp));
    PetscCall(KSPDestroy(&ctx->inner_Qc_ksp));
    PetscCall(KSPDestroy(&ctx->inner_Qi_ksp));
    PetscCall(VecDestroy(&ctx->tmp_1));
    PetscCall(VecDestroy(&ctx->tmp_2));
    PetscCall(VecDestroy(&ctx->rho));
    PetscCall(VecDestroy(&ctx->gamma));
    PetscCall(VecDestroy(&ctx->zeta));
    PetscCall(VecDestroy(&ctx->kappa));
    PetscCall(VecDestroy(&ctx->zeta_1));
    PetscCall(VecDestroy(&ctx->zeta_2));
    PetscCall(VecDestroy(&ctx->kappa_1));
    PetscCall(VecDestroy(&ctx->kappa_2));
    PetscCall(PetscFree(ctx));
    PetscFunctionReturn(0);
}




PetscErrorCode T_Omega_AMS_c::solve_AMS(
    AMS_Info& AMS_info,
    int algorithm_id,
    Mat P, Mat G, Mat I, Mat LpM, 
    Mat Qc, Mat Qi, 
    Mat H1, Mat H2,
    Mat J1, Mat J2,
    //Mat X,
    Block_Rack* br_system,
    Dirichlet_BC* bc_T_1_H1_s,
    Dirichlet_BC* bc_T_1_H1_v,
    Dirichlet_BC* bc_O_o,
    Dirichlet_BC* bc_O_i,
    Dirichlet_BC* global_H1,
    PetscReal rtol, PetscInt max_iters, bool enable_monitor)
{
    AMS_Context *ctx;
    KSP outer;
    PC  outer_pc;
    PetscFunctionBeginUser;

    

    // ---------- Allocate context ----------
    PetscCall(PetscNew(&ctx));
    ctx->P = P;
    ctx->G = G;
    ctx->I = I;
    ctx->LpM = LpM;
    ctx->Qc = Qc;
    ctx->Qi = Qi;
    ctx->H1 = H1;
    ctx->H2 = H2;
    ctx->J1 = J1;
    ctx->J2 = J2;
    ctx->br_system = br_system;
    ctx->bc_T_1_H1_s = bc_T_1_H1_s;
    ctx->bc_T_1_H1_v = bc_T_1_H1_v;
    ctx->bc_O_o = bc_O_o;
    ctx->bc_O_i = bc_O_i;

    PetscInt N_Vcycles = 1;  // number of V-cycles for Q

    // ---------- Inner KSP for L: PREONLY + BoomerAMG (single V-cycle) ----------
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_LpM_ksp));
    PetscCall(KSPSetType(ctx->inner_LpM_ksp, KSPPREONLY));  
    {
        PC pc_in;
        PetscCall(KSPGetPC(ctx->inner_LpM_ksp, &pc_in));
        PetscCall(PCSetType(pc_in, PCHYPRE));
        PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
        //PetscCall(PCSetType(pc_in, PCGAMG));
        //PetscCall(PCGAMGSetType(pc_in, PCGAMGAGG));
    }
    PetscCall(KSPSetOptionsPrefix(ctx->inner_LpM_ksp, "inner_LpM_"));
    //PetscCall(PetscOptionsSetValue(NULL, "-inner_LpM_pc_mg_cycle_type", "v"));
    //PetscCall(PetscOptionsSetValue(NULL, "-inner_LpM_pc_mg_type", "multiplicative"));

    char buf[32];
    PetscCall(PetscSNPrintf(buf, sizeof(buf), "%d", (int)N_Vcycles));
    //PetscCall(PetscOptionsSetValue(NULL, "-inner_LpM_pc_hypre_boomeramg_max_iter", buf));
    //PetscCall(PetscOptionsSetValue(NULL, "-inner_LpM_pc_hypre_boomeramg_tol", "0.0")); // 0 = never stop early


    PetscCall(KSPSetOperators(ctx->inner_LpM_ksp, LpM, LpM));
    PetscCall(KSPSetFromOptions(ctx->inner_LpM_ksp));
    PetscCall(KSPSetUp(ctx->inner_LpM_ksp)); 
    



    if(algorithm_id == 1){
        // ---------- Inner KSP for [Qc] ----------
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Qc_ksp));
        PetscCall(KSPSetType(ctx->inner_Qc_ksp, KSPPREONLY));
        {
            PC pc_in;
            PetscCall(KSPGetPC(ctx->inner_Qc_ksp, &pc_in));
            PetscCall(PCSetType(pc_in, PCHYPRE));
            PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
        }
        PetscCall(KSPSetOptionsPrefix(ctx->inner_Qc_ksp, "inner_Qc_"));

        // fixed number of V-cycles, done INSIDE hypre:
        //char buf[32];
        //PetscCall(PetscSNPrintf(buf, sizeof(buf), "%d", (int)N_Vcycles));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Qc_pc_hypre_boomeramg_max_iter", buf));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Qc_pc_hypre_boomeramg_tol", "0.0")); // 0 = never stop early

        PetscCall(KSPSetOperators(ctx->inner_Qc_ksp, Qc, Qc));
        PetscCall(KSPSetFromOptions(ctx->inner_Qc_ksp));
        PetscCall(KSPSetUp(ctx->inner_Qc_ksp));


        // ---------- Inner KSP for [Qi] ----------
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Qi_ksp));
        PetscCall(KSPSetType(ctx->inner_Qi_ksp, KSPPREONLY));
        {
            PC pc_in;
            PetscCall(KSPGetPC(ctx->inner_Qi_ksp, &pc_in));
            PetscCall(PCSetType(pc_in, PCHYPRE));
            PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
        }
        PetscCall(KSPSetOptionsPrefix(ctx->inner_Qi_ksp, "inner_Qi_"));

        PetscCall(PetscOptionsSetValue(NULL, "-inner_Qi_pc_hypre_boomeramg_max_iter", buf));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Qi_pc_hypre_boomeramg_tol", "0.0")); // 0 = never stop early

        PetscCall(KSPSetOperators(ctx->inner_Qi_ksp, Qi, Qi));
        PetscCall(KSPSetFromOptions(ctx->inner_Qi_ksp));
        PetscCall(KSPSetUp(ctx->inner_Qi_ksp));

        // ---------- Initialize vectors ----------
        PetscCall(MatCreateVecs(P, &ctx->gamma_1, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_1, NULL));    // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->gamma_2, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_2, NULL));    // size 3 * #node in T-field

        PetscCall(MatCreateVecs(Qc, &ctx->kappa_1, &ctx->zeta_1)); 
        PetscCall(MatCreateVecs(Qc, &ctx->kappa_2, &ctx->zeta_2)); 
        PetscCall(MatCreateVecs(Qi, &ctx->kappa_3, &ctx->zeta_3)); 
        PetscCall(MatCreateVecs(Qi, &ctx->kappa_4, &ctx->zeta_4)); 

    }if(algorithm_id == 2){ // global
        // ---------- Inner KSP for [H1] ----------
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Q1_ksp));
        PetscCall(KSPSetType(ctx->inner_Q1_ksp, KSPPREONLY));
        {
            PC pc_in;
            PetscCall(KSPGetPC(ctx->inner_Q1_ksp, &pc_in));
            PetscCall(PCSetType(pc_in, PCHYPRE));
            PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
        }
        PetscCall(KSPSetOptionsPrefix(ctx->inner_Q1_ksp, "inner_Q1_"));

        // fixed number of V-cycles, done INSIDE hypre:
        //char buf[32];
        //PetscCall(PetscSNPrintf(buf, sizeof(buf), "%d", (int)N_Vcycles));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q1_pc_hypre_boomeramg_max_iter", buf));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q1_pc_hypre_boomeramg_tol", "0.0")); // 0 = never stop early
        PetscCall(KSPSetOperators(ctx->inner_Q1_ksp, H1, H1));
        PetscCall(KSPSetFromOptions(ctx->inner_Q1_ksp));
        PetscCall(KSPSetUp(ctx->inner_Q1_ksp));

        // ---------- Inner KSP for [H2] ----------
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Q2_ksp));
        PetscCall(KSPSetType(ctx->inner_Q2_ksp, KSPPREONLY));
        {
            PC pc_in;
            PetscCall(KSPGetPC(ctx->inner_Q2_ksp, &pc_in));
            PetscCall(PCSetType(pc_in, PCHYPRE));
            PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
        }
        PetscCall(KSPSetOptionsPrefix(ctx->inner_Q2_ksp, "inner_Q2_"));

        // fixed number of V-cycles, done INSIDE hypre:
        //PetscCall(PetscSNPrintf(buf, sizeof(buf), "%d", (int)N_Vcycles));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q1_pc_hypre_boomeramg_max_iter", buf));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q1_pc_hypre_boomeramg_tol", "0.0")); // 0 = never stop early
        PetscCall(KSPSetOperators(ctx->inner_Q2_ksp, H2, H2));
        PetscCall(KSPSetFromOptions(ctx->inner_Q2_ksp));
        PetscCall(KSPSetUp(ctx->inner_Q2_ksp));

        // ---------- Initialize vectors ----------
        PetscCall(MatCreateVecs(P, &ctx->gamma_1, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_1, NULL));    // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->gamma_2, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_2, NULL));    // size 3 * #node in T-field
        PetscCall(MatCreateVecs(H1, &ctx->kappa, &ctx->zeta)); // size     #node in global field

    }else if(algorithm_id == 3){
        /*
        // ---------- Inner KSP for [J] ----------
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Q1_ksp));
        PetscCall(KSPSetType(ctx->inner_Q1_ksp, KSPPREONLY));
        {
            PC pc_in;
            PetscCall(KSPGetPC(ctx->inner_Q1_ksp, &pc_in));
            PetscCall(PCSetType(pc_in, PCHYPRE));
            PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
            //PetscCall(PCSetType(pc_in, PCGAMG));
            //PetscCall(PCGAMGSetType(pc_in, PCGAMGAGG));
        }
        PetscCall(KSPSetOptionsPrefix(ctx->inner_Q1_ksp, "inner_Q1_"));

        // fixed number of V-cycles, done INSIDE hypre:
        PetscCall(PetscSNPrintf(buf, sizeof(buf), "%d", (int)N_Vcycles));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q1_pc_hypre_boomeramg_max_iter", buf));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q1_pc_hypre_boomeramg_tol", "0.0")); // 0 = never stop early
    

        //PetscCall(PetscOptionsSetValue(NULL, "-inner_Q1_pc_mg_cycle_type", "v"));
        //PetscCall(PetscOptionsSetValue(NULL, "-inner_Q1_pc_mg_type", "multiplicative"));


        //PetscCall(KSPSetOperators(ctx->inner_Q1_ksp, J1, J1));
        PetscCall(KSPSetOperators(ctx->inner_Q1_ksp, J2, J2));
        PetscCall(KSPSetFromOptions(ctx->inner_Q1_ksp));
        PetscCall(KSPSetUp(ctx->inner_Q1_ksp));

        */


        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Q2_ksp));
        PetscCall(KSPSetType(ctx->inner_Q2_ksp, KSPPREONLY));
        {
            PC pc_in;
            PetscCall(KSPGetPC(ctx->inner_Q2_ksp, &pc_in));
            PetscCall(PCSetType(pc_in, PCHYPRE));
            PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
            //PetscCall(PCSetType(pc_in, PCGAMG));
            //PetscCall(PCGAMGSetType(pc_in, PCGAMGAGG));
        }
        PetscCall(KSPSetOptionsPrefix(ctx->inner_Q2_ksp, "inner_Q2_"));

        // fixed number of V-cycles, done INSIDE hypre:
        //PetscCall(PetscSNPrintf(buf, sizeof(buf), "%d", (int)N_Vcycles));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q2_pc_hypre_boomeramg_max_iter", buf));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q2_pc_hypre_boomeramg_tol", "0.0")); // 0 = never stop early

        //PetscCall(PetscOptionsSetValue(NULL, "-inner_Q2_pc_mg_cycle_type", "v"));
        //PetscCall(PetscOptionsSetValue(NULL, "-inner_Q2_pc_mg_type", "multiplicative"));

        PetscCall(KSPSetOperators(ctx->inner_Q2_ksp, J2, J2));
        PetscCall(KSPSetFromOptions(ctx->inner_Q2_ksp));
        PetscCall(KSPSetUp(ctx->inner_Q2_ksp));

        // ---------- Initialize vectors ----------
        PetscCall(MatCreateVecs(P, &ctx->gamma_1, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_1, NULL));    // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->gamma_2, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_2, NULL));    // size 3 * #node in T-field

        PetscCall(MatCreateVecs(G, &ctx->kappa_1, NULL)); 
        PetscCall(MatCreateVecs(G, &ctx->zeta_1, NULL)); 
        PetscCall(MatCreateVecs(G, &ctx->kappa_2, NULL)); 
        PetscCall(MatCreateVecs(G, &ctx->zeta_2, NULL)); 

        PetscCall(MatCreateVecs(I, &ctx->kappa_3, NULL)); 
        PetscCall(MatCreateVecs(I, &ctx->zeta_3, NULL)); 
        PetscCall(MatCreateVecs(I, &ctx->kappa_4, NULL)); 
        PetscCall(MatCreateVecs(I, &ctx->zeta_4, NULL)); 


        
    }

 
    // ---------- Initialize vectors ----------
    const Mat T   = ctx->br_system->get_block_lhs(0,0);
    const Mat O   = ctx->br_system->get_block_lhs(2,2);
    PetscCall(MatCreateVecs(T, NULL, &ctx->tmp_1));       // size #edge in T-field
    PetscCall(MatCreateVecs(T, NULL, &ctx->tmp_2));       // size #edge in T-field
    PetscCall(MatCreateVecs(O, NULL, &ctx->tmp_3));       // size #node in Omega-field
    PetscCall(MatCreateVecs(O, NULL, &ctx->tmp_4));       // size #node in Omega-field

 
    // ---------- Outer KSP + PCShell(AMS) ----------
    const Mat A = ctx->br_system->get_lhs();
    const Vec r = ctx->br_system->get_rhs();
    const Vec x = ctx->br_system->get_x();

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &outer));


    PetscCall(KSPSetType(outer, KSPGMRES));
    PetscCall(KSPGMRESSetRestart(outer, 1000));  


    PetscCall(KSPSetOperators(outer, A, A));
    PetscCall(KSPSetTolerances(outer, rtol, PETSC_DEFAULT, PETSC_DEFAULT, max_iters));

    PetscCall(KSPSetNormType(outer, KSP_NORM_UNPRECONDITIONED));  // true residual


    if (enable_monitor) 
    {
        PetscViewerAndFormat *vf;
        PetscCall(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,
                                            PETSC_VIEWER_DEFAULT, &vf));
        PetscCall(KSPMonitorSet(outer,
                                (KSPMonitorFn *)KSPMonitorTrueResidual,
                                vf,
                                (PetscCtxDestroyFn *)PetscViewerAndFormatDestroy));
    }
 
    PetscCall(KSPGetPC(outer, &outer_pc));
    PetscCall(PCSetType(outer_pc, PCSHELL));
    PetscCall(PCShellSetContext(outer_pc, ctx));
    switch (algorithm_id)
    {
    case 1: PetscCall(PCShellSetApply  (outer_pc, AMS_apply_decoupled)); break;
    case 2: PetscCall(PCShellSetApply  (outer_pc, AMS_apply_global)); break;
    case 3: PetscCall(PCShellSetApply  (outer_pc, AMS_apply_coupled)); break;
    case 4: PetscCall(PCShellSetApply  (outer_pc, AMS_apply_fully_coupled)); break;
    default:
        Logger::error("T_Omega_AMS_c::solve_AMS: unknown AMS algorithm id: "+std::to_string(algorithm_id));
        break;
    }
    //PetscCall(PCShellSetApply  (outer_pc, AMS_apply));
    PetscCall(PCShellSetDestroy(outer_pc, AMS_destroy));
    PetscCall(PCShellSetName   (outer_pc, "AMS"));
 
    PetscCall(KSPSetFromOptions(outer));

    PetscCall(KSPSetComputeSingularValues(outer, PETSC_TRUE));  // for computing condition number
 
    // ---------- Solve ----------
    PetscCall(VecSet(x, 0.0));
    PetscCall(KSPSolve(outer, r, x));
 
    // ---------- Report ----------
    KSPConvergedReason reason;
    PetscInt           iters;
    PetscReal          rnorm;
    PetscReal          smax, smin;
    PetscCall(KSPComputeExtremeSingularValues(outer, &smax, &smin));
    PetscCall(KSPGetConvergedReason(outer, &reason));
    PetscCall(KSPGetIterationNumber(outer, &iters));
    PetscCall(KSPGetResidualNorm   (outer, &rnorm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "[AMS] outer GMRES: iters=%" PetscInt_FMT
        ", final ||r||=%.3e, reason=%d\n",
        iters, (double)rnorm, (int)reason));
 
    PetscCall(KSPDestroy(&outer));   // also invokes AMS_destroy through PCShell

    AMS_info.n_iteration = iters;
    AMS_info.condition_number = smax / smin;

    PetscFunctionReturn(0);
}


