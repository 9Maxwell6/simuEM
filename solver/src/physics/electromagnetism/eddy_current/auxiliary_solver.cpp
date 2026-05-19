#include "physics/electromagnetism/eddy_current/auxiliary_solver.h"


using namespace simu;

PetscErrorCode T_Omega_AMS::AMS_apply(PC pc, Vec r, Vec x)
{
    AMS_Context *ctx;
    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));

    const Mat A = ctx->br_system->get_lhs();
    //const Vec br_x  = ctx->br_system->get_x();



    const Mat O   = ctx->br_system->get_block_lhs(0,0);
    const Mat O_c = ctx->br_system->get_block_lhs(0,1);
    const Mat T   = ctx->br_system->get_block_lhs(1,1);
    const Mat T_c = ctx->br_system->get_block_lhs(1,0);


    /**
     *  block system view:
     * 
     *      [   Omega    ] [ coupling_o ]  *  [x_o]  =  [r_o]
     *      [ coupling_t ] [     T      ]     [x_t]     [r_t]
     * 
     */
 
    PetscCall(VecZeroEntries(x));
    
    petsc_util::petsc_save_ascii_vec(x, "x_00.txt");

    // pre-smoother.  single Gauss-Seidel sweep
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    petsc_util::petsc_save_ascii_vec(x, "x_01.txt");

    std::vector<Vec> block_x(2);
    std::vector<Vec> block_r(2);
    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), x, block_x);
    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), r, block_r);

    const Vec xo = block_x[0];
    const Vec xt = block_x[1];

    const Vec ro = block_r[0];
    const Vec rt = block_r[1];
 
    // -- Step 2-3: residual after pre-smooth, then project to auxiliary spaces
    //    tmp_1 = rt - [T]*xt - [coupling_t]*xo
    PetscCall(MatMult(T, xt, ctx->tmp_1));                    // tmp_1 = [T]*xt
    PetscCall(MatMultAdd(T_c, xo, ctx->tmp_1, ctx->tmp_1));   // tmp_1 = tmp_1 + [coupling_t]*xo
    PetscCall(VecAYPX(ctx->tmp_1, -1., rt));                  // tmp_1 = rt - tmp_1

    //    tmp_2 = ro - [O]*xo - [coupling_o]*xt
    PetscCall(MatMult(O, xo, ctx->tmp_2));                    // tmp_2 = [O]*xo
    PetscCall(MatMultAdd(O_c, xt, ctx->tmp_2, ctx->tmp_2));   // tmp_2 = tmp_2 + [coupling_o]*xt
    PetscCall(VecAYPX(ctx->tmp_2, -1., ro));                  // tmp_2 = ro - tmp_2
 
    //    rho  = P^T * tmp_1   ((H1)^3)
    //    zeta = G^T * tmp_1   (H1)
    PetscCall(MatMultTranspose(ctx->P, ctx->tmp_1, ctx->rho));
    PetscCall(MatMultTranspose(ctx->G, ctx->tmp_1, ctx->zeta));

    //petsc_util::petsc_save_ascii_vec(ctx->zeta, "zeta_T0.txt");

    //ctx->bc_s3->apply_to_vec(ctx->zeta);

    //PetscCall(VecScale(ctx->zeta, 0.5));
    //PetscCall(VecScale(ctx->tmp_2, 0.5));



    //for(dof_idx id : ctx->bc_s3->bc_dofs)
    //    std::cout<<id<<std::endl;
    //std::cout<<"total T-BC size: "<<ctx->bc_s3->bc_dofs.size()<<std::endl;

    //petsc_util::petsc_save_ascii_vec(ctx->zeta, "zeta_T.txt");
    

    //PetscCall(MatMultTranspose(ctx->I, ctx->tmp_2, ctx->zeta));

    petsc_util::petsc_save_ascii_vec(ctx->zeta, "zeta_O.txt");

    //exit(0);

    //    zeta = zeta + I^T * tmp_2   (H1)
    PetscCall(MatMultTransposeAdd(ctx->I, ctx->tmp_2, ctx->zeta, ctx->zeta));

    
 
    // -- Step 4: enforce Dirichlet on inner-solve RHS and initial guess,
    //            then solve each inner system with one BoomerAMG V-cycle.
    PetscCall(VecSet(ctx->gamma, 0.));
    PetscCall(VecSet(ctx->kappa, 0.));

    petsc_util::petsc_save_ascii_vec(ctx->rho, "rho_00.txt");
 
    ctx->bc_v->apply_to_vec(ctx->rho);   // ??
    ctx->bc_s1->apply_to_vec(ctx->zeta);  // ??
    ctx->bc_s2->apply_to_vec(ctx->zeta);
    ctx->bc_s3->apply_to_vec(ctx->zeta);

    petsc_util::petsc_save_ascii_vec(rt, "rt.txt");

    petsc_util::petsc_save_ascii_vec(ctx->rho, "rho.txt");

    petsc_util::petsc_save_ascii_vec(ctx->zeta, "zeta.txt");
 
    //    Single BoomerAMG V-cycle for each (no inner GMRES; KSPPREONLY).
    PetscCall(KSPSolve(ctx->inner_L_ksp, ctx->rho,  ctx->gamma));
    PetscCall(KSPSolve(ctx->inner_Q_ksp, ctx->zeta, ctx->kappa));

    petsc_util::petsc_save_ascii_vec(ctx->gamma, "gamma.txt");
    petsc_util::petsc_save_ascii_vec(ctx->kappa, "kappa.txt");

    //ctx->bc_s1->apply_to_vec(ctx->kappa);
    //ctx->bc_s2->apply_to_vec(ctx->kappa);
    //ctx->bc_s3->apply_to_vec(ctx->kappa);

 
    // -- Step 5: transfer corrections back to E_h
    //    xt = xt + P*gamma_1 + G*kappa
    //    xt = xt + P*gamma_1
    PetscCall(MatMultAdd(ctx->P, ctx->gamma, xt, xt));
    //    xt = xt + G*kappa
    PetscCall(MatMultAdd(ctx->G, ctx->kappa, xt, xt));
    //    xo = xo + I*kappa
    PetscCall(MatMultAdd(ctx->I, ctx->kappa, xo, xo));

    petsc_util::petsc_save_ascii_vec(xo, "xo.txt");
    petsc_util::petsc_save_ascii_vec(xt, "xt.txt");

    /*
    PetscReal xt_norm, Gk_norm, kappa_norm, zeta_norm;
    PetscCall(VecNorm(xt, NORM_2, &xt_norm));
    PetscCall(VecNorm(ctx->kappa, NORM_2, &kappa_norm));
    PetscCall(VecNorm(ctx->zeta,  NORM_2, &zeta_norm));
    PetscCall(MatMult(ctx->G, ctx->kappa, ctx->tmp_1));  // borrow tmp_1
    PetscCall(VecNorm(ctx->tmp_1, NORM_2, &Gk_norm));
    PetscPrintf(PETSC_COMM_WORLD,
        "[AMS] ||xt||=%.3e  ||zeta||=%.3e  ||kappa||=%.3e  ||G*kappa||=%.3e\n",
        (double)xt_norm, (double)zeta_norm, (double)kappa_norm, (double)Gk_norm);
    //*/

    
    

    la_kernel::write_to_vec(block_x, x);

    petsc_util::petsc_save_ascii_vec(x, "x_02.txt");

 
    // -- Step 6: post-smoother (continues from updated z)
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    //PetscCall(VecCopy(x, br_x));
    
    petsc_util::petsc_save_ascii_vec(x, "x_03.txt");


    PetscFunctionReturn(0);
}



PetscErrorCode T_Omega_AMS::AMS_destroy(PC pc)
{
    AMS_Context *ctx;
    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));
    PetscCall(KSPDestroy(&ctx->inner_L_ksp));
    PetscCall(KSPDestroy(&ctx->inner_Q_ksp));
    PetscCall(VecDestroy(&ctx->tmp_1));
    PetscCall(VecDestroy(&ctx->tmp_2));
    PetscCall(VecDestroy(&ctx->rho));
    PetscCall(VecDestroy(&ctx->gamma));
    PetscCall(VecDestroy(&ctx->zeta));
    PetscCall(VecDestroy(&ctx->kappa));
    PetscCall(PetscFree(ctx));
    PetscFunctionReturn(0);
}




PetscErrorCode T_Omega_AMS::solve_AMS(
    AMS_Info& AMS_info,
    Mat P, Mat G, Mat I, Mat L, Mat Q,
    Block_Rack* br_system,
    Dirichlet_BC *bc_v, Dirichlet_BC *bc_s1, Dirichlet_BC *bc_s2, Dirichlet_BC *bc_s3,
    PetscReal rtol, PetscInt max_iters, bool enable_monitor)
{
    AMS_Context *ctx;
    KSP outer;
    PC  outer_pc;
    PetscFunctionBeginUser;

    bc_v->apply_to_mat(L);
    bc_s1->apply_to_mat(Q);
    bc_s2->apply_to_mat(Q);
    bc_s3->apply_to_mat(Q);

    petsc_util::petsc_save_ascii_mat(L, "L_initial.txt");
 
    // ---------- Allocate context ----------
    PetscCall(PetscNew(&ctx));
    ctx->P = P;
    ctx->G = G;
    ctx->I = I;
    ctx->L = L;
    ctx->Q = Q;
    ctx->br_system = br_system;
    ctx->bc_v = bc_v;
    ctx->bc_s1 = bc_s1;
    ctx->bc_s2 = bc_s2;
    ctx->bc_s3 = bc_s3;
 
    // ---------- Inner KSP for L: PREONLY + BoomerAMG (single V-cycle) ----------
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_L_ksp));
    PetscCall(KSPSetType(ctx->inner_L_ksp, KSPPREONLY));  
    {
        PC pc_in;
        PetscCall(KSPGetPC(ctx->inner_L_ksp, &pc_in));
        PetscCall(PCSetType(pc_in, PCHYPRE));
        PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
    }
    PetscCall(KSPSetOptionsPrefix(ctx->inner_L_ksp, "inner_L_"));
    PetscCall(KSPSetOperators(ctx->inner_L_ksp, L, L));
    PetscCall(KSPSetFromOptions(ctx->inner_L_ksp));
    PetscCall(KSPSetUp(ctx->inner_L_ksp)); 
 
    // ---------- Inner KSP for Q ----------
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Q_ksp));
    PetscCall(KSPSetType(ctx->inner_Q_ksp, KSPPREONLY));
    {
        PC pc_in;
        PetscCall(KSPGetPC(ctx->inner_Q_ksp, &pc_in));
        PetscCall(PCSetType(pc_in, PCHYPRE));
        PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
    }

    PetscCall(KSPSetOptionsPrefix(ctx->inner_Q_ksp, "inner_Q_"));
    PetscCall(KSPSetOperators(ctx->inner_Q_ksp, Q, Q));
    PetscCall(KSPSetFromOptions(ctx->inner_Q_ksp));
    PetscCall(KSPSetUp(ctx->inner_Q_ksp));
 
    // ---------- Scratch vectors ----------
    const Mat O   = ctx->br_system->get_block_lhs(0,0);
    const Mat T   = ctx->br_system->get_block_lhs(1,1);
    PetscCall(MatCreateVecs(T, NULL, &ctx->tmp_1));       // size #edge in T-field
    PetscCall(MatCreateVecs(O, NULL, &ctx->tmp_2));       // size #node in Omega-field
    PetscCall(MatCreateVecs(L, &ctx->gamma, &ctx->rho));  // size 3 * #node in T-field
    PetscCall(MatCreateVecs(Q, &ctx->kappa, &ctx->zeta)); // size     #node in global field
 
    // ---------- Outer KSP + PCShell(AMS) ----------
    const Mat A = ctx->br_system->get_lhs();
    const Vec r = ctx->br_system->get_rhs();
    const Vec x = ctx->br_system->get_x();

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &outer));
    PetscCall(KSPSetType(outer, KSPGMRES));     

    //PetscCall(KSPSetType(outer, KSPFGMRES));
    //PetscCall(KSPSetPCSide(outer, PC_RIGHT));
    //PetscCall(KSPSetNormType(outer, KSP_NORM_UNPRECONDITIONED));
    //PetscCall(KSPGMRESSetRestart(outer, 200));

    PetscCall(KSPSetOperators(outer, A, A));
    PetscCall(KSPSetTolerances(outer, rtol, PETSC_DEFAULT, PETSC_DEFAULT, max_iters));

    //PetscCall(KSPSetNormType(outer, KSP_NORM_UNPRECONDITIONED));  // true residual
    //PetscCall(KSPSetNormType(outer, KSP_NORM_DEFAULT));

    if(enable_monitor)
    {
        {
            PetscViewerAndFormat *vf;
            PetscCall(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,
                                                PETSC_VIEWER_DEFAULT, &vf));
            PetscCall(KSPMonitorSet(outer,
                (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*))KSPMonitorTrueResidual,
                vf,
                (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy)
            );
        }
    }
 
    PetscCall(KSPGetPC(outer, &outer_pc));
    PetscCall(PCSetType(outer_pc, PCSHELL));
    PetscCall(PCShellSetContext(outer_pc, ctx));
    PetscCall(PCShellSetApply  (outer_pc, AMS_apply));
    PetscCall(PCShellSetDestroy(outer_pc, AMS_destroy));
    PetscCall(PCShellSetName   (outer_pc, "AMS"));
 
    PetscCall(KSPSetFromOptions(outer));
 
    // ---------- Solve ----------
    PetscCall(VecSet(x, 0.0));
    PetscCall(KSPSolve(outer, r, x));
 
    // ---------- Report ----------
    KSPConvergedReason reason;
    PetscInt           iters;
    PetscReal          rnorm;
    PetscCall(KSPGetConvergedReason(outer, &reason));
    PetscCall(KSPGetIterationNumber(outer, &iters));
    PetscCall(KSPGetResidualNorm   (outer, &rnorm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "[AMS] outer FGMRES: iters=%" PetscInt_FMT
        ", final ||r||=%.3e, reason=%d\n",
        iters, (double)rnorm, (int)reason));
 
    PetscCall(KSPDestroy(&outer));   // also invokes AMS_destroy through PCShell

    AMS_info.n_iteration = iters;

    PetscFunctionReturn(0);
}


