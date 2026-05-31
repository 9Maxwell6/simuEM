#include "physics/electromagnetism/eddy_current/auxiliary_solver.h"


using namespace simu;

PetscErrorCode T_Omega_AMS::AMS_apply(PC pc, Vec r, Vec x)
{
    //petsc_util::petsc_save_ascii_vec(r, "r.txt");
    //petsc_util::petsc_save_ascii_vec(x, "x.txt");
    //exit(0);
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

    //PetscScalar value = 1.;
    //PetscInt    ix[1] = {0};

    // pre-smoother.  single Gauss-Seidel sweep
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    std::vector<Vec> block_x(2);
    std::vector<Vec> block_r(2);
    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), x, block_x);
    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), r, block_r);


    const Vec xo = block_x[0];
    const Vec xt = block_x[1];

    const Vec ro = block_r[0];
    const Vec rt = block_r[1];
    
    for(int n=0; n<100; ++n){

 
        //    tmp_1 = rt - [T]*xt - [coupling_t]*xo
        PetscCall(MatMult(T, xt, ctx->tmp_1));                    // tmp_1 = [T]*xt
        PetscCall(MatMultAdd(T_c, xo, ctx->tmp_1, ctx->tmp_1));   // tmp_1 = tmp_1 + [coupling_t]*xo
        PetscCall(VecAYPX(ctx->tmp_1, -1., rt));                  // tmp_1 = rt - tmp_1

        
        //    rho  = P^T * tmp_1   ((H1)^3)
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_1, ctx->rho));

        ctx->bc_v->apply_to_vec(ctx->rho);  
        //ctx->bc_v->apply_to_system(ctx->L, ctx->rho, ctx->gamma);


        PetscCall(VecSet(ctx->gamma, 0.));
        PetscCall(KSPSolve(ctx->inner_L_ksp, ctx->rho,  ctx->gamma));

        //    xt = xt + P*gamma_1
        PetscCall(MatMultAdd(ctx->P, ctx->gamma, xt, xt));


        
        PetscCall(MatMult(T, xt, ctx->tmp_1));                    // tmp_1 = [T]*xt
        PetscCall(MatMultAdd(T_c, xo, ctx->tmp_1, ctx->tmp_1));   // tmp_1 = tmp_1 + [coupling_t]*xo
        PetscCall(VecAYPX(ctx->tmp_1, -1., rt));                  // tmp_1 = rt - tmp_1


        

        //    tmp_2 = ro - [O]*xo - [coupling_o]*xt
        PetscCall(MatMult(O, xo, ctx->tmp_2));                    // tmp_2 = [O]*xo
        PetscCall(MatMultAdd(O_c, xt, ctx->tmp_2, ctx->tmp_2));   // tmp_2 = tmp_2 + [coupling_o]*xt
        PetscCall(VecAYPX(ctx->tmp_2, -1., ro));                  // tmp_2 = ro - tmp_2

        /*


        PetscCall(MatMultTranspose(ctx->G_T, ctx->tmp_1, ctx->zeta_1));
        PetscCall(MatMultTranspose(ctx->I_O, ctx->tmp_2, ctx->zeta_2));

        petsc_util::petsc_save_ascii_vec(ctx->zeta_2, "zeta_d.txt");
        petsc_util::petsc_save_ascii_vec(ctx->zeta_1, "zeta_1d.txt");

        PetscCall(VecSet(ctx->kappa_1, 0.));
        PetscCall(VecSet(ctx->kappa_2, 0.));

        ctx->bc_T__H1->apply_to_vec(ctx->zeta_1);
        ctx->bc_Oo_H1->apply_to_vec(ctx->zeta_2);
        ctx->bc_Oi_H1->apply_to_vec(ctx->zeta_2);

        petsc_util::petsc_save_ascii_vec(ctx->zeta_2, "zeta_d1.txt");
        petsc_util::petsc_save_ascii_vec(ctx->zeta_1, "zeta_1d1.txt");
        
        //ctx->bc_T__H1->apply_to_system(ctx->Q_T, ctx->zeta_1, ctx->kappa_1);
        //ctx->bc_Oo_H1->apply_to_system(ctx->Q_O, ctx->zeta_2, ctx->kappa_2);
        //ctx->bc_Oi_H1->apply_to_system(ctx->Q_O, ctx->zeta_2, ctx->kappa_2);

        PetscCall(KSPSolve(ctx->inner_Q_T_ksp, ctx->zeta_1, ctx->kappa_1));
        PetscCall(KSPSolve(ctx->inner_Q_O_ksp, ctx->zeta_2, ctx->kappa_2));

        PetscCall(MatMultAdd(ctx->G_T, ctx->kappa_1, xt, xt));
        PetscCall(MatMultAdd(ctx->I_O, ctx->kappa_2, xo, xo));

        petsc_util::petsc_save_ascii_vec(ctx->kappa_2, "kappa_d.txt");
        petsc_util::petsc_save_ascii_vec(ctx->kappa_1, "kappa_1d.txt");

        //*/

        //*

        //    zeta = G^T * tmp_1   (H1)
        PetscCall(MatMultTranspose(ctx->G, ctx->tmp_1, ctx->zeta));



        //PetscCall(VecScale(ctx->zeta, -1.));
        //ctx->bc_s1->apply_to_vec(ctx->zeta);


        //    zeta = zeta + I^T * tmp_2   (H1)
        PetscCall(MatMultTransposeAdd(ctx->I, ctx->tmp_2, ctx->zeta, ctx->zeta));
        //PetscCall(MatMultTranspose(ctx->I, ctx->tmp_2, ctx->zeta));
        
        //petsc_util::petsc_save_ascii_vec(ctx->zeta, "zeta_c.txt");

        
        //     enforce Dirichlet on inner-solve RHS and initial guess,
        //            then solve each inner system with one BoomerAMG V-cycle.
        PetscCall(VecSet(ctx->kappa, 0.));


        ctx->bc_s2->apply_to_vec(ctx->zeta); 

        //ctx->bc_s1->apply_to_vec(ctx->zeta); 
        //ctx->bc_s3->apply_to_vec(ctx->zeta); 

        //ctx->bc_s3->apply_to_vec(ctx->zeta); 
        //ctx->bc_s2->apply_to_system(ctx->Q, ctx->zeta, ctx->kappa);
        //petsc_util::petsc_save_ascii_vec(ctx->zeta, "zeta_c_1.txt");

        //petsc_util::petsc_save_ascii_mat(ctx->Q, "Q.txt");
        //petsc_util::petsc_save_ascii_mat(ctx->P, "P.txt");
        //petsc_util::petsc_save_ascii_mat(ctx->L, "L.txt");
        //petsc_util::petsc_save_ascii_mat(ctx->G, "G.txt");
        //petsc_util::petsc_save_ascii_mat(ctx->I, "I.txt");
        //petsc_util::petsc_save_ascii_mat(T, "T.txt");
        //petsc_util::petsc_save_ascii_mat(T_c, "T_c.txt");
        //petsc_util::petsc_save_ascii_mat(O, "O.txt");
        //petsc_util::petsc_save_ascii_mat(O_c, "O_c.txt");

        //petsc_util::petsc_save_ascii_vec(ro, "ro.txt");
        //petsc_util::petsc_save_ascii_vec(rt, "rt.txt");

        //petsc_util::petsc_save_ascii_mat(ctx->Q_O, "Q_O.txt");


    
        //    Single BoomerAMG V-cycle for each (no inner GMRES; KSPPREONLY).
        PetscCall(KSPSolve(ctx->inner_Q_ksp, ctx->zeta, ctx->kappa));

        //petsc_util::petsc_save_ascii_vec(ctx->kappa, "kappa_c.txt");



        //    xt = xt + P*gamma_1 + G*kappa
        //    xt = xt + P*gamma_1
        //PetscCall(MatMultAdd(ctx->P, ctx->gamma, xt, xt));
        //    xt = xt + G*kappa
        PetscCall(MatMultAdd(ctx->G, ctx->kappa, xt, xt));

        
        //ctx->bc_s3->apply_to_vec(ctx->kappa);

        //    xo = xo + I*kappa
        PetscCall(MatMultAdd(ctx->I, ctx->kappa, xo, xo));

        ctx->bc_T__Hcurl->apply_to_vec(xt); 
        ctx->bc_Oo_H1->apply_to_vec(xo); 
        ctx->bc_Oi_H1->apply_to_vec(xo); 

        //*/


    }



    la_kernel::write_to_vec(block_x, x);
    //petsc_util::petsc_save_ascii_vec(x, "x_d.txt");
    //petsc_util::petsc_save_ascii_vec(r, "r_d.txt");
    //petsc_util::petsc_save_ascii_vec(x, "x_c.txt");
    //petsc_util::petsc_save_ascii_vec(r, "r_c.txt");

    
    
    // post-smoother (continues from updated z)
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    //petsc_util::petsc_save_ascii_vec(x, "x_d_1.txt");
    //petsc_util::petsc_save_ascii_vec(x, "x_c_1.txt");


    //exit(0);

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
    Mat Q_T, Mat Q_O, Mat G_T, Mat I_O,
    Mat P, Mat G, Mat I, Mat L, Mat Q,
    //Mat X,
    Block_Rack* br_system,
    Dirichlet_BC *bc_v, Dirichlet_BC *bc_s1, Dirichlet_BC *bc_s2, Dirichlet_BC *bc_s3,
    Dirichlet_BC *bc_T__H1,                 
    Dirichlet_BC *bc_Oo_H1,                  
    Dirichlet_BC *bc_Oi_H1, 
    Dirichlet_BC *bc_T__Hcurl,
    PetscReal rtol, PetscInt max_iters, bool enable_monitor)
{
    AMS_Context *ctx;
    KSP outer;
    PC  outer_pc;
    PetscFunctionBeginUser;

    bc_v->apply_to_mat(L);
    bc_s2->apply_to_mat(Q);
    //bc_s3->apply_to_mat(Q);
    //bc_s1->apply_to_mat(Q);

    bc_T__H1->apply_to_mat(Q_T);
    bc_Oo_H1->apply_to_mat(Q_O);
    bc_Oi_H1->apply_to_mat(Q_O);
 
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

    ctx->Q_T = Q_T;
    ctx->Q_O = Q_O;
    ctx->G_T = G_T;
    ctx->I_O = I_O;
    ctx->bc_T__H1 = bc_T__H1;
    ctx->bc_Oo_H1 = bc_Oo_H1;
    ctx->bc_Oi_H1 = bc_Oi_H1;

    ctx->bc_T__Hcurl = bc_T__Hcurl;

    //ctx->X = X;
 
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

    Mat M, M1, M2, M3, M4;
    Mat tM, oM;
    const Mat O_   = ctx->br_system->get_block_lhs(0,0);
    const Mat O_c_ = ctx->br_system->get_block_lhs(0,1);
    const Mat T_   = ctx->br_system->get_block_lhs(1,1);
    const Mat T_c_ = ctx->br_system->get_block_lhs(1,0);
    PetscCall(MatPtAP(T_,G,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M1));
    PetscCall(MatPtAP(O_,I,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M2));

    //bc_s1->apply_to_mat(M1);
    
    //bc_s2->apply_to_mat(M2);
    //bc_s3->apply_to_mat(M2);
    //exit(0);

    PetscCall(MatTransposeMatMult(ctx->G,T_c_,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&tM));
    PetscCall(MatTransposeMatMult(ctx->I,O_c_,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&oM));
    PetscCall(MatMatMult(tM, ctx->I,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M3));
    PetscCall(MatMatMult(oM, ctx->G,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M4));

    //petsc_util::petsc_save_ascii_mat(M1, "M1.txt");
    //petsc_util::petsc_save_ascii_mat(M2, "M2.txt");
    //petsc_util::petsc_save_ascii_mat(M3, "M3.txt");
    //petsc_util::petsc_save_ascii_mat(M4, "M4.txt");

    //petsc_util::petsc_save_ascii_mat(Q, "Q.txt");

    //exit(0);
    //petsc_util::petsc_save_ascii_mat(Q_T, "Q_T.txt");

    PetscCall(MatDuplicate(M1, MAT_COPY_VALUES,&M));
    PetscCall(MatAXPY(M, 1.0, M2, DIFFERENT_NONZERO_PATTERN));
    PetscCall(MatAXPY(M, 1.0, M3, DIFFERENT_NONZERO_PATTERN));
    PetscCall(MatAXPY(M, 1.0, M4, DIFFERENT_NONZERO_PATTERN));
    
    //bc_s2->apply_to_mat(M);

    //bc_s1->apply_to_mat(M);
    //bc_s3->apply_to_mat(M);

    //petsc_util::petsc_save_ascii_mat(M, "M.txt");

    //exit(0);

    PetscCall(KSPSetOptionsPrefix(ctx->inner_Q_ksp, "inner_Q_"));
    PetscCall(KSPSetOperators(ctx->inner_Q_ksp, Q, Q));
    PetscCall(KSPSetOperators(ctx->inner_Q_ksp, M, M));
    //PetscCall(KSPSetFromOptions(ctx->inner_Q_ksp));
    PetscCall(KSPSetUp(ctx->inner_Q_ksp));

    // ---------- Inner KSP for Q_T ----------
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Q_T_ksp));
    PetscCall(KSPSetType(ctx->inner_Q_T_ksp, KSPCG));
    {
        PC pc_in;
        PetscCall(KSPGetPC(ctx->inner_Q_T_ksp, &pc_in));
        PetscCall(PCSetType(pc_in, PCHYPRE));
        PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
    }

    Mat Q_T_, Q_T1, Q_T2;
    Mat Q_O_, Q_O1, Q_O2;
    Mat tQ, oQ;
    
    PetscCall(MatPtAP(T_,G_T,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&Q_T1));
    PetscCall(MatPtAP(O_,I_O,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&Q_O1));

    //bc_s1->apply_to_mat(M1);
    //bc_s2->apply_to_mat(M2);
    //bc_s3->apply_to_mat(M2);

    PetscCall(MatTransposeMatMult(ctx->G_T,T_c_,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&tQ));
    PetscCall(MatTransposeMatMult(ctx->I_O,O_c_,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&oQ));
    PetscCall(MatMatMult(tQ, ctx->I_O,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&Q_T2));
    PetscCall(MatMatMult(oQ, ctx->G_T,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&Q_O2));
    PetscCall(MatDuplicate(Q_T1, MAT_COPY_VALUES,&Q_T_));
    PetscCall(MatDuplicate(Q_O1, MAT_COPY_VALUES,&Q_O_));
    //PetscCall(MatAXPY(Q_T_, 1.0, Q_T2, DIFFERENT_NONZERO_PATTERN));
    //PetscCall(MatAXPY(Q_O_, 1.0, Q_O2, DIFFERENT_NONZERO_PATTERN));
    bc_T__H1->apply_to_mat(Q_T_);
    bc_Oo_H1->apply_to_mat(Q_O_);
    bc_Oi_H1->apply_to_mat(Q_O_);


    PetscCall(KSPSetOptionsPrefix(ctx->inner_Q_T_ksp, "inner_Q_T_"));
    //PetscCall(KSPSetOperators(ctx->inner_Q_T_ksp, Q_T, Q_T));
    PetscCall(KSPSetOperators(ctx->inner_Q_T_ksp, Q_T_, Q_T_));

    PetscCall(KSPSetFromOptions(ctx->inner_Q_T_ksp));
    PetscCall(KSPSetUp(ctx->inner_Q_T_ksp));
 
    // ---------- Inner KSP for Q_O ----------
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Q_O_ksp));
    PetscCall(KSPSetType(ctx->inner_Q_O_ksp, KSPPREONLY));
    {
        PC pc_in;
        PetscCall(KSPGetPC(ctx->inner_Q_O_ksp, &pc_in));
        PetscCall(PCSetType(pc_in, PCHYPRE));
        PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
    }

    PetscCall(KSPSetOptionsPrefix(ctx->inner_Q_O_ksp, "inner_Q_O_"));
    //PetscCall(KSPSetOperators(ctx->inner_Q_O_ksp, Q_O, Q_O));
    PetscCall(KSPSetOperators(ctx->inner_Q_O_ksp, Q_O_, Q_O_));
    PetscCall(KSPSetFromOptions(ctx->inner_Q_O_ksp));
    PetscCall(KSPSetUp(ctx->inner_Q_O_ksp));

    // ---------- Inner KSP for X ----------
    /*
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_X_ksp));
    PetscCall(KSPSetType(ctx->inner_X_ksp, KSPPREONLY));
    {
        PC pc_in;
        PetscCall(KSPGetPC(ctx->inner_X_ksp, &pc_in));
        PetscCall(PCSetType(pc_in, PCHYPRE));
        PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
    }

    PetscCall(KSPSetOptionsPrefix(ctx->inner_X_ksp, "inner_X_"));
    PetscCall(KSPSetOperators(ctx->inner_X_ksp, X, X));
    PetscCall(KSPSetFromOptions(ctx->inner_X_ksp));
    PetscCall(KSPSetUp(ctx->inner_X_ksp));
    */
 
    // ---------- Initialize vectors ----------
    const Mat O   = ctx->br_system->get_block_lhs(0,0);
    const Mat T   = ctx->br_system->get_block_lhs(1,1);
    PetscCall(MatCreateVecs(T, NULL, &ctx->tmp_1));       // size #edge in T-field
    PetscCall(MatCreateVecs(O, NULL, &ctx->tmp_2));       // size #node in Omega-field
    PetscCall(MatCreateVecs(L, &ctx->gamma, &ctx->rho));  // size 3 * #node in T-field
    PetscCall(MatCreateVecs(Q, &ctx->kappa, &ctx->zeta)); // size     #node in global field

    PetscCall(MatCreateVecs(Q_T, &ctx->kappa_1, &ctx->zeta_1)); 
    PetscCall(MatCreateVecs(Q_O, &ctx->kappa_2, &ctx->zeta_2)); 
    //PetscCall(MatCreateVecs(X, NULL, &ctx->x)); 

    //petsc_util::petsc_save_ascii_vec(ctx->zeta, "zeta_initial.txt");
 
    // ---------- Outer KSP + PCShell(AMS) ----------
    const Mat A = ctx->br_system->get_lhs();
    const Vec r = ctx->br_system->get_rhs();
    const Vec x = ctx->br_system->get_x();

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &outer));
    //PetscCall(KSPSetType(outer, KSPGMRES));   
    //PetscCall(KSPGMRESSetRestart(outer, 1000));  

    KSPSetType(outer, KSPCG);
    //KSPSetType(outer, KSPCGS);

    //PetscCall(KSPSetType(outer, KSPFGMRES));
    //PetscCall(KSPSetPCSide(outer, PC_RIGHT));
    //PetscCall(KSPSetNormType(outer, KSP_NORM_UNPRECONDITIONED));
    //PetscCall(KSPGMRESSetRestart(outer, 200));

    PetscCall(KSPSetOperators(outer, A, A));
    PetscCall(KSPSetTolerances(outer, rtol, PETSC_DEFAULT, PETSC_DEFAULT, max_iters));

    PetscCall(KSPSetNormType(outer, KSP_NORM_UNPRECONDITIONED));  // true residual
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


