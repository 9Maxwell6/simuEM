#include "physics/electromagnetism/eddy_current/auxiliary_solver_complex.h"


using namespace simu;


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
    const Mat Mc  = ctx->br_system->get_block_lhs(0,0);
    const Mat X_  = ctx->br_system->get_block_lhs(0,1);
    const Mat Xt  = ctx->br_system->get_block_lhs(1,0);
    const Mat Ki  = ctx->br_system->get_block_lhs(1,1);
    const Mat Kc_r  = ctx->br_system->get_block_lhs(0,2);
    const Mat Kc_c  = ctx->br_system->get_block_lhs(2,0);


 
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

    std::vector<Vec> rho_c_l = {ctx->rho_1, ctx->rho_3};
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, nullptr, rho_c_l.data(), &rho_c));

    std::vector<Vec> gamma_c_l = {ctx->gamma_1, ctx->gamma_3};
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, nullptr, gamma_c_l.data(), &gamma_c));

    
    for(int n=0; n<1; ++n){
        //  rho_1 = Pt * (r1 - Mc*x1 - X*x2 - Kc*x3)
        PetscCall(MatMult(Mc, x1, ctx->tmp_1));                   
        PetscCall(MatMultAdd(X_, x2, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(MatMultAdd(Kc_r, x3, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(VecAYPX(ctx->tmp_1, -1., r1));                  
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_1, ctx->rho_1));

        //  rho_3 = Pt * (r3 - Kc*x1 - Mc*x3 - X*x4)
        PetscCall(MatMult(Kc_c, x1, ctx->tmp_3));                   
        PetscCall(MatMultAdd(Mc, x3, ctx->tmp_3, ctx->tmp_3));   
        PetscCall(MatMultAdd(X_, x4, ctx->tmp_3, ctx->tmp_3));   
        PetscCall(VecAYPX(ctx->tmp_3, -1., r3));                  
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_3, ctx->rho_3));

        ctx->bc_T_1_H1_v->apply_to_vec(ctx->rho_1);  
        ctx->bc_T_1_H1_v->apply_to_vec(ctx->rho_3);  

        PetscCall(VecSet(ctx->gamma_1, 0.));
        PetscCall(VecSet(ctx->gamma_3, 0.));
        //PetscCall(KSPSolve(ctx->inner_LpM_ksp, ctx->rho_1,  ctx->gamma_1));
        //PetscCall(KSPSolve(ctx->inner_LpM_ksp, ctx->rho_3,  ctx->gamma_3));
        PetscCall(KSPSolve(ctx->inner_LpM_ksp, rho_c, gamma_c));

        KSPConvergedReason reason;
        PetscCall(KSPGetConvergedReason(ctx->inner_LpM_ksp, &reason));
        PetscCheck(reason > 0, PETSC_COMM_SELF, PETSC_ERR_CONV_FAILED,
                "inner_LpM_ksp failed, reason %d", (int)reason);

        //    xt = xt + P*gamma_1
        PetscCall(MatMultAdd(ctx->P, ctx->gamma_1, x1, x1));
        PetscCall(MatMultAdd(ctx->P, ctx->gamma_3, x3, x3));

        //  zeta_1 = Gt * (r1 - Mc*x1 - X*x2 - Kc*x3)
        PetscCall(MatMult(Mc, x1, ctx->tmp_1));                   
        PetscCall(MatMultAdd(X_, x2, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(MatMultAdd(Kc_r, x3, ctx->tmp_1, ctx->tmp_1));   
        PetscCall(VecAYPX(ctx->tmp_1, -1., r1));                  
        PetscCall(MatMultTranspose(ctx->G, ctx->tmp_1, ctx->zeta_1));

        //  zeta_3 = Gt * (r3 - Kc*x1 - Mc*x3 - X*x4)
        PetscCall(MatMult(Kc_c, x1, ctx->tmp_3));                   
        PetscCall(MatMultAdd(Mc, x3, ctx->tmp_3, ctx->tmp_3));   
        PetscCall(MatMultAdd(X_, x4, ctx->tmp_3, ctx->tmp_3));   
        PetscCall(VecAYPX(ctx->tmp_3, -1., r3));                  
        PetscCall(MatMultTranspose(ctx->G, ctx->tmp_3, ctx->zeta_3));

        //  zeta_2 = It * (r2 - Xt*x1 - Ki*x2)
        PetscCall(MatMult(Xt, x1, ctx->tmp_2));                   
        PetscCall(MatMultAdd(Ki, x2, ctx->tmp_2, ctx->tmp_2));   
        PetscCall(VecAYPX(ctx->tmp_2, -1., r2));                  
        PetscCall(MatMultTranspose(ctx->I, ctx->tmp_2, ctx->zeta_2));

        //  zeta_4 = It * (r4 - Xt*x3 - Ki*x4)
        PetscCall(MatMult(Xt, x3, ctx->tmp_4));                   
        PetscCall(MatMultAdd(Ki, x4, ctx->tmp_4, ctx->tmp_4));   
        PetscCall(VecAYPX(ctx->tmp_4, -1., r4));                  
        PetscCall(MatMultTranspose(ctx->I, ctx->tmp_4, ctx->zeta_4));


        PetscCall(VecSet(ctx->kappa_1, 0.));
        PetscCall(VecSet(ctx->kappa_2, 0.));
        PetscCall(VecSet(ctx->kappa_3, 0.));
        PetscCall(VecSet(ctx->kappa_4, 0.));


        ctx->bc_T_1_H1_s->apply_to_vec(ctx->zeta_1);
        ctx->bc_T_1_H1_s->apply_to_vec(ctx->zeta_3);

        ctx->bc_O_o->apply_to_vec(ctx->zeta_2);
        ctx->bc_O_i->apply_to_vec(ctx->zeta_2);
        ctx->bc_O_o->apply_to_vec(ctx->zeta_4);
        ctx->bc_O_i->apply_to_vec(ctx->zeta_4);
        
        PetscCall(KSPSolve(ctx->inner_Qc_ksp, ctx->zeta_1, ctx->kappa_1));
        PetscCall(KSPSolve(ctx->inner_Qc_ksp, ctx->zeta_3, ctx->kappa_3));

        PetscCall(KSPSolve(ctx->inner_Qi_ksp, ctx->zeta_2, ctx->kappa_2));
        PetscCall(KSPSolve(ctx->inner_Qi_ksp, ctx->zeta_4, ctx->kappa_4));

        PetscCall(MatMultAdd(ctx->G, ctx->kappa_1, x1, x1));
        PetscCall(MatMultAdd(ctx->I, ctx->kappa_2, x2, x2));
        PetscCall(MatMultAdd(ctx->G, ctx->kappa_3, x3, x3));
        PetscCall(MatMultAdd(ctx->I, ctx->kappa_4, x4, x4));
    }

    la_kernel::write_to_vec(ctx->block_x, x);
    // post-smoother (continues from updated z)
    //PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    for (auto& b : ctx->block_x) PetscCall(VecDestroy(&b));
    for (auto& b : ctx->block_r) PetscCall(VecDestroy(&b));
    PetscCall(VecDestroy(&rho_c));
    PetscCall(VecDestroy(&gamma_c));
    //*/
    

    PetscFunctionReturn(0);
}

PetscErrorCode T_Omega_AMS_c::AMS_apply_coupled(PC pc, Vec r, Vec x)
{
    PetscFunctionBeginUser;

    PetscFunctionReturn(0);
}

PetscErrorCode T_Omega_AMS_c::AMS_apply_global(PC pc, Vec r, Vec x)
{

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
    PetscCall(KSPDestroy(&ctx->inner_Q_ksp));
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
    Mat Qc, Mat Qi, Mat H, Mat J,
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
    ctx->H = H;
    ctx->J = J;
    ctx->br_system = br_system;
    ctx->bc_T_1_H1_s = bc_T_1_H1_s;
    ctx->bc_T_1_H1_v = bc_T_1_H1_v;
    ctx->bc_O_o = bc_O_o;
    ctx->bc_O_i = bc_O_i;




 
    // ---------- Inner KSP for L: PREONLY + BoomerAMG (single V-cycle) ----------
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_LpM_ksp));
    PetscCall(KSPSetType(ctx->inner_LpM_ksp, KSPPREONLY));  
    {
        PC pc_in;
        PetscCall(KSPGetPC(ctx->inner_LpM_ksp, &pc_in));
        PetscCall(PCSetType(pc_in, PCHYPRE));
        PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
    }
    PetscCall(KSPSetOptionsPrefix(ctx->inner_LpM_ksp, "inner_LpM_"));
    PetscCall(KSPSetOperators(ctx->inner_LpM_ksp, LpM, LpM));
    PetscCall(KSPSetFromOptions(ctx->inner_LpM_ksp));
    PetscCall(KSPSetUp(ctx->inner_LpM_ksp)); 
    

    PetscInt N_Vcycles = 5;  // number of V-cycles for Q


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
        PetscCall(KSPSetOperators(ctx->inner_Qi_ksp, Qi, Qi));
        PetscCall(KSPSetFromOptions(ctx->inner_Qi_ksp));
        PetscCall(KSPSetUp(ctx->inner_Qi_ksp));

        // ---------- Initialize vectors ----------
        PetscCall(MatCreateVecs(P, &ctx->gamma_1, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_1, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->gamma_3, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_3, NULL));  // size 3 * #node in T-field

        PetscCall(MatCreateVecs(Qc, &ctx->kappa_1, &ctx->zeta_1)); 
        PetscCall(MatCreateVecs(Qi, &ctx->kappa_2, &ctx->zeta_2)); 
        PetscCall(MatCreateVecs(Qc, &ctx->kappa_3, &ctx->zeta_3)); 
        PetscCall(MatCreateVecs(Qi, &ctx->kappa_4, &ctx->zeta_4)); 

    }if(algorithm_id == 2){ // global
        // ---------- Inner KSP for [H] ----------
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Q_ksp));
        PetscCall(KSPSetType(ctx->inner_Q_ksp, KSPPREONLY));
        {
            PC pc_in;
            PetscCall(KSPGetPC(ctx->inner_Q_ksp, &pc_in));
            PetscCall(PCSetType(pc_in, PCHYPRE));
            PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
        }
        PetscCall(KSPSetOptionsPrefix(ctx->inner_Q_ksp, "inner_Q_"));

        // fixed number of V-cycles, done INSIDE hypre:
        char buf[32];
        PetscCall(PetscSNPrintf(buf, sizeof(buf), "%d", (int)N_Vcycles));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q_pc_hypre_boomeramg_max_iter", buf));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q_pc_hypre_boomeramg_tol", "0.0")); // 0 = never stop early
        PetscCall(KSPSetOperators(ctx->inner_Q_ksp, H, H));
        PetscCall(KSPSetFromOptions(ctx->inner_Q_ksp));
        PetscCall(KSPSetUp(ctx->inner_Q_ksp));

        // ---------- Initialize vectors ----------
        PetscCall(MatCreateVecs(P, &ctx->gamma_1, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_1, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->gamma_3, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_3, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(H, &ctx->kappa, &ctx->zeta)); // size     #node in global field

    }else if(algorithm_id == 3){
        // ---------- Inner KSP for [J] ----------
        // ---------- Inner KSP for [H] ----------
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Q_ksp));
        PetscCall(KSPSetType(ctx->inner_Q_ksp, KSPPREONLY));
        {
            PC pc_in;
            PetscCall(KSPGetPC(ctx->inner_Q_ksp, &pc_in));
            PetscCall(PCSetType(pc_in, PCHYPRE));
            PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
        }
        PetscCall(KSPSetOptionsPrefix(ctx->inner_Q_ksp, "inner_Q_"));

        // fixed number of V-cycles, done INSIDE hypre:
        char buf[32];
        PetscCall(PetscSNPrintf(buf, sizeof(buf), "%d", (int)N_Vcycles));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q_pc_hypre_boomeramg_max_iter", buf));
        PetscCall(PetscOptionsSetValue(NULL, "-inner_Q_pc_hypre_boomeramg_tol", "0.0")); // 0 = never stop early
        PetscCall(KSPSetOperators(ctx->inner_Q_ksp, J, J));
        PetscCall(KSPSetFromOptions(ctx->inner_Q_ksp));
        PetscCall(KSPSetUp(ctx->inner_Q_ksp));

        // ---------- Initialize vectors ----------
        PetscCall(MatCreateVecs(P, &ctx->gamma_1, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_1, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->gamma_3, NULL));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(P, &ctx->rho_3, NULL));  // size 3 * #node in T-field

        PetscCall(MatCreateVecs(G, &ctx->kappa_1, NULL)); 
        PetscCall(MatCreateVecs(G, &ctx->zeta_1, NULL)); 

        PetscCall(MatCreateVecs(I, &ctx->kappa_2, NULL)); 
        PetscCall(MatCreateVecs(I, &ctx->zeta_2, NULL)); 


        
    }

 
    // ---------- Initialize vectors ----------
    const Mat T   = ctx->br_system->get_block_lhs(0,0);
    const Mat O   = ctx->br_system->get_block_lhs(1,1);
    PetscCall(MatCreateVecs(T, NULL, &ctx->tmp_1));       // size #edge in T-field
    PetscCall(MatCreateVecs(O, NULL, &ctx->tmp_2));       // size #node in Omega-field
    PetscCall(MatCreateVecs(T, NULL, &ctx->tmp_3));       // size #edge in T-field
    PetscCall(MatCreateVecs(O, NULL, &ctx->tmp_4));       // size #node in Omega-field

 
    // ---------- Outer KSP + PCShell(AMS) ----------
    const Mat A = ctx->br_system->get_lhs();
    const Vec r = ctx->br_system->get_rhs();
    const Vec x = ctx->br_system->get_x();

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &outer));


    PetscCall(KSPSetType(outer, KSPGMRES));
    //PetscCall(KSPGMRESSetRestart(outer, 1000));  


    PetscCall(KSPSetOperators(outer, A, A));
    PetscCall(KSPSetTolerances(outer, rtol, PETSC_DEFAULT, PETSC_DEFAULT, max_iters));

    PetscCall(KSPSetNormType(outer, KSP_NORM_UNPRECONDITIONED));  // true residual

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


