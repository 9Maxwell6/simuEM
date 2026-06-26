#include "physics/electromagnetism/eddy_current/auxiliary_solver_complex.h"


using namespace simu;


PetscErrorCode T_Omega_AMS_c::AMS_apply_decoupled(PC pc, Vec r, Vec x)
{

    PetscFunctionReturn(0);
}

PetscErrorCode T_Omega_AMS_c::AMS_apply_coupled(PC pc, Vec r, Vec x)
{

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
    PetscReal rtol, PetscInt max_iters, bool enable_monitor)
{
    AMS_Context *ctx;
    KSP outer;
    PC  outer_pc;
    PetscFunctionBeginUser;

    bc_T_1_H1_v->apply_to_mat(LpM);

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
        PetscCall(MatCreateVecs(LpM, &ctx->gamma, &ctx->rho));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(Qc, &ctx->kappa_1, &ctx->zeta_1)); 
        PetscCall(MatCreateVecs(Qi, &ctx->kappa_2, &ctx->zeta_2)); 

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
        PetscCall(MatCreateVecs(LpM, &ctx->gamma, &ctx->rho));  // size 3 * #node in T-field
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
        PetscCall(MatCreateVecs(LpM, &ctx->gamma, &ctx->rho));  // size 3 * #node in T-field

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

 
    // ---------- Outer KSP + PCShell(AMS) ----------
    const Mat A = ctx->br_system->get_lhs();
    const Vec r = ctx->br_system->get_rhs();
    const Vec x = ctx->br_system->get_x();

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &outer));
    //PetscCall(KSPSetType(outer, KSPGMRES));   
    //PetscCall(KSPGMRESSetRestart(outer, 1000));  

    KSPSetType(outer, KSPGMRES);
    //KSPSetType(outer, KSPCGS);

    //PetscCall(KSPSetType(outer, KSPFGMRES));
    //PetscCall(KSPSetPCSide(outer, PC_RIGHT));
    //PetscCall(KSPSetNormType(outer, KSP_NORM_UNPRECONDITIONED));
    //PetscCall(KSPGMRESSetRestart(outer, 200));

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


