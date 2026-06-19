#include "physics/electromagnetism/eddy_current/auxiliary_solver.h"


using namespace simu;



PetscErrorCode T_Omega_AMS::AMS_apply_global(PC pc, Vec r, Vec x)
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
    
    for(int n=0; n<2; ++n){
        //    tmp_1 = rt - [T]*xt - [coupling_t]*xo
        PetscCall(MatMult(T, xt, ctx->tmp_1));                    // tmp_1 = [T]*xt
        PetscCall(MatMultAdd(T_c, xo, ctx->tmp_1, ctx->tmp_1));   // tmp_1 = tmp_1 + [coupling_t]*xo
        PetscCall(VecAYPX(ctx->tmp_1, -1., rt));                  // tmp_1 = rt - tmp_1

        //    rho  = P^T * tmp_1   ((H1)^3)
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_1, ctx->rho));

        ctx->bc_v->apply_to_vec(ctx->rho);  

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

        //*

        //    zeta = G^T * tmp_1   (H1)
        PetscCall(MatMultTranspose(ctx->G, ctx->tmp_1, ctx->zeta));

        //    zeta = zeta + I^T * tmp_2   (H1)
        PetscCall(MatMultTransposeAdd(ctx->I, ctx->tmp_2, ctx->zeta, ctx->zeta));
        
        //     enforce Dirichlet on inner-solve RHS and initial guess,
        //            then solve each inner system with one BoomerAMG V-cycle.
        PetscCall(VecSet(ctx->kappa, 0.));

        ctx->bc_s2->apply_to_vec(ctx->zeta); 
    
        //    Single BoomerAMG V-cycle for each (no inner GMRES; KSPPREONLY).
        PetscCall(KSPSolve(ctx->inner_Q_ksp, ctx->zeta, ctx->kappa));

        //    xt = xt + P*gamma_1 + G*kappa
        //    xt = xt + P*gamma_1
        //PetscCall(MatMultAdd(ctx->P, ctx->gamma, xt, xt));
        //    xt = xt + G*kappa
        PetscCall(MatMultAdd(ctx->G, ctx->kappa, xt, xt));

        //    xo = xo + I*kappa
        PetscCall(MatMultAdd(ctx->I, ctx->kappa, xo, xo));
    }
    la_kernel::write_to_vec(block_x, x);

    // post-smoother
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    PetscFunctionReturn(0);
}




PetscErrorCode T_Omega_AMS::AMS_apply_decoupled(PC pc, Vec r, Vec x)
{

    AMS_Context *ctx;
    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));

    const Mat A = ctx->br_system->get_lhs();

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
    
    for(int n=0; n<2; ++n){
        //    tmp_1 = rt - [T]*xt - [coupling_t]*xo
        PetscCall(MatMult(T, xt, ctx->tmp_1));                    // tmp_1 = [T]*xt
        PetscCall(MatMultAdd(T_c, xo, ctx->tmp_1, ctx->tmp_1));   // tmp_1 = tmp_1 + [coupling_t]*xo
        PetscCall(VecAYPX(ctx->tmp_1, -1., rt));                  // tmp_1 = rt - tmp_1
        
        //    rho  = P^T * tmp_1   ((H1)^3)
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_1, ctx->rho));

        ctx->bc_v->apply_to_vec(ctx->rho);  

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


        PetscCall(MatMultTranspose(ctx->G_T, ctx->tmp_1, ctx->zeta_1));
        PetscCall(MatMultTranspose(ctx->I_O, ctx->tmp_2, ctx->zeta_2));

        PetscCall(VecSet(ctx->kappa_1, 0.));
        PetscCall(VecSet(ctx->kappa_2, 0.));


        ctx->bc_T__H1->apply_to_vec(ctx->zeta_1);
        ctx->bc_Oo_H1->apply_to_vec(ctx->zeta_2);
        ctx->bc_Oi_H1->apply_to_vec(ctx->zeta_2);

        
        PetscCall(KSPSolve(ctx->inner_Q_T_ksp, ctx->zeta_1, ctx->kappa_1));
        PetscCall(KSPSolve(ctx->inner_Q_O_ksp, ctx->zeta_2, ctx->kappa_2));

        PetscCall(MatMultAdd(ctx->G_T, ctx->kappa_1, xt, xt));
        PetscCall(MatMultAdd(ctx->I_O, ctx->kappa_2, xo, xo));
    }

    la_kernel::write_to_vec(block_x, x);
    
    // post-smoother (continues from updated z)
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    PetscFunctionReturn(0);
}





PetscErrorCode T_Omega_AMS::AMS_apply_coupled(PC pc, Vec r, Vec x)
{

    AMS_Context *ctx;
    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));

    const Mat A = ctx->br_system->get_lhs();

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

    Vec zeta_l[2] = {ctx->zeta_1, ctx->zeta_2};
    Vec kappa_l[2] = {ctx->kappa_1, ctx->kappa_2};

    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, zeta_l, &ctx->zeta));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, kappa_l, &ctx->kappa));

    
    for(int n=0; n<30; ++n){
        //    tmp_1 = rt - [T]*xt - [coupling_t]*xo
        PetscCall(MatMult(T, xt, ctx->tmp_1));                    // tmp_1 = [T]*xt
        PetscCall(MatMultAdd(T_c, xo, ctx->tmp_1, ctx->tmp_1));   // tmp_1 = tmp_1 + [coupling_t]*xo
        PetscCall(VecAYPX(ctx->tmp_1, -1., rt));                  // tmp_1 = rt - tmp_1
        
        //    rho  = P^T * tmp_1   ((H1)^3)
        PetscCall(MatMultTranspose(ctx->P, ctx->tmp_1, ctx->rho));

        ctx->bc_v->apply_to_vec(ctx->rho);  

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



        PetscCall(MatMultTranspose(ctx->G_T, ctx->tmp_1, ctx->zeta_1));
        PetscCall(MatMultTranspose(ctx->I_O, ctx->tmp_2, ctx->zeta_2));


        PetscCall(VecZeroEntries(ctx->kappa)); 


        //ctx->bc_T__H1->apply_to_vec(ctx->zeta_1);
        ctx->bc_Oo_H1->apply_to_vec(ctx->zeta_2);
        //ctx->bc_Oi_H1->apply_to_vec(ctx->zeta_2);

        
        PetscCall(KSPSolve(ctx->inner_Q_ksp, ctx->zeta, ctx->kappa));

        PetscCall(MatMultAdd(ctx->G_T, ctx->kappa_1, xt, xt));
        PetscCall(MatMultAdd(ctx->I_O, ctx->kappa_2, xo, xo));
    }

    la_kernel::write_to_vec(block_x, x);
    
    // post-smoother (continues from updated z)
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    PetscFunctionReturn(0);
}






PetscErrorCode T_Omega_AMS::AMS_apply_fully_coupled(PC pc, Vec r, Vec x)
{

    AMS_Context *ctx;
    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));

    const Mat A = ctx->br_system->get_lhs();

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

    Vec rho_l[2] = {ctx->rho_1, ctx->rho_2};
    Vec gamma_l[2] = {ctx->gamma_1, ctx->gamma_2};

    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, rho_l, &ctx->rho));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, gamma_l, &ctx->gamma));

    Vec zeta_l[2] = {ctx->zeta_1, ctx->zeta_2};
    Vec kappa_l[2] = {ctx->kappa_1, ctx->kappa_2};

    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, zeta_l, &ctx->zeta));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, kappa_l, &ctx->kappa));

    
    for(int n=0; n<5; ++n){
        //    tmp_1 = rt - [T]*xt - [coupling_t]*xo
        PetscCall(MatMult(T, xt, ctx->tmp_1));                    // tmp_1 = [T]*xt
        PetscCall(MatMultAdd(T_c, xo, ctx->tmp_1, ctx->tmp_1));   // tmp_1 = tmp_1 + [coupling_t]*xo
        PetscCall(VecAYPX(ctx->tmp_1, -1., rt));                  // tmp_1 = rt - tmp_1

        //    tmp_2 = ro - [O]*xo - [coupling_o]*xt
        PetscCall(MatMult(O, xo, ctx->tmp_2));                    // tmp_2 = [O]*xo
        PetscCall(MatMultAdd(O_c, xt, ctx->tmp_2, ctx->tmp_2));   // tmp_2 = tmp_2 + [coupling_o]*xt
        PetscCall(VecAYPX(ctx->tmp_2, -1., ro));                  // tmp_2 = ro - tmp_2
        
        //    rho_1  = P^T * tmp_1   ((H1)^3)
        PetscCall(MatMultTranspose(ctx->P,   ctx->tmp_1, ctx->rho_1));
        //    rho_2  = I^T * tmp_2     H1
        PetscCall(MatMultTranspose(ctx->I_O, ctx->tmp_2, ctx->rho_2));

        //ctx->bc_v->apply_to_vec(ctx->rho);  
        ctx->bc_Oo_H1->apply_to_vec(ctx->rho_2);

        PetscCall(VecZeroEntries(ctx->gamma)); 

        PetscCall(KSPSolve(ctx->inner_L_ksp, ctx->rho,  ctx->gamma));

        //    xt = xt + P*gamma_1
        PetscCall(MatMultAdd(ctx->P, ctx->gamma_1, xt, xt));
        //    xo = xo + I*gamma_2
        PetscCall(MatMultAdd(ctx->I_O, ctx->gamma_2, xo, xo));

        PetscCall(MatMult(T, xt, ctx->tmp_1));                    // tmp_1 = [T]*xt
        PetscCall(MatMultAdd(T_c, xo, ctx->tmp_1, ctx->tmp_1));   // tmp_1 = tmp_1 + [coupling_t]*xo
        PetscCall(VecAYPX(ctx->tmp_1, -1., rt));                  // tmp_1 = rt - tmp_1


        //    tmp_2 = ro - [O]*xo - [coupling_o]*xt
        PetscCall(MatMult(O, xo, ctx->tmp_2));                    // tmp_2 = [O]*xo
        PetscCall(MatMultAdd(O_c, xt, ctx->tmp_2, ctx->tmp_2));   // tmp_2 = tmp_2 + [coupling_o]*xt
        PetscCall(VecAYPX(ctx->tmp_2, -1., ro));                  // tmp_2 = ro - tmp_2



        PetscCall(MatMultTranspose(ctx->G_T, ctx->tmp_1, ctx->zeta_1));
        PetscCall(MatMultTranspose(ctx->I_O, ctx->tmp_2, ctx->zeta_2));


        PetscCall(VecZeroEntries(ctx->kappa)); 


        //ctx->bc_T__H1->apply_to_vec(ctx->zeta_1);
        ctx->bc_Oo_H1->apply_to_vec(ctx->zeta_2);
        //ctx->bc_Oi_H1->apply_to_vec(ctx->zeta_2);

        
        PetscCall(KSPSolve(ctx->inner_Q_ksp, ctx->zeta, ctx->kappa));

        PetscCall(MatMultAdd(ctx->G_T, ctx->kappa_1, xt, xt));
        PetscCall(MatMultAdd(ctx->I_O, ctx->kappa_2, xo, xo));
    }

    la_kernel::write_to_vec(block_x, x);
    
    // post-smoother (continues from updated z)
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

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
    PetscCall(VecDestroy(&ctx->zeta_1));
    PetscCall(VecDestroy(&ctx->zeta_2));
    PetscCall(VecDestroy(&ctx->kappa_1));
    PetscCall(VecDestroy(&ctx->kappa_2));
    PetscCall(PetscFree(ctx));
    PetscFunctionReturn(0);
}




PetscErrorCode T_Omega_AMS::solve_AMS(
    AMS_Info& AMS_info,
    int algorithm_id,
    Mat Q_T, Mat Q_O, Mat G_T, Mat I_O,
    Mat P, Mat G, Mat I, Mat L, Mat Q,
    //Mat X,
    Block_Rack* br_system,
    Dirichlet_BC *bc_v, Dirichlet_BC *bc_s1, Dirichlet_BC *bc_s2, Dirichlet_BC *bc_s3,
    Dirichlet_BC *bc_T__H1,                 
    Dirichlet_BC *bc_Oo_H1,                  
    Dirichlet_BC *bc_Oi_H1, 
    Dirichlet_BC *bc_T__Hcurl,
    Dirichlet_BC *bc_Ti_H1, 
    Dirichlet_BC *bc_Ot_H1, 
    PetscReal rtol, PetscInt max_iters, bool enable_monitor)
{
    AMS_Context *ctx;
    KSP outer;
    PC  outer_pc;
    PetscFunctionBeginUser;

    bc_v->apply_to_mat(L);

    //bc_s2->apply_to_mat(Q);
    //bc_s3->apply_to_mat(Q);
    //bc_s1->apply_to_mat(Q);

    //bc_T__H1->apply_to_mat(Q_T);
    //bc_Oo_H1->apply_to_mat(Q_O);
    //bc_Oi_H1->apply_to_mat(Q_O);
 
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

    ctx->bc_Ti_H1 = bc_Ti_H1;
    ctx->bc_Ot_H1 = bc_Ot_H1;

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
    

    if(algorithm_id == 2 || algorithm_id == 3 || algorithm_id == 4)
    {
        // ---------- Inner KSP for global Q ----------
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
        //PetscCall(KSPSetOperators(ctx->inner_Q_ksp, M, M));
        //PetscCall(KSPSetFromOptions(ctx->inner_Q_ksp));
        PetscCall(KSPSetUp(ctx->inner_Q_ksp));

        // ---------- Initialize vectors ----------
        if(algorithm_id == 2){
            PetscCall(MatCreateVecs(L, &ctx->gamma, &ctx->rho));  // size 3 * #node in T-field
            PetscCall(MatCreateVecs(Q, &ctx->kappa, &ctx->zeta)); // size     #node in global field
        }

        if(algorithm_id == 3){
            PetscCall(MatCreateVecs(L, &ctx->gamma, &ctx->rho));  // size 3 * #node in T-field
        }


        if(algorithm_id ==4){
            PetscCall(MatCreateVecs(P, &ctx->rho_1, NULL)); 
            PetscCall(MatCreateVecs(P, &ctx->gamma_1, NULL)); 
            PetscCall(MatCreateVecs(Q_O, &ctx->gamma_2, &ctx->rho_2)); 
        }
        
        
        PetscCall(MatCreateVecs(Q_T, &ctx->kappa_1, &ctx->zeta_1)); 
        PetscCall(MatCreateVecs(Q_O, &ctx->kappa_2, &ctx->zeta_2)); 
    }

    if(algorithm_id == 1)
    {
        // ---------- Inner KSP for Q_T ----------
        PetscCall(KSPCreate(PETSC_COMM_WORLD, &ctx->inner_Q_T_ksp));
        PetscCall(KSPSetType(ctx->inner_Q_T_ksp, KSPCG));
        {
            PC pc_in;
            PetscCall(KSPGetPC(ctx->inner_Q_T_ksp, &pc_in));
            PetscCall(PCSetType(pc_in, PCHYPRE));
            PetscCall(PCHYPRESetType(pc_in, "boomeramg"));
        }

        /*

        Mat Q_T_, Q_T1, Q_T2;
        Mat Q_O_, Q_O1, Q_O2;
        Mat tQ, oQ;
        
        PetscCall(MatPtAP(T_,G_T,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&Q_T1));
        PetscCall(MatPtAP(O_,I_O,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&Q_O1));



        PetscCall(MatTransposeMatMult(ctx->G_T,T_c_,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&tQ));
        PetscCall(MatTransposeMatMult(ctx->I_O,O_c_,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&oQ));
        PetscCall(MatMatMult(tQ, ctx->I_O,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&Q_T2));
        PetscCall(MatMatMult(oQ, ctx->G_T,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&Q_O2));
        PetscCall(MatDuplicate(Q_T1, MAT_COPY_VALUES,&Q_T_));
        PetscCall(MatDuplicate(Q_O1, MAT_COPY_VALUES,&Q_O_));

        bc_T__H1->apply_to_mat(Q_T_);
        bc_Oo_H1->apply_to_mat(Q_O_);
        bc_Oi_H1->apply_to_mat(Q_O_);

        */


        PetscCall(KSPSetOptionsPrefix(ctx->inner_Q_T_ksp, "inner_Q_T_"));
        PetscCall(KSPSetOperators(ctx->inner_Q_T_ksp, Q_T, Q_T));
        //PetscCall(KSPSetOperators(ctx->inner_Q_T_ksp, Q_T_, Q_T_));

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
        PetscCall(KSPSetOperators(ctx->inner_Q_O_ksp, Q_O, Q_O));
        //PetscCall(KSPSetOperators(ctx->inner_Q_O_ksp, Q_O_, Q_O_));
        PetscCall(KSPSetFromOptions(ctx->inner_Q_O_ksp));
        PetscCall(KSPSetUp(ctx->inner_Q_O_ksp));

        // ---------- Initialize vectors ----------
        PetscCall(MatCreateVecs(L, &ctx->gamma, &ctx->rho));  // size 3 * #node in T-field
        PetscCall(MatCreateVecs(Q_T, &ctx->kappa_1, &ctx->zeta_1)); 
        PetscCall(MatCreateVecs(Q_O, &ctx->kappa_2, &ctx->zeta_2)); 
    }

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
    
    //PetscCall(MatCreateVecs(M, &ctx->kappa, &ctx->zeta)); // size     #node in global field
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
    switch (algorithm_id)
    {
    case 1: PetscCall(PCShellSetApply  (outer_pc, AMS_apply_decoupled)); break;
    case 2: PetscCall(PCShellSetApply  (outer_pc, AMS_apply_global)); break;
    case 3: PetscCall(PCShellSetApply  (outer_pc, AMS_apply_coupled)); break;
    case 4: PetscCall(PCShellSetApply  (outer_pc, AMS_apply_fully_coupled)); break;
    default:
        Logger::error("T_Omega_AMS::solve_AMS: unknown AMS algorithm id: "+std::to_string(algorithm_id));
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


