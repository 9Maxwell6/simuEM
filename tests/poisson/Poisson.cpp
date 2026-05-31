#include "Poisson.h"

//#include "petsc.h"
#include "gmsh.h"
//#include "mfem.hpp"


#include "math/fem/fem_space.h"
#include "math/fem/space_H1.h"
#include "math/fem/space_Hcurl.h"

#include "utils/logger.h"
#include "utils/util_la.h"

#include <Eigen/Dense>
#include <stdio.h>
//#include <config.h>
#include "world/mesh/mesh_parser.h"

// for test
#include <fstream>

using namespace simu;

Poisson::Poisson(Mesh& mesh) : mesh_(mesh), fe_system_(mesh)
{

    dim_ = mesh.get_mesh_dimension();

    space_ = H1_Space(mesh.get_mesh_dimension(), 1);

    H1_s_ = H1_Space(mesh.get_mesh_dimension(), 1);

    key_true_boundary_ = mesh_.get_keys_true_boundary()[0];  // there must be only one true boundary.
    std::string true_boundary_description = mesh.get_group_description(key_true_boundary_);
    Logger::info("Poisson - Found simulation boundary: " + true_boundary_description);
    
    
    
    for(const Key& mesh_key : mesh_.get_keys_domain())
    {
        std::string description = mesh.get_group_description(mesh_key);

        key_interior_ = mesh_key;
        mesh_.set_group_property_id(mesh_key, Domain::INTERIOR);
        break;
    }

    
    Logger::info("[Poisson] - Assign space H1 to the interior.");
    dof_field_ = fe_system_.register_FE_space(space_, key_interior_);
    
    Logger::info("[Poisson] - register Dirichlet BC to the true boundary.");
    bc_ = fe_system_.register_Dirichlet_BC(dof_field_, key_true_boundary_, Dirichlet_Type::HOMOGENEOUS);
    
    Logger::block_info(dof_field_.id, dof_field_.row_offset, dof_field_.col_offset, dof_field_.row_size, dof_field_.col_size);


    Logger::info("[Poisson preconditioner] - initialize preconditioner in H1 global field. ");
    pc_Q_ = fe_system_.register_dual_FE_space(H1_s_, H1_s_, key_interior_, nullptr, nullptr);
    Logger::block_info(pc_Q_.id, 
                       pc_Q_.row_offset, 
                       pc_Q_.col_offset, 
                       pc_Q_.row_size, 
                       pc_Q_.col_size);


    Logger::info("[Poisson preconditioner] - initialize dof mapping from global H1 to omega H1 field.");
    pc_I_ = fe_system_.register_FE_space_coupling(dof_field_, pc_Q_, key_interior_);
    Logger::block_info(pc_I_.id, 
                       pc_I_.row_offset, 
                       pc_I_.col_offset, 
                       pc_I_.row_size, 
                       pc_I_.col_size);




    Logger::info("[Poisson] - delete temporary block hash.");
    fe_system_.delete_block_hash();

    br_system_ = fe_system_.initialize_block_rack(1, 1);
    
    br_system_.insert_block(dof_field_,         0, 0);

    br_system_.compute_block_offset();
    Logger::info("[Poisson] - initialize block rack: \n"+br_system_.print_block_rack());

    
}





bool Poisson::assemble_system()
{

    // auto& e_data - per-element data shared across integrators.
    // 
    // see assemble_data.h -> Element_Data<phy_dim, ref_dim>
    //
    // Available members:
    //   e             - const Element*, current element
    //   J             - Matrix<phy_dim, ref_dim>, Jacobian at current quad point
    //   inv_J         - Matrix<ref_dim, phy_dim>, inverse/pseudo-inverse of Jacobian
    //   det_J         - double, determinant of Jacobian
    //   b_shape       - Basis_Shape, element geometry
    //   shape_space_1 - const FEM_Space*, trial function space
    //   shape_space_2 - const FEM_Space*, test function space
    //
    //
    // Template parameters:
    //   phy_dim - physical space dimension (1, 2, or 3)
    //   ref_dim - reference element dimension (1, 2, or 3), ref_dim <= phy_dim
    //
    // Note: accessed via auto& in user callbacks. Use e_data.e, e_data.J, etc.

    
    Logger::info("[Poisson] - assemble Poisson-field block matrix.");
    assemble_mat(fe_system_.assemble_mat_data(dof_field_), [&](auto& e_data, auto& mat) {


        Integrator__s_grad_S__grad_S::assemble_element_matrix(1., e_data, mat);

    });

    

    const double pi = CONST::PI;

    
    S_Field_function source_field(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd) {
        return (pi*pi/3)*std::cos(pi*x(0)/3)*std::cos(pi*x(1)/3)*std::cos(pi*x(2)/3);
    });


    Logger::info("[Poisson] - assemble RHS vector for Poisson block.");
    assemble_vec(fe_system_.assemble_vec_data(dof_field_), [&](auto& e_data, auto& vec) {
        Integrator__s__S::assemble_element_vector(source_field, e_data, vec);
    });



   

    


    Logger::info("[Poisson] - build linear system.");
    br_system_.build_linear_system();

    Logger::info("[Poisson] - apply Dirichlet BC to the linear system.");
    bc_.apply_to_system(br_system_);
    
    return true;
}



bool Poisson::assemble_pc_system()
{
    Logger::info("[Poisson - preconditioner] - assemble preconditioner in H1 Omega field.");
    assemble_mat(fe_system_.assemble_mat_data(pc_Q_), [&](auto& e_data, auto& mat) {
        Integrator__s_grad_S__grad_S::assemble_element_matrix(1, e_data, mat);
    });


    Logger::info("[Poisson - preconditioner] - assemble dof mapping from global H1 to omega H1 field.");
    assemble_mat(fe_system_.assemble_mat_data(pc_I_), [&](auto& e_data, auto& mat) {
        Identity_Mapping::direct_mapping(e_data, mat);
    });

    return true;
}



int Poisson::solve_system()
{
    const G_Matrix lhs = br_system_.get_lhs();
    const G_Vector rhs = br_system_.get_rhs();
    const G_Vector x   = br_system_.get_x();

#ifdef LOAD_PETSC
    // for text
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);

    // Set the operator (matrix A)
    KSPSetOperators(ksp, lhs, lhs);

    // Set solver type
    //KSPSetType(ksp, KSPMINRES);
    KSPSetType(ksp, KSPGMRES);

    // Optionally configure the preconditioner (e.g., Jacobi)
    PC pc;
    KSPGetPC(ksp, &pc);    
    //PCSetType(pc, PCHYPRE);
    //PCHYPRESetType(pc, "boomeramg");

    //PCSetType(pc, PCNONE);

    PCSetType(pc, PCSOR);           // SOR family (includes SGS)
    PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP);
    PCSORSetOmega(pc, 1.0);

    // Optionally set tolerances
    KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    // Allow command-line overrides (e.g., -ksp_type, -pc_type)
    KSPSetFromOptions(ksp);

    // for computing condition number
    PetscCall(KSPSetComputeSingularValues(ksp, PETSC_TRUE));

    // Solve
    KSPSolve(ksp, rhs, x);

    PetscReal          smax, smin;
    PetscCall(KSPComputeExtremeSingularValues(ksp, &smax, &smin));
    system_condition_ = (double) smax / smin;

    // Check convergence
    KSPConvergedReason reason;
    PetscInt its;
    KSPGetConvergedReason(ksp, &reason);
    if (reason < 0) {
        PetscPrintf(PETSC_COMM_WORLD, "KSP did not converge: reason %d\n", reason);
    } else {
        KSPGetIterationNumber(ksp, &its);
        PetscPrintf(PETSC_COMM_WORLD, "Converged in %d iterations\n", its);
    }
    petsc_util::petsc_save_ascii_mat(lhs, "poisson_lhs_mat.txt");
    petsc_util::petsc_save_ascii_vec(x, "poisson_x_vec.txt");
    petsc_util::petsc_save_ascii_vec(rhs, "poisson_rhs_vec.txt");
    // Clean up
    KSPDestroy(&ksp);

#else
    Logger::error("[Poisson] - this solver require petsc support!");
#endif

    return its;
}



typedef struct AMSContext_s {
    Mat A, I, Q;
    KSP inner_Q_ksp;     
    Vec tmp;                          // size #edge
    Vec zeta,  kappa;                 // size   #node (Q-space RHS / sol)
    Dirichlet_BC *bc_s;               // for Q (scalar nodal space)
} AMSContext;



static PetscErrorCode AMS_apply(PC pc, Vec r, Vec z)
{
    AMSContext *ctx;
    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));
 
    PetscCall(VecZeroEntries(z));

 
    // single Gauss-Seidel
    PetscCall(MatSOR(ctx->A, r, 1.0, SOR_SYMMETRIC_SWEEP, 0.0, 1, 1, z));

    for(int n=0; n<100; ++n){
 
    //    tmp = r - A z
    PetscCall(MatMult(ctx->A, z, ctx->tmp));
    PetscCall(VecAYPX(ctx->tmp, -1.0, r));
    
    //    zeta = I^T * tmp   (H1)
    PetscCall(MatMultTranspose(ctx->I, ctx->tmp, ctx->zeta));
 
    PetscCall(VecSet(ctx->kappa, 0.0));
 
    ctx->bc_s->apply_to_system(ctx->Q, ctx->zeta, ctx->kappa);
 
    //    Single BoomerAMG V-cycle for each
    PetscCall(KSPSolve(ctx->inner_Q_ksp, ctx->zeta, ctx->kappa));
 
    //    z = z + I*kappa
    PetscCall(MatMultAdd(ctx->I, ctx->kappa, z, z));

    }

 
    // -- post-smoother (continues from updated z)
    PetscCall(MatSOR(ctx->A, r, 1.0, SOR_SYMMETRIC_SWEEP, 0.0, 1, 1, z));
 
    PetscFunctionReturn(0);
}



static PetscErrorCode AMS_destroy(PC pc)
{
    AMSContext *ctx;
    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(pc, (void**)&ctx));
    PetscCall(KSPDestroy(&ctx->inner_Q_ksp));
    PetscCall(VecDestroy(&ctx->tmp));
    PetscCall(VecDestroy(&ctx->zeta));
    PetscCall(VecDestroy(&ctx->kappa));
    PetscCall(PetscFree(ctx));
    PetscFunctionReturn(0);
}



PetscErrorCode solve_AMS(PetscInt& n_iter, double& system_condition,
    Mat A, Vec b, Vec x,
    Mat I, Mat Q,
    Dirichlet_BC *bc_s,
    PetscReal rtol = 1e-10, PetscInt max_iters = PETSC_DEFAULT)
{
    AMSContext *ctx;
    KSP outer;
    PC  outer_pc;
    PetscFunctionBeginUser;
 
    // ---------- Allocate context ----------
    PetscCall(PetscNew(&ctx));
    ctx->A = A;  
    ctx->I = I;  
    ctx->Q = Q;
    ctx->bc_s = bc_s;
 
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
    PetscCall(MatCreateVecs(A, NULL, &ctx->tmp));    // size #edge
    PetscCall(MatCreateVecs(Q, &ctx->kappa, &ctx->zeta)); // size   #node
 
    // ---------- Outer KSP: FGMRES + PCShell(AMS) ----------
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &outer));
    //PetscCall(KSPSetType(outer, KSPGMRES));  
    PetscCall(KSPSetType(outer, KSPCG));    
    PetscCall(KSPSetOperators(outer, A, A));
    PetscCall(KSPSetTolerances(outer, rtol, PETSC_DEFAULT, PETSC_DEFAULT, max_iters));

    /*
    PetscCall(KSPSetNormType(outer, KSP_NORM_UNPRECONDITIONED));  // monitor true ||r||
    {
        PetscViewerAndFormat *vf;
        PetscCall(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,
                                             PETSC_VIEWER_DEFAULT, &vf));
        PetscCall(KSPMonitorSet(outer,
            (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*))KSPMonitorTrueResidual,
            vf,
            (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy));
    }
    */
 
    PetscCall(KSPGetPC(outer, &outer_pc));
    PetscCall(PCSetType(outer_pc, PCSHELL));
    PetscCall(PCShellSetContext(outer_pc, ctx));
    PetscCall(PCShellSetApply  (outer_pc, AMS_apply));
    PetscCall(PCShellSetDestroy(outer_pc, AMS_destroy));
    PetscCall(PCShellSetName   (outer_pc, "AMS"));
 
    PetscCall(KSPSetFromOptions(outer)); 

    // for computing condition number
    PetscCall(KSPSetComputeSingularValues(outer, PETSC_TRUE));

    // ---------- Solve ----------
    PetscCall(VecSet(x, 0.0));
    PetscCall(KSPSolve(outer, b, x));

    PetscReal          smax, smin;
    PetscCall(KSPComputeExtremeSingularValues(outer, &smax, &smin));
    system_condition = (double) smax / smin;
 
    // ---------- Report ----------
    KSPConvergedReason reason;
    PetscInt           iters;
    PetscReal          rnorm;
    PetscCall(KSPGetConvergedReason(outer, &reason));
    PetscCall(KSPGetIterationNumber(outer, &iters));
    PetscCall(KSPGetResidualNorm   (outer, &rnorm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "[AMS] outer GMRES: iters=%" PetscInt_FMT
        ", final ||r||=%.3e, reason=%d\n",
        iters, (double)rnorm, (int)reason));
 
    PetscCall(KSPDestroy(&outer));
    n_iter = iters;

    PetscFunctionReturn(0);
}


int Poisson::solve_pc_system()
{
    int n_iteration_ = 0;
    int iteration_buffer = 0;

#ifdef LOAD_PETSC

    br_system_.extract_block_system();

    Mat A = br_system_.get_block_lhs(0,0);
    Vec b = br_system_.get_block_rhs(0);
    Vec x = br_system_.get_block_x(0);
 
    // pc_P_, pc_G_, pc_L_, pc_Q_ are already assembled.
    // L (= pc_L_) and Q (= pc_Q_) must have Dirichlet BC applied to rows/cols.
    // P and G stay raw.
    
    int n_iteration = 0;
    PetscCall(solve_AMS(n_iteration, system_condition_,
                        A, b, x,
                        pc_I_.mat, pc_Q_.mat,
                        &bc_,
                        1e-10, PETSC_DEFAULT));

    br_system_.assemble_block_x();

#else
    Logger::error("[CurlCurl] - this solver require petsc support!");
#endif

    return n_iteration;
}



scalar_t Poisson::compute_L2_error()
{
    const double pi = CONST::PI;

    // manufactured solution
    S_Field_function u_exact(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd) {
        return std::cos(pi*x(0)/3)*std::cos(pi*x(1)/3)*std::cos(pi*x(2)/3);

    });


    std::string filepath = std::string(TEST_DATA_OUTPUT_DIR) + "/test_poisson.txt";
    std::ofstream file(filepath, std::ios::trunc);  
    

    Logger::info("[Poisson] - compute L2 error.");
    scalar_t l2_error = integrate_element(br_system_, fe_system_, [&](Element_Data<3, 3>& e_data, scalar_t& result) {
        scalar_t local_integral = 0.;

        const std::vector<const FEM_Space*>& space_list = *e_data.space_list;
        const std::vector<std::vector<scalar_t>>& dof_value_list = *e_data.dof_value_list; 

        // compute T - grad Omega
        const std::vector<Integration_Point>& i_p_list = Integration::get_rule(e_data.b_shape, 5);

        Matrix<4, 1> H1_basis;


        scalar_t solved_field = 0.;
        scalar_t solution_field = 0.;


        scalar_t last_solve_f = 0.; // for test
        scalar_t last_exact_f = 0.; // for test

        for(const Integration_Point& i_p : i_p_list)
        {
            solved_field = 0.;
            solution_field = 0.;
            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));

            for(int i=0; i<space_list.size(); ++i)
            {
                const FEM_Space* space = space_list[i];
                const std::vector<scalar_t>& dof_value = dof_value_list[i];
                if(space->get_function_space()==Space::H_1){

                    space->get_basis_s(i_p.coord, H1_basis);

                    for (int j = 0; j < dof_value.size(); ++j) {
                        solved_field += dof_value[j] * H1_basis[j];
                    }
                    
                    
                }else{
                    Logger::error("Poisson::compute_L2_error - integrate_element: detect non H1 space, not expected.");
                }

                last_solve_f = solved_field;  // for test
            }
            
            solution_field = u_exact.eval(i_p.coord, *e_data.e);

            last_solve_f = solved_field;
            last_exact_f = solution_field;

            double diff = solved_field - solution_field;
            local_integral += i_p.weight * abs_det_J * diff * diff;


        }

        result += local_integral;

        Vector<3> node_phys = e_data.physical_point(i_p_list.back().coord);
        size_t property_id = e_data.e->get_property_id();
        file <<"e id: "<<e_data.e->get_id()<<";   p id: "<<property_id<<std::endl;
        file << "phy coord: " << node_phys(0) << ", " << node_phys(1) << ", " << node_phys(2) << ", " << std::endl;
        file << "my result: " <<last_solve_f  << std::endl;
        file << "solution: "  <<last_exact_f<< std::endl;
        file << "==============================" << std::endl;

    });

    file.close();

    return l2_error;
}







using namespace simu;

int main(int argc, char** argv) {

    std::vector<char*> petsc_argv_list;
    petsc_argv_list.push_back(argv[0]);

    std::string mesh_file = "test_cube_0.msh";
    bool enable_preconditioner = false;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if(arg.rfind("--mesh=", 0) == 0) {
            mesh_file = arg.substr(7);
        }else if(arg.rfind("--pc", 0) == 0){
            enable_preconditioner = true;
        }else{
            petsc_argv_list.push_back(argv[i]);  // leave it for PETSc
        }
    }


    char** petsc_argv = petsc_argv_list.data();
    int petsc_argc = petsc_argv_list.size();
    la_kernel::initialize(&petsc_argc, &petsc_argv);

    
    const std::vector<std::pair<std::string, double>> mesh_sweep = {
        {"test_cube_0.geo", 0.500000000},
        {"test_cube_1.geo", 0.353553391},
        {"test_cube_2.geo", 0.250000000},
        {"test_cube_3.geo", 0.176776695},
        {"test_cube_4.geo", 0.125000000},
        {"test_cube_5.geo", 0.088388348},
        {"test_cube_6.geo", 0.062500000},
    };

    std::string file_str = enable_preconditioner ? "/poisson_l2_pc.dat" : "/poisson_l2.dat";

    const std::string dat_path = TEST_DATA_OUTPUT_DIR + file_str;
    std::ofstream l2_convergence(dat_path);
    l2_convergence << "# h                   #element      L2_error      #iteration      condition_number   assemble[s]      solve[s]\n";
    l2_convergence << std::scientific << std::setprecision(15);
    

    for (const auto& [mesh_file, h] : mesh_sweep) {
        scalar_t l2_error;
        size_t n_element = 0;
        int n_iteration = 0;
        double system_condition = 0.;
        double assemble_time = 0.;
        double solving_time = 0.;
        {  
            Logger::start_timer("Loading mesh");
            Mesh_Parser mp(Mesh_Format::GMSH);
            Mesh mesh = mp.load_mesh(SCRIPT_PATH + mesh_file);
            n_element = mesh.get_mesh_elements().size();
            Logger::stop_timer("Loading mesh");

            Logger::start_timer("Initialize Poisson solver");
            Poisson P(mesh);
            Logger::stop_timer("Initialize Poisson solver");

            Logger::start_timer("Assemble Poisson matrix system");
            P.assemble_system();
            if(enable_preconditioner) P.assemble_pc_system();
            assemble_time = Logger::stop_timer("Assemble Poisson matrix system");

            Logger::start_timer("Solve Poisson matrix system");
            if(enable_preconditioner){
                n_iteration = P.solve_pc_system();
            }else{
                n_iteration = P.solve_system();
            }
            solving_time = Logger::stop_timer("Solve Poisson matrix system");

            system_condition = P.get_condition_number();

            Logger::start_timer("Compute L2 error.");
            l2_error = P.compute_L2_error();
            Logger::stop_timer("Compute L2 error.");
        }

        std::ostringstream ss;
        ss << std::scientific << std::setprecision(15) << l2_error;
        Logger::info("[Poisson] h = " + std::to_string(h) + "  L2 error: " + ss.str());

        l2_convergence << h << "  " << n_element
                            << "  " << l2_error
                            << "  " << n_iteration
                            << "  " << system_condition
                            << "  " << assemble_time
                            << "  " << solving_time << "\n";
        l2_convergence.flush();  
    }

    l2_convergence.close();

    la_kernel::finalize();


    Logger::export_to_file("simuEM.log");
    return 0;
}
