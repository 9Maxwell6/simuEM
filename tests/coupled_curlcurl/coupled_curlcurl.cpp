#include "coupled_curlcurl.h"

// for test
#include <fstream>

using namespace simu;

Coupled_CurlCurl::Coupled_CurlCurl(Mesh& mesh, bool enable_preconditioner) : mesh_(mesh), enable_pc_(enable_preconditioner), fe_system_(mesh)
{

    dim_ = mesh.get_mesh_dimension();

    Omega_space_ = H1_Space(mesh.get_mesh_dimension(), 1);
    T_space_1_   = Hcurl_Space(mesh.get_mesh_dimension(),1);

    T_H1_v_    = H1_Space(mesh.get_mesh_dimension(), 1, true, 0);
    T_H1_s_    = H1_Space(mesh.get_mesh_dimension(), 1);

    global_H1  = H1_Space(mesh.get_mesh_dimension(), 1);



    key_true_boundary_ = mesh_.get_keys_true_boundary()[0];  // there must be only one true boundary.
    std::string true_boundary_description = mesh.get_group_description(key_true_boundary_);
    Logger::info("Coupled_CurlCurl - Found simulation boundary: " + true_boundary_description);
    
    

    for(const Key& mesh_key : mesh_.get_keys_internal_surface())
    {
        std::string description = mesh.get_group_description(mesh_key);
        int id = util::extract_last_int(description);
        if(util::a_contains_b(description, {{"Conductor", "Boundary", "1"}, {"Conducting", "Boundary", "1"}})){
            key_conductor_interface_1_ = mesh_key;
            Logger::debug("Coupled_CurlCurl - Found conductor boundary: " + description);
        }
    }
    
    for(const Key& mesh_key : mesh_.get_keys_domain())
    {
        std::string description = mesh.get_group_description(mesh_key);
        int id = util::extract_last_int(description);
        Region new_region{.id=id, .description=description, .r_group=mesh_key };

        if(util::a_contains_b(description, {{"Global"}})) key_global_ = mesh_key;

        if(util::a_contains_b(description, {{"Source, Current"}}))
        {   
            key_source_ = mesh_key;
            mesh_.set_group_property_id(mesh_key, Domain::SOURCE);
            Logger::debug("Coupled_CurlCurl - Found source current: " + description);
        }
        else if(util::a_contains_b(description, {{"Insulator"}, {"Insulating"}}))
        {
            key_insulator_=mesh_key;
            mesh_.set_group_property_id(mesh_key, Domain::EMPTY);
            Logger::debug("Coupled_CurlCurl - Found insulating region: " + description);
        }
        else if(util::a_contains_b(description, {{"Conductor", "1"}, {"Conducting", "1"}}))
        {
            key_conductor_1_ = mesh_key;
            mesh_.set_group_property_id(mesh_key, Domain::CONDUCTOR);
            Logger::debug("Coupled_CurlCurl - Found conductor: " + description);
        }
    }

    // mark all elements in conductors' boundary. (in 3D mesh, they are the 3D elements at boundary layer inside conductor)
    
    std::string description = mesh_.get_group_description(key_conductor_1_);
    key_conductor_interface_layer_1_ =  mesh_.mark_elements(conductor_interface_layer_filter(key_conductor_interface_1_), key_conductor_1_, "Conductor Outerlayer | "+description);

    std::string new_description_conductor_interface_layer = mesh_.get_group_description(key_conductor_interface_layer_1_);
    
    Logger::info("Coupled_CurlCurl - Create conductor interface layer, marked as: " + new_description_conductor_interface_layer + ", #element = " + std::to_string(mesh_.get_element_group(key_conductor_interface_layer_1_).size()));
    Logger::mesh_entity(key_conductor_interface_layer_1_.dim, -1, key_conductor_interface_layer_1_.id, 
                                                                    mesh.get_element_group(key_conductor_interface_layer_1_).size(), 
                                                                    mesh.get_element_geometry_size_group(key_conductor_interface_layer_1_).size(),
                                                                    mesh.get_element_size_group(key_conductor_interface_layer_1_),
                                                                    mesh.get_group_description(key_conductor_interface_layer_1_));
    
    mesh_.set_group_property_id(key_conductor_interface_layer_1_, Domain::CONDUCTOR_OUTER_LAYER);

    key_Omega_inner_boundary_1_ =  mesh_.mark_new_elements(scalar_field_Omega_inner_boundary_filter(key_conductor_interface_1_), dim_-1, key_conductor_interface_layer_1_, "Omega field inner boundary | "+description);
    
    std::string new_description_Omega_field_inner_boundary = mesh_.get_group_description(key_Omega_inner_boundary_1_);
    
    Logger::info("Coupled_CurlCurl - Create Omega-field inner boundary, marked as: " + new_description_Omega_field_inner_boundary + ", #element = " + std::to_string(mesh_.get_element_group(key_Omega_inner_boundary_1_).size()));
    Logger::mesh_entity(key_Omega_inner_boundary_1_.dim, -1, key_Omega_inner_boundary_1_.id, 
                                                                    mesh.get_element_group(key_Omega_inner_boundary_1_).size(), 
                                                                    mesh.get_element_geometry_size_group(key_Omega_inner_boundary_1_).size(),
                                                                    mesh.get_element_size_group(key_Omega_inner_boundary_1_),
                                                                    mesh.get_group_description(key_Omega_inner_boundary_1_));


    Logger::info("Coupled_CurlCurl - Create Omega-field group.");
    key_Omega_ = mesh_.group_union(key_insulator_, key_conductor_interface_layer_1_, "Omega field");
    
    Logger::mesh_entity(key_Omega_.dim, -1, key_Omega_.id, 
                                                     mesh.get_element_group(key_Omega_).size(), 
                                                     mesh.get_element_geometry_size_group(key_Omega_).size(),
                                                     mesh.get_element_size_group(key_Omega_),
                                                     mesh.get_group_description(key_Omega_));




    Logger::info("[Coupled_CurlCurl] - Assign space H1 to Omega-field region.");
    //Block dof_Omega = fe_system_.register_FE_space(field_Omega, new_key_Omega_field);
    dof_Omega_ = fe_system_.register_FE_space(Omega_space_, key_insulator_);
    
    Logger::info("[Coupled_CurlCurl] - register Dirichlet BC to Omega-field.");
    bc_Omega_out_ = fe_system_.register_Dirichlet_BC(dof_Omega_, key_true_boundary_, Dirichlet_Type::HOMOGENEOUS);
    //bc_Omega_in_ = fe_system_.register_Dirichlet_BC(dof_Omega_, key_Omega_inner_boundary_1_, Dirichlet_Type::HOMOGENEOUS);
    
    Logger::block_info(dof_Omega_.id, dof_Omega_.row_offset, dof_Omega_.col_offset, dof_Omega_.row_size, dof_Omega_.col_size);



    Logger::info("[Coupled_CurlCurl] - Assign space Hcurl to T-field region.");
    dof_T_1_ = fe_system_.register_FE_space(T_space_1_, key_conductor_1_);
    
    Logger::info("[Coupled_CurlCurl] - register Dirichlet BC to T-field.");
    bc_T_1_ = fe_system_.register_Dirichlet_BC(dof_T_1_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);

    Logger::block_info(dof_T_1_.id, dof_T_1_.row_offset, dof_T_1_.col_offset, dof_T_1_.row_size, dof_T_1_.col_size);  
    
    //Block dof_T = fe_system_.register_FE_space(field_T, key_conductor[0]);
    
    //Logger::info("[Coupled_CurlCurl] - register Dirichlet BC to T-field.");
    //for(Key& key : key_Omega_field_boundary) bc_Omega_.push_back(_fe_system_.register_Dirichlet_BC(dof_Omega_, key, Dirichlet_Type::HOMOGENEOUS));
    
    



    Logger::info("[Coupled_CurlCurl] - initialize coupling between T-field and Omega-field.");
    dof_coupling_1_ = fe_system_.register_FE_space_coupling(dof_Omega_, dof_T_1_, key_conductor_interface_layer_1_);
    Logger::block_info(dof_coupling_1_.id, 
                       dof_coupling_1_.row_offset, 
                       dof_coupling_1_.col_offset, 
                       dof_coupling_1_.row_size, 
                       dof_coupling_1_.col_size);
    
    //Block dof_coupling = fe_system_.register_FE_space_coupling(dof_T_[i], dof_Omega, key_conductor_interface_layer[i]);



    dof_coupling_tp_1_ = fe_system_.transpose_block(dof_coupling_1_);
    Logger::block_info(dof_coupling_tp_1_.id, 
                       dof_coupling_tp_1_.row_offset, 
                       dof_coupling_tp_1_.col_offset, 
                       dof_coupling_tp_1_.row_size, 
                       dof_coupling_tp_1_.col_size);
    
    //Block dof_coupling_tp = fe_system_.transpose_block(dof_coupling);
    
    

    Logger::info("[Coupled_CurlCurl preconditioner] - initialize edge interpolation matrix. ");
    pc_P_ = fe_system_.register_dual_FE_space(T_space_1_, T_H1_v_, key_conductor_1_, &dof_T_1_, nullptr);
    Logger::block_info(pc_P_.id, 
                       pc_P_.row_offset, 
                       pc_P_.col_offset, 
                       pc_P_.row_size, 
                       pc_P_.col_size);

    Logger::info("[Coupled_CurlCurl preconditioner] - initialize preconditioner in (H1)^3. ");
    pc_L_ = fe_system_.register_dual_FE_space(T_H1_v_, T_H1_v_, key_conductor_1_, &pc_P_, &pc_P_);
    Logger::block_info(pc_L_.id, 
                       pc_L_.row_offset, 
                       pc_L_.col_offset, 
                       pc_L_.row_size, 
                       pc_L_.col_size);

    Logger::info("[Coupled_CurlCurl preconditioner] - initialize preconditioner in H1 global field. ");
    //pc_Q_ = fe_system_.register_FE_space(Omega_space_, key_Omega_, &dof_Omega_);  // TODO debug this
    pc_Q_ = fe_system_.register_dual_FE_space(global_H1, global_H1, key_global_, nullptr, nullptr);
    Logger::block_info(pc_Q_.id, 
                       pc_Q_.row_offset, 
                       pc_Q_.col_offset, 
                       pc_Q_.row_size, 
                       pc_Q_.col_size);

    Logger::info("[Coupled_CurlCurl preconditioner] - initialize descrete gradient matrix from global H1 to T Hcurl field.");
    pc_G_ = fe_system_.register_FE_space_coupling(dof_T_1_, pc_Q_, key_conductor_1_);
    Logger::block_info(pc_G_.id, 
                       pc_G_.row_offset, 
                       pc_G_.col_offset, 
                       pc_G_.row_size, 
                       pc_G_.col_size);

    Logger::info("[Coupled_CurlCurl preconditioner] - initialize dof mapping from global H1 to omega H1 field.");
    pc_I_ = fe_system_.register_FE_space_coupling(dof_Omega_, pc_Q_, key_insulator_);
    Logger::block_info(pc_I_.id, 
                       pc_I_.row_offset, 
                       pc_I_.col_offset, 
                       pc_I_.row_size, 
                       pc_I_.col_size);
    

    pc_bc_O_s_ = fe_system_.register_Dirichlet_BC(pc_Q_, key_true_boundary_, Dirichlet_Type::HOMOGENEOUS);
    pc_bc_O_s_in_ = fe_system_.register_Dirichlet_BC(pc_Q_, key_Omega_inner_boundary_1_, Dirichlet_Type::HOMOGENEOUS);

    pc_bc_T_1_v_ = fe_system_.register_Dirichlet_BC(pc_L_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);
    pc_bc_T_1_s_ = fe_system_.register_Dirichlet_BC(pc_Q_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);

    



    Logger::info("[Coupled_CurlCurl] - delete temporary block hash.");
    fe_system_.delete_block_hash();

    br_system_ = fe_system_.initialize_block_rack(2, 2);
    
    br_system_.insert_block(dof_Omega_,         0, 0);
    br_system_.insert_block(dof_T_1_,           1, 1);
    br_system_.insert_block(dof_coupling_1_,    0, 1);
    br_system_.insert_block(dof_coupling_tp_1_, 1, 0);
    br_system_.compute_block_offset();
    Logger::info("[Coupled_CurlCurl] - initialize block rack: \n"+br_system_.print_block_rack());

    
}





bool Coupled_CurlCurl::assemble_system()
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

    
    Logger::info("[Coupled_CurlCurl] - assemble Omega-field block matrix.");
    assemble_mat(fe_system_.assemble_mat_data(dof_Omega_), [&](auto& e_data, auto& mat) {
        double mu = 1.;
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR) mu = 1.;
        else if(property_id == Domain::EMPTY) mu = 1.;
        //std::cout<<e_data.e->get_property_id()<<std::endl;

        //Integrator__s_S__S::assemble_element_matrix(sigma, e_data, mat);

        Integrator__s_grad_S__grad_S::assemble_element_matrix(-mu, e_data, mat);

    });

    //*
    Logger::info("[Coupled_CurlCurl] - assemble T-field block matrix.");
    assemble_mat(fe_system_.assemble_mat_data(dof_T_1_), [&](auto& e_data, auto& mat) {
        double sigma = 1.;
        double mu = 1.;
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR) { mu = 1.; sigma = 1.; }

        Integrator__s_curl_V__curl_V::assemble_element_matrix(1/sigma, e_data, mat);
        Integrator__s_V__V::assemble_element_matrix(-mu, e_data, mat);

    });

    
    //*/

    //*
    Logger::info("[Coupled_CurlCurl] - assemble coupling block matrix.");
    assemble_mat(fe_system_.assemble_mat_data(dof_coupling_1_), [&](auto& e_data, auto& mat) {
        double mu = 1.;
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR) mu = 1.;
        else if(property_id == Domain::EMPTY) mu = 1.;

        //Integrator__s_grad_S__V::assemble_element_matrix(mu, e_data, mat);
        Integrator__s_V__grad_S::assemble_element_matrix(mu, e_data, mat);

    });
    
    dof_coupling_tp_1_.block_transpose(dof_coupling_1_);
    
    //*/

    const double pi = CONST::PI;

    double sigma = 1.;
    double mu = 1.;

    // 3D vector field.  Hs = ∇×∇×T - T + ∇Ω
    V_Field_function f_T(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        auto sx = std::sin(2*pi*x(0));
        auto sy = std::sin(2*pi*x(1));
        auto sz = std::sin(2*pi*x(2));

        auto cx = std::cos(2*pi*x(0));
        auto cy = std::cos(2*pi*x(1));
        auto cz = std::cos(2*pi*x(2));

        auto Sx = std::sin(pi*x(0));
        auto Sy = std::sin(pi*x(1));
        auto Sz = std::sin(pi*x(2));

        auto Cx = std::cos(pi*x(0));
        auto Cy = std::cos(pi*x(1));
        auto Cz = std::cos(pi*x(2));

        auto ex = std::exp(x(0));
        auto ey = std::exp(x(1));
        auto ez = std::exp(x(2));

        v(0) = 2.0*pi*cx*(ey*sz + ez*sy) - ex*sy*sz*(1.0 - 8*pi*pi) - pi*Sx*Cy*Cz;
        v(1) = 2.0*pi*cy*(ex*sz + ez*sx) - ey*sx*sz*(1.0 - 8*pi*pi) - pi*Cx*Sy*Cz;
        v(2) = 2.0*pi*cz*(ex*sy + ey*sx) - ez*sx*sy*(1.0 - 8*pi*pi) - pi*Cx*Cy*Sz;
    });

    // 3D vector field.  -Hs = T - ∇Ω
    V_Field_function f_Omega_conductor_outer_layer(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        auto sx = std::sin(2*pi*x(0));
        auto sy = std::sin(2*pi*x(1));
        auto sz = std::sin(2*pi*x(2));

        auto Sx = std::sin(pi*x(0));
        auto Sy = std::sin(pi*x(1));
        auto Sz = std::sin(pi*x(2));

        auto Cx = std::cos(pi*x(0));
        auto Cy = std::cos(pi*x(1));
        auto Cz = std::cos(pi*x(2));

        auto ex = std::exp(x(0));
        auto ey = std::exp(x(1));
        auto ez = std::exp(x(2));


        v(0) = ex*sy*sz + pi*Sx*Cy*Cz;
        v(1) = ey*sx*sz + pi*Cx*Sy*Cz;
        v(2) = ez*sx*sy + pi*Cx*Cy*Sz;


    });

    // 3D vector field.  -Hs = -∇Ω
    V_Field_function f_empty(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        if(std::abs(x(0))==0.5 && std::abs(x(1))<=0.5 && std::abs(x(2))<=0.5)
        {
            std::cout<<"WOCAO"<<std::endl;
        }
        v(0) = pi*std::sin(pi*x(0))*std::cos(pi*x(1))*std::cos(pi*x(2));
        v(1) = pi*std::cos(pi*x(0))*std::sin(pi*x(1))*std::cos(pi*x(2));
        v(2) = pi*std::cos(pi*x(0))*std::cos(pi*x(1))*std::sin(pi*x(2));
    });

    Logger::info("[Coupled_CurlCurl] - assemble RHS vector for Omega block.");
    assemble_vec(fe_system_.assemble_vec_data(dof_Omega_), [&](auto& e_data, auto& vec) {    
        Integrator__v__grad_S::assemble_element_vector(f_empty, e_data, vec);                         //   -Hs = -∇Ω  
    });



    Logger::info("[Coupled_CurlCurl] - assemble RHS vector for T-field block.");
    assemble_vec(fe_system_.assemble_vec_data(dof_T_1_), [&](auto& e_data, auto& vec) {
        size_t property_id = e_data.e->get_property_id();
        Integrator__v__V::assemble_element_vector(f_T, e_data, vec);    //   Hs = ∇×∇×T - T + ∇Ω

    
    });

    


    Logger::info("[Coupled_CurlCurl] - build linear system.");
    br_system_.build_linear_system();

    Logger::info("[Coupled_CurlCurl] - apply Dirichlet BC to the linear system.");
    bc_Omega_out_.apply_to_system(br_system_);
    //bc_Omega_in_.apply_to_system(br_system_);
    bc_T_1_.apply_to_system(br_system_);

    //petsc_util::petsc_save_ascii_mat(br_system_.get_lhs(), "lhs.txt");
    
    return true;
}



bool Coupled_CurlCurl::assemble_preconditioner()
{

    Logger::info("[Coupled_CurlCurl - preconditioner] - assemble edge interpolation matrix.");
    assemble_mat(fe_system_.assemble_mat_data(pc_P_), [&](auto& e_data, auto& mat) {
        Interpolator__H1_to_Hcurl::interpolate_element(e_data, mat);
    });

    Logger::info("[Coupled_CurlCurl - preconditioner] - assemble preconditioner in (H1)^3.");
    assemble_mat(fe_system_.assemble_mat_data(pc_L_), [&](auto& e_data, auto& mat) {
        double inv_sigma = 1;
        double mu = 1;
        Integrator__s_grad_V__grad_V::assemble_element_matrix(inv_sigma, e_data, mat);
        Integrator_H1__s_V__V::assemble_element_matrix(-mu, e_data, mat);
    });

    Logger::info("[Coupled_CurlCurl - preconditioner] - assemble preconditioner in H1 global field.");
    assemble_mat(fe_system_.assemble_mat_data(pc_Q_), [&](auto& e_data, auto& mat) {
        double mu = 1;
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR) mu = 1.;
        else if(property_id == Domain::EMPTY) mu = 1.;
        Integrator__s_grad_S__grad_S::assemble_element_matrix(-mu, e_data, mat);
        //Integrator__s_S__S::assemble_element_matrix(mu, e_data, mat);
    });

    Logger::info("[Coupled_CurlCurl - preconditioner] - assemble discrete gradient matrix.");
    assemble_mat(fe_system_.assemble_mat_data(pc_G_), [&](auto& e_data, auto& mat) {
        Interpolator__grad_H1_to_Hcurl::interpolate_element(e_data, mat);
    });

    Logger::info("[Coupled_CurlCurl - preconditioner] - assemble dof mapping from global H1 to omega H1 field.");
    assemble_mat(fe_system_.assemble_mat_data(pc_I_), [&](auto& e_data, auto& mat) {

        Identity_Mapping::direct_mapping(e_data, mat);
    });


    // for debug
    //petsc_util::petsc_save_ascii_mat(pc_P_.mat, "P_mat.txt");
    //petsc_util::petsc_save_ascii_mat(pc_G_.mat, "G_mat.txt");
    //petsc_util::petsc_save_ascii_mat(pc_L_.mat, "L_mat.txt");
    //petsc_util::petsc_save_ascii_mat(pc_G_.mat, "G_T_mat.txt");
    //petsc_util::petsc_save_ascii_mat(pc_I_.mat, "pc_I_mat.txt");


    return true;
}



typedef struct AMS_Context 
{
    Mat P, G, I, L, Q;
    KSP inner_L_ksp, inner_Q_ksp;     // built once, reused per PCApply
    Vec tmp_1, tmp_2;
    Vec rho,   gamma;                 // size 3*nV (L-space RHS / sol)
    Vec zeta,  kappa;                 // size   nV (Q-space RHS / sol)
    Block_Rack* br_system;
    Dirichlet_BC *bc_v;                   // for L (vector / 3D nodal space)
    Dirichlet_BC *bc_s1;                   // for Q (scalar nodal space)
    Dirichlet_BC *bc_s2;                   // for Q (scalar nodal space)
    Dirichlet_BC *bc_s3;                   // for Q (scalar nodal space)
} AMS_Context;

typedef struct AMS_Info 
{
    PetscInt  n_iteration = 0;
    PetscReal condition_number = 0.;
} AMS_Info;



PetscErrorCode AMS_apply(PC pc, Vec r, Vec x)
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

    //PetscScalar value = 1.;
    //PetscInt    ix[1] = {0};

    // pre-smoother.  single Gauss-Seidel sweep
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));
    
    for(int n=0; n<1; ++n){

    std::vector<Vec> block_x(2);
    std::vector<Vec> block_r(2);
    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), x, block_x);
    la_kernel::extract_block_vec(ctx->br_system->get_local_row_size(), r, block_r);


    const Vec xo = block_x[0];
    const Vec xt = block_x[1];

    const Vec ro = block_r[0];
    const Vec rt = block_r[1];
 
    //    tmp_1 = rt - [T]*xt - [coupling_t]*xo
    PetscCall(MatMult(T, xt, ctx->tmp_1));                    // tmp_1 = [T]*xt
    PetscCall(MatMultAdd(T_c, xo, ctx->tmp_1, ctx->tmp_1));   // tmp_1 = tmp_1 + [coupling_t]*xo
    PetscCall(VecAYPX(ctx->tmp_1, -1., rt));                  // tmp_1 = rt - tmp_1

 
    //    rho  = P^T * tmp_1   ((H1)^3)
    PetscCall(MatMultTranspose(ctx->P, ctx->tmp_1, ctx->rho));

    //ctx->bc_v->apply_to_vec(ctx->rho);  
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

    //    zeta = G^T * tmp_1   (H1)
    PetscCall(MatMultTranspose(ctx->G, ctx->tmp_1, ctx->zeta));




    //ctx->bc_s1->apply_to_vec(ctx->zeta);


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


    //ctx->bc_s3->apply_to_vec(ctx->kappa);

    //    xo = xo + I*kappa
    PetscCall(MatMultAdd(ctx->I, ctx->kappa, xo, xo));


    
    //if(n==60) exit(0);

    la_kernel::write_to_vec(block_x, x);



    }
    
    // post-smoother (continues from updated z)
    PetscCall(MatSOR(A, r, 1., SOR_SYMMETRIC_SWEEP, 0., 1, 1, x));

    PetscFunctionReturn(0);
}



PetscErrorCode AMS_destroy(PC pc)
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




PetscErrorCode solve_AMS(
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

    //bc_v->apply_to_mat(L);
    bc_s2->apply_to_mat(Q);
 
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
 
    // ---------- Initialize vectors ----------
    const Mat O   = ctx->br_system->get_block_lhs(0,0);
    const Mat T   = ctx->br_system->get_block_lhs(1,1);
    PetscCall(MatCreateVecs(T, NULL, &ctx->tmp_1));       // size #edge in T-field
    PetscCall(MatCreateVecs(O, NULL, &ctx->tmp_2));       // size #node in Omega-field
    PetscCall(MatCreateVecs(L, &ctx->gamma, &ctx->rho));  // size 3 * #node in T-field
    PetscCall(MatCreateVecs(Q, &ctx->kappa, &ctx->zeta)); // size     #node in global field

    petsc_util::petsc_save_ascii_vec(ctx->zeta, "zeta_initial.txt");
 
    // ---------- Outer KSP + PCShell(AMS) ----------
    const Mat A = ctx->br_system->get_lhs();
    const Vec r = ctx->br_system->get_rhs();
    const Vec x = ctx->br_system->get_x();

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &outer));
    PetscCall(KSPSetType(outer, KSPGMRES));   
    //PetscCall(KSPGMRESSetRestart(outer, 1000));  

    //KSPSetType(outer, KSPCG);

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






bool Coupled_CurlCurl::solve_system()
{
    const G_Matrix lhs = br_system_.get_lhs();
    const G_Vector rhs = br_system_.get_rhs();
    const G_Vector x   = br_system_.get_x();

    petsc_util::petsc_save_ascii_mat(lhs, "lhs_B.txt");
    petsc_util::petsc_save_ascii_vec(rhs, "rhs_B.txt");

    bool successful_flag = false;
#ifdef LOAD_PETSC
    if(!enable_pc_){
        // direct solver
        // Krilov subspace solver
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);

        // Set the operator (lhs matrix)
        KSPSetOperators(ksp, lhs, lhs);

        // Set solver type
        //KSPSetType(ksp, KSPMINRES);
        //KSPSetType(ksp, KSPCG);
        KSPSetType(ksp, KSPGMRES);
        //PetscCall(KSPGMRESSetRestart(ksp, 10000)); 

        // Configure the preconditioner (e.g., Jacobi)
        PC pc;
        KSPGetPC(ksp, &pc);   
        PCSetType(pc, PCHYPRE);
        PCHYPRESetType(pc, "boomeramg");

        //PCSetType(pc, PCNONE);

        //PCSetType(pc, PCSOR);          
        //PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP);
        //PCSORSetOmega(pc, 1.0);

        // Optionally set tolerances
        KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

        // Allow command-line overrides (e.g., -ksp_type, -pc_type)
        KSPSetFromOptions(ksp);

        // for computing condition number
        PetscCall(KSPSetComputeSingularValues(ksp, PETSC_TRUE));

        // Solve
        KSPSolve(ksp, rhs, x);

        // Convergence iterations
        petsc_util::petsc_ksp_convergence(ksp, nullptr, &n_iteration_, "direct solver");

        PetscReal          smax, smin;
        PetscCall(KSPComputeExtremeSingularValues(ksp, &smax, &smin));
        system_condition_ = (double) smax / smin;

        Logger::info("[T-Omega] - total interations: "+std::to_string(n_iteration_)+".");

        // for debug, load matrix to txt file.
        //petsc_util::petsc_save_ascii_mat(lhs, "lhs_mat.txt");
        //petsc_util::petsc_save_ascii_vec(x, "x_vec.txt");
        //petsc_util::petsc_save_ascii_vec(rhs, "rhs_vec.txt");

        // Clean up
        KSPDestroy(&ksp);
    
    }else{

        br_system_.extract_block_system();

        //la_kernel::zero_row_mat(bc_T_1_.get_dof_idx(), 0, pc_G_.mat, NULL, NULL);
        //la_kernel::zero_row_mat(bc_Omega_in_.get_dof_idx(), 0, pc_I_.mat, NULL, NULL);

        AMS_Info ams_info;

        PetscCall(solve_AMS(
            ams_info,
            pc_P_.mat, pc_G_.mat, pc_I_.mat, pc_L_.mat, pc_Q_.mat, 
            &br_system_,
            &pc_bc_T_1_v_, &pc_bc_T_1_s_, &pc_bc_O_s_, &pc_bc_O_s_in_,
            1e-10, PETSC_DEFAULT, true));

        n_iteration_ = ams_info.n_iteration;
        system_condition_ = ams_info.condition_number;
            
    }
    

#else
    Logger::error("[Coupled_CurlCurl] - this solver require petsc support!");
#endif

    return successful_flag;
}







scalar_t Coupled_CurlCurl::compute_L2_error()
{
    const double pi = CONST::PI;

    double sigma = 1.;
    double mu = 1.;

    // manufactured solution
    // 3D vector field.  u = T - ∇Ω
    V_Field_function x_conductor(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        v(0) = std::exp(x(0))*std::sin(2*pi*x(1))*std::sin(2*pi*x(2)) + pi*std::sin(pi*x(0))*std::cos(pi*x(1))*std::cos(pi*x(2));
        v(1) = std::sin(2*pi*x(0))*std::exp(x(1))*std::sin(2*pi*x(2)) + pi*std::cos(pi*x(0))*std::sin(pi*x(1))*std::cos(pi*x(2));
        v(2) = std::sin(2*pi*x(0))*std::sin(2*pi*x(1))*std::exp(x(2)) + pi*std::cos(pi*x(0))*std::cos(pi*x(1))*std::sin(pi*x(2));

    });

    // 3D vector field.  u = -∇Ω
    V_Field_function x_empty(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        v(0) = pi*std::sin(pi*x(0))*std::cos(pi*x(1))*std::cos(pi*x(2));
        v(1) = pi*std::cos(pi*x(0))*std::sin(pi*x(1))*std::cos(pi*x(2));
        v(2) = pi*std::cos(pi*x(0))*std::cos(pi*x(1))*std::sin(pi*x(2));
    });

    std::string debug_file_str = enable_pc_ ? "/pc_debug.txt" : "/debug.txt";

    std::string filepath = std::string(DEBUG_DATA_OUTPUT_DIR) + debug_file_str;
    std::ofstream file(filepath, std::ios::trunc);  
    

    Logger::info("[Coupled_CurlCurl] - compute L2 error.");
    scalar_t l2_error = integrate_element(br_system_, fe_system_, [&](Element_Data<3, 3>& e_data, scalar_t& result) {
        scalar_t local_integral = 0.;

        //if(e_data.e->get_property_id()!=5) return;
        //if(e_data.e->get_property_id()==1) return;

        const std::vector<const FEM_Space*>& space_list = *e_data.space_list;
        const std::vector<std::vector<scalar_t>>& dof_value_list = *e_data.dof_value_list; 

        // compute T - grad Omega
        const std::vector<Integration_Point>& i_p_list = Integration::get_rule(e_data.b_shape, 3);

        Matrix<4, 3> H1_grad_basis;
        Matrix<4, 3> H1_phy_grad_basis;

        Matrix<6, 3> Hcurl_basis;
        Matrix<6, 3> Hcurl_phy_basis;

        Matrix<6, 6> dof_transform;

        Vector<3> solved_field = Vector<3>::Zero();
        Vector<3> solution_field = Vector<3>::Zero();

        Vector<3> last_solve_f; // for test
        Vector<3> last_exact_f; // for test

        Vector<3> temp;

        for(const Integration_Point& i_p : i_p_list)
        {
            solved_field.setZero();
            solution_field.setZero();
            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));
            const Matrix<3, 3>& J_inv = e_data.get_inv_J(i_p.coord); 

            for(int i=0; i<space_list.size(); ++i)
            {
                const FEM_Space* space = space_list[i];
                const std::vector<scalar_t>& dof_value = dof_value_list[i];
                if(space->get_function_space()==Space::H_1){

                    space->get_ED_basis_v(i_p.coord, H1_grad_basis);
                    H1_phy_grad_basis = H1_grad_basis * J_inv;
                    for (int j = 0; j < dof_value.size(); ++j) {
                        solved_field -= dof_value[j] * H1_phy_grad_basis.row(j).transpose();
                    }


                }else if(space->get_function_space()==Space::H_curl){
                    space->dof_transformation(e_data.e->get_node_idx(), dof_transform);
                    space->get_basis_v(i_p.coord, Hcurl_basis);
                    Hcurl_phy_basis = dof_transform*Hcurl_basis * J_inv;
                    for (int j = 0; j < dof_value.size(); ++j) {
                        solved_field += dof_value[j] * Hcurl_phy_basis.row(j).transpose();
                    }
                }
            }

            size_t property_id = e_data.e->get_property_id();
            if(property_id == Domain::CONDUCTOR){
                x_conductor.eval(i_p.coord, *e_data.e, temp);
                solution_field += temp;
            }else if(property_id == Domain::CONDUCTOR_OUTER_LAYER){
                x_conductor.eval(i_p.coord, *e_data.e, temp);
                solution_field += temp;
            }else if(property_id == Domain::EMPTY){
                x_empty.eval(i_p.coord, *e_data.e, temp);
                solution_field += temp;
            }   

            last_solve_f = solved_field;
            last_exact_f = solution_field;

            double diff_sq = (solved_field - solution_field).squaredNorm();
            local_integral += i_p.weight * abs_det_J * diff_sq;

        }

        result += local_integral;

        Vector<3> node_phys = e_data.physical_point(i_p_list.back().coord);
        size_t property_id = e_data.e->get_property_id();
        file <<"e id: "<<e_data.e->get_id()<<";   p id: "<<property_id<<std::endl;
        file << "phy coord: " << node_phys(0) << ", " << node_phys(1) << ", " << node_phys(2) << ", " << std::endl;
        file << "my result: " <<last_solve_f.transpose()  << std::endl;
        file << "solution: "  <<last_exact_f.transpose()<< std::endl;
        file << "==============================" << std::endl;

    });

    file.close();

    return l2_error;

}


void Coupled_CurlCurl::finalize()
{
    la_kernel::destroy_mat(dof_Omega_.mat);
    la_kernel::destroy_mat(dof_T_1_.mat);
    la_kernel::destroy_mat(dof_coupling_1_.mat);
    la_kernel::destroy_mat(dof_coupling_tp_1_.mat);

    la_kernel::destroy_mat(pc_P_.mat);
    la_kernel::destroy_mat(pc_G_.mat);
    la_kernel::destroy_mat(pc_L_.mat);
    la_kernel::destroy_mat(pc_I_.mat);
    la_kernel::destroy_mat(pc_Q_.mat);

    la_kernel::destroy_vec(dof_Omega_.vec);
    la_kernel::destroy_vec(dof_T_1_.vec);
    la_kernel::destroy_vec(dof_coupling_1_.vec);
    la_kernel::destroy_vec(dof_coupling_tp_1_.vec);

    la_kernel::destroy_vec(pc_P_.vec);
    la_kernel::destroy_vec(pc_L_.vec);
    la_kernel::destroy_vec(pc_G_.vec);
    la_kernel::destroy_vec(pc_I_.vec);
    la_kernel::destroy_vec(pc_Q_.vec);
    
    br_system_.finalize();
}

Coupled_CurlCurl::~Coupled_CurlCurl() 
{
    la_kernel::destroy_mat(dof_Omega_.mat);
    la_kernel::destroy_mat(dof_T_1_.mat);
    la_kernel::destroy_mat(dof_coupling_1_.mat);
    la_kernel::destroy_mat(dof_coupling_tp_1_.mat);

    la_kernel::destroy_mat(pc_P_.mat);
    la_kernel::destroy_mat(pc_G_.mat);
    la_kernel::destroy_mat(pc_L_.mat);
    la_kernel::destroy_mat(pc_I_.mat);
    la_kernel::destroy_mat(pc_Q_.mat);

    la_kernel::destroy_vec(dof_Omega_.vec);
    la_kernel::destroy_vec(dof_T_1_.vec);
    la_kernel::destroy_vec(dof_coupling_1_.vec);
    la_kernel::destroy_vec(dof_coupling_tp_1_.vec);

    la_kernel::destroy_vec(pc_P_.vec);
    la_kernel::destroy_vec(pc_G_.vec);
    la_kernel::destroy_vec(pc_L_.vec);
    la_kernel::destroy_vec(pc_I_.vec);
    la_kernel::destroy_vec(pc_Q_.vec);

    br_system_.finalize();

    //if(dof_Omega_.vec) VecDestroy(&dof_Omega_.vec);
    //if(dof_T_1_.vec) VecDestroy(&dof_T_1_.vec);
    //if(dof_coupling_1_.vec) VecDestroy(&dof_coupling_1_.vec);
    //if(dof_coupling_tp_1_.vec) VecDestroy(&dof_coupling_tp_1_.vec);

}





struct Option 
{
    std::string_view long_name; 
    char             short_name; 
    bool             takes_value;
    std::function<void(std::string_view)> handler;
};

using namespace simu;

int main(int argc, char** argv) 
{
    std::string mesh_file    = "test_cc.msh";
    bool l2_test_flag = false;
    bool enable_preconditioner = false;


    const std::vector<Option> options = {
        {"mesh",    'm',   true,  [&](std::string_view v) { mesh_file = std::string(v);    }},
        {"l2-test", '\0',  false, [&](std::string_view)   { l2_test_flag = true;           }},
        {"pc",      '\0',  false, [&](std::string_view)   { enable_preconditioner = true;  }},
        // add more here — one line each
    };


    std::vector<char*> petsc_argv_list{ argv[0] };

    for (int i = 1; i < argc; ++i) 
    {
        std::string_view arg = argv[i];
        bool matched = false;

        for (const auto& opt : options) 
        {
            std::string long_flag  = "--" + std::string(opt.long_name);
            std::string short_flag = opt.short_name ? std::string("-") + opt.short_name: std::string{};

            if (opt.takes_value) {
                if (arg.rfind(long_flag + "=", 0) == 0) {
                    opt.handler(arg.substr(long_flag.size() + 1));
                    matched = true;
                    break;
                }
                if (opt.short_name && arg.rfind(short_flag + "=", 0) == 0) {
                    opt.handler(arg.substr(short_flag.size() + 1));
                    matched = true;
                    break;
                }
            } else {
                if (arg == long_flag || (opt.short_name && arg == short_flag)) {
                    opt.handler("");
                    matched = true;
                    break;
        }}}
        if (!matched) petsc_argv_list.push_back(argv[i]);
    }

    int    petsc_argc = static_cast<int>(petsc_argv_list.size());
    char** petsc_argv = petsc_argv_list.data();
    la_kernel::initialize(&petsc_argc, &petsc_argv);


    if(l2_test_flag)
    {
        const std::vector<std::pair<std::string, double>> mesh_sweep = {
            //{"test_cc_0.geo", 0.500000000},
            {"test_cc_1.geo", 0.353553391},
            {"test_cc_2.geo", 0.250000000},
            {"test_cc_3.geo", 0.176776695},
            {"test_cc_4.geo", 0.125000000},
            {"test_cc_5.geo", 0.088388348},
            {"test_cc_6.geo", 0.062500000},
            {"test_cc_7.geo", 0.044200000},
            {"test_cc_8.geo", 0.031250000},
            
            //{"test_cc_9.geo", 0.022100000},
            //{"test_cc_10.geo", 0.01562500}
        };

        std::string file_str = enable_preconditioner ? "/Coupled_CurlCurl_l2_pc.dat" : "/Coupled_CurlCurl_l2.dat";

        const std::string dat_path = TEST_DATA_OUTPUT_DIR + file_str;

        std::ofstream l2_convergence(dat_path);
        l2_convergence << "# h                    #cell      L2_error           #iteration      condition_number     assemble[s]            solve[s]\n";
        l2_convergence << std::scientific << std::setprecision(15);
        

        for (const auto& [mesh_file, h] : mesh_sweep) {
            petsc_util::petsc_print_memory_usage("iter N start");
            scalar_t l2_error;
            size_t n_element = 0;
            int n_iteration = -1;
            double condition_numbder = 0.;
            double assemble_time = 0.;
            double solving_time = 0.;
            {
                Logger::start_timer("Loading mesh");
                Mesh_Parser mp(Mesh_Format::GMSH);
                Mesh mesh = mp.load_mesh(SCRIPT_PATH + mesh_file);
                n_element = mesh.get_mesh_elements().size();
                Logger::stop_timer("Loading mesh");

                Logger::start_timer("Initialize T-Omega solver");
                Coupled_CurlCurl T_O(mesh, enable_preconditioner);
                Logger::stop_timer("Initialize T-Omega solver");

                Logger::start_timer("Assemble T-Omega matrix system");
                T_O.assemble_system();
                assemble_time += Logger::stop_timer("Assemble T-Omega matrix system");

                if(enable_preconditioner){
                    Logger::start_timer("Assemble T-Omega preconditioner");
                    T_O.assemble_preconditioner();
                    assemble_time += Logger::stop_timer("Assemble T-Omega preconditioner");
                }

                Logger::start_timer("Solve T-Omega matrix system");
                T_O.solve_system();
                solving_time = Logger::stop_timer("Solve T-Omega matrix system");

                Logger::start_timer("Compute L2 error.");
                l2_error = T_O.compute_L2_error();
                Logger::stop_timer("Compute L2 error.");

                n_iteration = T_O.get_n_iteration();
                condition_numbder = T_O.get_system_condition();

                T_O.finalize();
            }

            std::ostringstream ss;
            ss << std::scientific << std::setprecision(15) << l2_error;
            Logger::info("[T-Omega] h = " + std::to_string(h) + "  L2 error: " + ss.str());

            l2_convergence << h << "  " << n_element 
                                << "  " << l2_error 
                                << "  " << n_iteration
                                << "  " << condition_numbder
                                << "  " << assemble_time
                                << "  " << solving_time << "\n";
            l2_convergence.flush();   // persist after every run — a crash on the
                                // finest mesh won't lose the earlier points
            petsc_util::petsc_print_memory_usage("iter N end");
            //PetscLogView(PETSC_VIEWER_STDOUT_WORLD);
        }

        l2_convergence.close();

    }else{
        PetscLogDefaultBegin();
        Logger::start_timer("Load mesh");
        Mesh_Parser mp(Mesh_Format::GMSH);
        Mesh mesh = mp.load_mesh(SCRIPT_PATH + mesh_file);
        Logger::stop_timer("Load mesh");


        Logger::start_timer("Initialize T-Omega solver");
        Coupled_CurlCurl T_O(mesh, enable_preconditioner);
        Logger::stop_timer("Initialize T-Omega solver");

        Logger::start_timer("Assemble T-Omega matrix system");
        T_O.assemble_system();
        Logger::stop_timer("Assemble T-Omega matrix system");

        if(enable_preconditioner){
            Logger::start_timer("Assemble T-Omega preconditioner");
            T_O.assemble_preconditioner();
            Logger::stop_timer("Assemble T-Omega preconditioner");
        }

        Logger::start_timer("Solve T-Omega matrix system");
        T_O.solve_system();
        Logger::stop_timer("Solve T-Omega matrix system");

        Logger::start_timer("Compute L2 error.");
        scalar_t l2_error = T_O.compute_L2_error();
        Logger::stop_timer("Compute L2 error.");

        std::ostringstream ss;
        ss << std::scientific << std::setprecision(15) << l2_error;
        Logger::info("[T-Omega] - test case L2 error: " + ss.str());

        T_O.finalize();
    }

    //PetscLogView(PETSC_VIEWER_STDOUT_WORLD);

    la_kernel::finalize();
    
    Logger::export_to_file("simuEM.log");
    return 0;
}