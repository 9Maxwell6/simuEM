#include "physics/electromagnetism/formulation/T_Omega.h"

// for test
#include <fstream>

using namespace simu;

T_Omega::T_Omega(Mesh& mesh, bool enable_preconditioner) : mesh_(mesh), enable_pc_(enable_preconditioner), fe_system_(mesh)
{

    dim_ = mesh.get_mesh_dimension();

    Omega_space_ = H1_Space(mesh.get_mesh_dimension(), 1);
    T_space_1_   = Hcurl_Space(mesh.get_mesh_dimension(),1);

    T_H1_v_    = H1_Space(mesh.get_mesh_dimension(), 1, true, 0);
    T_H1_s_    = H1_Space(mesh.get_mesh_dimension(), 1);

    global_H1  = H1_Space(mesh.get_mesh_dimension(), 1);



    key_true_boundary_ = mesh_.get_keys_true_boundary()[0];  // there must be only one true boundary.
    std::string true_boundary_description = mesh.get_group_description(key_true_boundary_);
    Logger::info("T_Omega - Found simulation boundary: " + true_boundary_description);
    
    

    for(const Key& mesh_key : mesh_.get_keys_internal_surface())
    {
        std::string description = mesh.get_group_description(mesh_key);
        int id = util::extract_last_int(description);
        if(util::a_contains_b(description, {{"Conductor", "Boundary", "1"}, {"Conducting", "Boundary", "1"}})){
            key_conductor_interface_1_ = mesh_key;
            Logger::debug("T_Omega - Found conductor boundary: " + description);
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
            Logger::debug("T_Omega - Found source current: " + description);
        }
        else if(util::a_contains_b(description, {{"Insulator"}, {"Insulating"}}))
        {
            key_insulator_=mesh_key;
            mesh_.set_group_property_id(mesh_key, Domain::EMPTY);
            Logger::debug("T_Omega - Found insulating region: " + description);
        }
        else if(util::a_contains_b(description, {{"Conductor", "1"}, {"Conducting", "1"}}))
        {
            key_conductor_1_ = mesh_key;
            mesh_.set_group_property_id(mesh_key, Domain::CONDUCTOR);
            Logger::debug("T_Omega - Found conductor: " + description);
        }
    }

    // mark all elements in conductors' boundary. (in 3D mesh, they are the 3D elements at boundary layer inside conductor)
    
    std::string description = mesh_.get_group_description(key_conductor_1_);
    key_conductor_interface_layer_1_ =  mesh_.mark_elements(conductor_interface_layer_filter(key_conductor_interface_1_), key_conductor_1_, "Conductor Outerlayer | "+description);

    std::string new_description_conductor_interface_layer = mesh_.get_group_description(key_conductor_interface_layer_1_);
    
    Logger::info("T_Omega - Create conductor interface layer, marked as: " + new_description_conductor_interface_layer + ", #element = " + std::to_string(mesh_.get_element_group(key_conductor_interface_layer_1_).size()));
    Logger::mesh_entity(key_conductor_interface_layer_1_.dim, -1, key_conductor_interface_layer_1_.id, 
                                                                    mesh.get_element_group(key_conductor_interface_layer_1_).size(), 
                                                                    mesh.get_element_geometry_size_group(key_conductor_interface_layer_1_).size(),
                                                                    mesh.get_element_size_group(key_conductor_interface_layer_1_),
                                                                    mesh.get_group_description(key_conductor_interface_layer_1_));
    
    mesh_.set_group_property_id(key_conductor_interface_layer_1_, Domain::CONDUCTOR_OUTER_LAYER);

    key_Omega_inner_boundary_1_ =  mesh_.mark_new_elements(scalar_field_Omega_inner_boundary_filter(key_conductor_interface_1_), dim_-1, key_conductor_interface_layer_1_, "Omega field inner boundary | "+description);
    
    std::string new_description_Omega_field_inner_boundary = mesh_.get_group_description(key_Omega_inner_boundary_1_);
    
    Logger::info("T_Omega - Create Omega-field inner boundary, marked as: " + new_description_Omega_field_inner_boundary + ", #element = " + std::to_string(mesh_.get_element_group(key_Omega_inner_boundary_1_).size()));
    Logger::mesh_entity(key_Omega_inner_boundary_1_.dim, -1, key_Omega_inner_boundary_1_.id, 
                                                                    mesh.get_element_group(key_Omega_inner_boundary_1_).size(), 
                                                                    mesh.get_element_geometry_size_group(key_Omega_inner_boundary_1_).size(),
                                                                    mesh.get_element_size_group(key_Omega_inner_boundary_1_),
                                                                    mesh.get_group_description(key_Omega_inner_boundary_1_));


    Logger::info("T_Omega - Create Omega-field group.");
    key_Omega_ = mesh_.group_union(key_insulator_, key_conductor_interface_layer_1_, "Omega field");
    
    Logger::mesh_entity(key_Omega_.dim, -1, key_Omega_.id, 
                                                     mesh.get_element_group(key_Omega_).size(), 
                                                     mesh.get_element_geometry_size_group(key_Omega_).size(),
                                                     mesh.get_element_size_group(key_Omega_),
                                                     mesh.get_group_description(key_Omega_));




    Logger::info("[T_Omega] - Assign space H1 to Omega-field region.");
    //Block dof_Omega = fe_system_.register_FE_space(field_Omega, new_key_Omega_field);
    dof_Omega_ = fe_system_.register_FE_space(Omega_space_, key_Omega_);
    
    Logger::info("[T_Omega] - register Dirichlet BC to Omega-field.");
    bc_Omega_out_ = fe_system_.register_Dirichlet_BC(dof_Omega_, key_true_boundary_, Dirichlet_Type::HOMOGENEOUS);
    bc_Omega_in_ = fe_system_.register_Dirichlet_BC(dof_Omega_, key_Omega_inner_boundary_1_, Dirichlet_Type::HOMOGENEOUS);
    
    Logger::block_info(dof_Omega_.id, dof_Omega_.row_offset, dof_Omega_.col_offset, dof_Omega_.row_size, dof_Omega_.col_size);



    Logger::info("[T_Omega] - Assign space Hcurl to T-field region.");
    dof_T_1_ = fe_system_.register_FE_space(T_space_1_, key_conductor_1_);
    
    Logger::info("[T_Omega] - register Dirichlet BC to T-field.");
    bc_T_1_ = fe_system_.register_Dirichlet_BC(dof_T_1_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);

    Logger::block_info(dof_T_1_.id, dof_T_1_.row_offset, dof_T_1_.col_offset, dof_T_1_.row_size, dof_T_1_.col_size);  
    
    //Block dof_T = fe_system_.register_FE_space(field_T, key_conductor[0]);
    
    //Logger::info("[T_Omega] - register Dirichlet BC to T-field.");
    //for(Key& key : key_Omega_field_boundary) bc_Omega_.push_back(_fe_system_.register_Dirichlet_BC(dof_Omega_, key, Dirichlet_Type::HOMOGENEOUS));
    
    



    Logger::info("[T_Omega] - initialize coupling between T-field and Omega-field.");
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
    
    

    Logger::info("[T_Omega preconditioner] - initialize edge interpolation matrix. ");
    pc_P_ = fe_system_.register_dual_FE_space(T_space_1_, T_H1_v_, key_conductor_1_, &dof_T_1_, nullptr);
    Logger::block_info(pc_P_.id, 
                       pc_P_.row_offset, 
                       pc_P_.col_offset, 
                       pc_P_.row_size, 
                       pc_P_.col_size);

    

    Logger::info("[T_Omega preconditioner] - initialize descrete gradient matrix. ");
    pc_G_ = fe_system_.register_dual_FE_space(T_space_1_, T_H1_s_, key_conductor_1_, &dof_T_1_, nullptr);
    Logger::block_info(pc_G_.id, 
                       pc_G_.row_offset, 
                       pc_G_.col_offset, 
                       pc_G_.row_size, 
                       pc_G_.col_size);

    Logger::info("[T_Omega preconditioner] - initialize preconditioner in (H1)^3. ");
    pc_L_ = fe_system_.register_dual_FE_space(T_H1_v_, T_H1_v_, key_conductor_1_, &pc_P_, &pc_P_);
    Logger::block_info(pc_L_.id, 
                       pc_L_.row_offset, 
                       pc_L_.col_offset, 
                       pc_L_.row_size, 
                       pc_L_.col_size);

    Logger::info("[T_Omega preconditioner] - initialize preconditioner in H1. ");
    pc_Q_ = fe_system_.register_dual_FE_space(T_H1_s_, T_H1_s_, key_conductor_1_, &pc_G_, &pc_G_);
    Logger::block_info(pc_Q_.id, 
                       pc_Q_.row_offset, 
                       pc_Q_.col_offset, 
                       pc_Q_.row_size, 
                       pc_Q_.col_size);

    //Logger::info("[T_Omega preconditioner] - initialize descrete gradient matrix for Omega field. ");
    //pc_I_Omega_ = fe_system_.register_FE_space(Omega_space_, T_H1_s_, key_Omega_, &dof_Omega_);
    //Logger::block_info(pc_I_Omega_.id, 
    //                   pc_I_Omega_.row_offset, 
    //                   pc_I_Omega_.col_offset, 
    //                   pc_I_Omega_.row_size, 
    //                   pc_I_Omega_.col_size);

    Logger::info("[T_Omega preconditioner] - initialize preconditioner in H1 Omega field. ");
    //pc_Q_Omega_ = fe_system_.register_FE_space(Omega_space_, key_Omega_, &dof_Omega_);  // TODO debug this
    pc_Q_Omega_ = fe_system_.register_dual_FE_space(Omega_space_, Omega_space_, key_Omega_, &dof_Omega_, &dof_Omega_);
    Logger::block_info(pc_Q_Omega_.id, 
                       pc_Q_Omega_.row_offset, 
                       pc_Q_Omega_.col_offset, 
                       pc_Q_Omega_.row_size, 
                       pc_Q_Omega_.col_size);

    Logger::info("[T_Omega preconditioner] - initialize preconditioner in H1 global field. ");
    //pc_Q_Omega_ = fe_system_.register_FE_space(Omega_space_, key_Omega_, &dof_Omega_);  // TODO debug this
    pc_global_Q_ = fe_system_.register_dual_FE_space(global_H1, global_H1, key_global_, nullptr, nullptr);
    Logger::block_info(pc_global_Q_.id, 
                       pc_global_Q_.row_offset, 
                       pc_global_Q_.col_offset, 
                       pc_global_Q_.row_size, 
                       pc_global_Q_.col_size);

    Logger::info("[T_Omega preconditioner] - initialize descrete gradient matrix from global H1 to T Hcurl field.");
    pc_G_T_ = fe_system_.register_FE_space_coupling(dof_T_1_, pc_global_Q_, key_conductor_1_);
    Logger::block_info(pc_G_T_.id, 
                       pc_G_T_.row_offset, 
                       pc_G_T_.col_offset, 
                       pc_G_T_.row_size, 
                       pc_G_T_.col_size);

    Logger::info("[T_Omega preconditioner] - initialize dof mapping from global H1 to omega H1 field.");
    pc_I_ = fe_system_.register_FE_space_coupling(dof_Omega_, pc_global_Q_, key_Omega_);
    Logger::block_info(pc_I_.id, 
                       pc_I_.row_offset, 
                       pc_I_.col_offset, 
                       pc_I_.row_size, 
                       pc_I_.col_size);
    
    pc_bc_global_ = fe_system_.register_Dirichlet_BC(pc_global_Q_, key_true_boundary_, Dirichlet_Type::HOMOGENEOUS);
    
    pc_bc_T_1_v_ = fe_system_.register_Dirichlet_BC(pc_L_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);
    pc_bc_T_1_s_ = fe_system_.register_Dirichlet_BC(pc_Q_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);

    pc_bc_Omega_out_ = fe_system_.register_Dirichlet_BC(pc_Q_Omega_, key_true_boundary_, Dirichlet_Type::HOMOGENEOUS);

    pc_bc_Omega_in_ = fe_system_.register_Dirichlet_BC(pc_Q_Omega_, key_Omega_inner_boundary_1_, Dirichlet_Type::HOMOGENEOUS);



    Logger::info("[T_Omega] - delete temporary block hash.");
    fe_system_.delete_block_hash();

    br_system_ = fe_system_.initialize_block_rack(2, 2);
    
    br_system_.insert_block(dof_Omega_,         0, 0);
    br_system_.insert_block(dof_T_1_,           1, 1);
    br_system_.insert_block(dof_coupling_1_,    0, 1);
    br_system_.insert_block(dof_coupling_tp_1_, 1, 0);
    br_system_.compute_block_offset();
    Logger::info("[T_Omega] - initialize block rack: \n"+br_system_.print_block_rack());

    
}





bool T_Omega::assemble_system()
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

    
    Logger::info("[T_Omega] - assemble Omega-field block matrix.");
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
    Logger::info("[T_Omega] - assemble T-field block matrix.");
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
    Logger::info("[T_Omega] - assemble coupling block matrix.");
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
        v(0) = pi*std::sin(pi*x(0))*std::cos(pi*x(1))*std::cos(pi*x(2));
        v(1) = pi*std::cos(pi*x(0))*std::sin(pi*x(1))*std::cos(pi*x(2));
        v(2) = pi*std::cos(pi*x(0))*std::cos(pi*x(1))*std::sin(pi*x(2));
    });

    Logger::info("[T_Omega] - assemble RHS vector for Omega block.");
    assemble_vec(fe_system_.assemble_vec_data(dof_Omega_), [&](auto& e_data, auto& vec) {
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR_OUTER_LAYER){
            Integrator__v__grad_S::assemble_element_vector(f_Omega_conductor_outer_layer, e_data, vec);   //  -Hs = T - ∇Ω 
        }else if(property_id == Domain::EMPTY){
            Integrator__v__grad_S::assemble_element_vector(f_empty, e_data, vec);                         //   -Hs = -∇Ω
        }       
    });



    Logger::info("[T_Omega] - assemble RHS vector for T-field block.");
    assemble_vec(fe_system_.assemble_vec_data(dof_T_1_), [&](auto& e_data, auto& vec) {
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR){
            Integrator__v__V::assemble_element_vector(f_T, e_data, vec);    //   Hs = ∇×∇×T - T + ∇Ω
        }else if(property_id == Domain::CONDUCTOR_OUTER_LAYER){
            Integrator__v__V::assemble_element_vector(f_T, e_data, vec);    //   Hs = ∇×∇×T - T + ∇Ω
        }else if(property_id == Domain::EMPTY){
            Logger::error("[T_Omega] - T-field only defined inside conductor!");
        }  
        

    });

    


    Logger::info("[T_Omega] - build linear system.");
    br_system_.build_linear_system();

    Logger::info("[T_Omega] - apply Dirichlet BC to the linear system.");
    bc_Omega_out_.apply_to_system(br_system_);
    bc_Omega_in_.apply_to_system(br_system_);
    bc_T_1_.apply_to_system(br_system_);

    
    return true;
}



bool T_Omega::assemble_preconditioner()
{

    Logger::info("[T_Omega - preconditioner] - assemble edge interpolation matrix.");
    assemble_mat(fe_system_.assemble_mat_data(pc_P_), [&](auto& e_data, auto& mat) {
        Interpolator__H1_to_Hcurl::interpolate_element(e_data, mat);
    });

    Logger::info("[T_Omega - preconditioner] - assemble discrete gradient matrix.");
    assemble_mat(fe_system_.assemble_mat_data(pc_G_), [&](auto& e_data, auto& mat) {
        Interpolator__grad_H1_to_Hcurl::interpolate_element(e_data, mat);
    });

    Logger::info("[T_Omega - preconditioner] - assemble preconditioner in (H1)^3.");
    assemble_mat(fe_system_.assemble_mat_data(pc_L_), [&](auto& e_data, auto& mat) {
        double inv_sigma = 1;
        double mu = 1;
        Integrator__s_grad_V__grad_V::assemble_element_matrix(inv_sigma, e_data, mat);
        Integrator_H1__s_V__V::assemble_element_matrix(-mu, e_data, mat);
    });

    Logger::info("[T_Omega - preconditioner] - assemble preconditioner in H1.");
    assemble_mat(fe_system_.assemble_mat_data(pc_Q_), [&](auto& e_data, auto& mat) {
        double mu = 1;
        Integrator__s_grad_S__grad_S::assemble_element_matrix(-mu, e_data, mat);
        //Integrator__s_S__S::assemble_element_matrix(mu, e_data, mat);
    });

    Logger::info("[T_Omega - preconditioner] - assemble preconditioner in H1 Omega field.");
    assemble_mat(fe_system_.assemble_mat_data(pc_Q_Omega_), [&](auto& e_data, auto& mat) {
        double mu = 1;
        Integrator__s_grad_S__grad_S::assemble_element_matrix(-mu, e_data, mat);
        //Integrator__s_S__S::assemble_element_matrix(mu, e_data, mat);
    });

    Logger::info("[T_Omega - preconditioner] - assemble preconditioner in H1 global field.");
    assemble_mat(fe_system_.assemble_mat_data(pc_global_Q_), [&](auto& e_data, auto& mat) {
        double mu = 1;
        Integrator__s_grad_S__grad_S::assemble_element_matrix(-mu, e_data, mat);
        //Integrator__s_S__S::assemble_element_matrix(mu, e_data, mat);
    });

    Logger::info("[T_Omega - preconditioner] - assemble discrete gradient matrix.");
    assemble_mat(fe_system_.assemble_mat_data(pc_G_T_), [&](auto& e_data, auto& mat) {
        Interpolator__grad_H1_to_Hcurl::interpolate_element(e_data, mat);
    });

    Logger::info("[T_Omega - preconditioner] - assemble dof mapping from global H1 to omega H1 field.");
    assemble_mat(fe_system_.assemble_mat_data(pc_I_), [&](auto& e_data, auto& mat) {

        Identity_Mapping::direct_mapping(e_data, mat);
    });


    // for debug
    petsc_util::petsc_save_ascii_mat(pc_P_.mat, "P_mat.txt");
    petsc_util::petsc_save_ascii_mat(pc_G_.mat, "G_mat.txt");
    petsc_util::petsc_save_ascii_mat(pc_L_.mat, "L_mat.txt");
    petsc_util::petsc_save_ascii_mat(pc_Q_.mat, "Q_mat.txt");
    petsc_util::petsc_save_ascii_mat(pc_G_T_.mat, "G_T_mat.txt");
    petsc_util::petsc_save_ascii_mat(pc_I_.mat, "pc_I_mat.txt");


    return true;
}



bool T_Omega::solve_system()
{
    const G_Matrix lhs = br_system_.get_lhs();
    const G_Vector rhs = br_system_.get_rhs();
    const G_Vector x   = br_system_.get_x();

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
        KSPSetType(ksp, KSPGMRES);

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

        // Solve
        KSPSolve(ksp, rhs, x);

        // Check convergence
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        if (reason < 0) {
            PetscPrintf(PETSC_COMM_WORLD, "KSP did not converge: reason %d\n", reason);
        } else {
            PetscInt its;
            KSPGetIterationNumber(ksp, &its);
            PetscPrintf(PETSC_COMM_WORLD, "Converged in %d iterations\n", its);
            successful_flag = true;
        }

        // for debug, load matrix to txt file.
        petsc_util::petsc_save_ascii_mat(lhs, "lhs_mat.txt");
        petsc_util::petsc_save_ascii_vec(x, "x_vec.txt");
        petsc_util::petsc_save_ascii_vec(rhs, "rhs_vec.txt");

        // Clean up
        KSPDestroy(&ksp);
    
    }else{
        //return solve_pc_system();
        int total_iterations = 0;
        int iteration_buffer = 0;

        // single symmetric Gauss-Seidel sweep to the global system
        PetscCall(MatSOR(lhs, rhs, 1.0, SOR_SYMMETRIC_SWEEP, 0.0, 1, 1, x));

        // switch to block-wise operations
        br_system_.extract_block_system();

        /**
         *  block system view:
         * 
         *      [   Omega    ] [ coupling_o ]  *  [x_o]  =  [r_o]
         *      [ coupling_t ] [     T      ]     [x_t]     [r_t]    
         * 
         */

        Vec tmp;
        // ρ_1 = P^T * ( [r_t] - [T]*[x_t] )
        // ρ_1 = P^T * ( [r_t] - [T]*[x_t] - [coupling_t]*[x_o]) 
        // ζ_1 = G^T * ( [r_t] - [T]*[x_t] - [coupling_t]*[x_o])
        MatCreateVecs(br_system_.get_block_lhs(1,1), NULL, &tmp); 
        MatMult(br_system_.get_block_lhs(1,1), br_system_.get_block_x(1), tmp);         // tmp=[T] * [x_t]
        MatMultAdd(br_system_.get_block_lhs(1,0), br_system_.get_block_x(0), tmp, tmp); // tmp=tmp + [coupling_t]*[x_o]
        VecAYPX(tmp, -1.0, br_system_.get_block_rhs(1));                                // tmp=[r_t] - tmp

        Vec rho_1;
        MatCreateVecs(pc_P_.mat, &rho_1, NULL); 
        MatMultTranspose(pc_P_.mat, tmp, rho_1);      // ρ_1 = P^T * tmp

        Vec zeta_1;
        MatCreateVecs(pc_G_.mat, &zeta_1, NULL); 
        MatMultTranspose(pc_G_.mat, tmp, zeta_1);     // ζ_1 = G^T * tmp

        // ζ_2 = [r_o] - [Omega]*[x_o] - [coupling_o]*[x_t]
        Vec zeta_2;
        MatCreateVecs(br_system_.get_block_lhs(0,0), NULL, &zeta_2); 
        MatMult(br_system_.get_block_lhs(0,0), br_system_.get_block_x(0), zeta_2);            // ζ_2=[Omega] * [x_o]
        MatMultAdd(br_system_.get_block_lhs(0,1), br_system_.get_block_x(1), zeta_2, zeta_2); // ζ_2=ζ_2 + [coupling_o]*[x_t]
        VecAYPX(zeta_2, -1.0, br_system_.get_block_rhs(0));                                   // ζ_2=[r_o] - ζ_2

        // ready to solve
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetType(ksp, KSPGMRES);
        //KSPSetType(ksp, KSPCG);            // need SPD
        PC pc;
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCHYPRE);
        PCHYPRESetType(pc, "boomeramg");     // enable hypre AMG
        KSPSetFromOptions(ksp);

        // L*γ_1  = ρ_1
        KSPSetOperators(ksp, pc_L_.mat, pc_L_.mat);

        Vec gamma_1;
        MatCreateVecs(pc_L_.mat, &gamma_1,  NULL);
        
        pc_bc_T_1_v_.apply_to_system(pc_L_.mat, rho_1, gamma_1);
        
        KSPSolve(ksp, rho_1, gamma_1);       // L*γ_1=ρ_1
        petsc_util::petsc_ksp_convergence(ksp, nullptr, &iteration_buffer, "L*γ_1  = ρ_1");
        total_iterations += iteration_buffer;

        petsc_util::petsc_save_ascii_vec(gamma_1, "gamma_1.txt");

        
        KSPReset(ksp);
        // Q*κ_1  = ζ_1 
        KSPSetOperators(ksp, pc_Q_.mat, pc_Q_.mat);
    
        Vec kappa_1;
        MatCreateVecs(pc_Q_.mat, &kappa_1, NULL);

        pc_bc_T_1_s_.apply_to_system(pc_Q_.mat, zeta_1, kappa_1);

        KSPSolve(ksp, zeta_1, kappa_1);       // Q*κ_1=ζ_1 
        petsc_util::petsc_ksp_convergence(ksp, nullptr, &iteration_buffer, "Q*κ_1  = ζ_1 ");
        total_iterations += iteration_buffer;

        petsc_util::petsc_save_ascii_vec(kappa_1, "kappa_1.txt");


        KSPReset(ksp);
        // Q_Omega*κ_2  = ζ_2
        KSPSetOperators(ksp, pc_Q_Omega_.mat, pc_Q_Omega_.mat);
    
        Vec kappa_2;
        MatCreateVecs(pc_Q_Omega_.mat, &kappa_2, NULL);
        
        pc_bc_Omega_out_.apply_to_system(pc_Q_Omega_.mat, zeta_2, kappa_2);
        pc_bc_Omega_in_.apply_to_system(pc_Q_Omega_.mat, zeta_2, kappa_2);

        KSPSolve(ksp, zeta_2, kappa_2);       // Q*κ_1=ζ_1 
        petsc_util::petsc_ksp_convergence(ksp, nullptr, &iteration_buffer, "Q_Omega*κ_2  = ζ_2");
        total_iterations += iteration_buffer;

        // update
        // [x_t] = [x_t] + P*γ_1 + G*κ_1
        //  1. [x_t] =  [x_t] + P*γ_1
        MatMultAdd(pc_P_.mat, gamma_1, br_system_.get_block_x(1), br_system_.get_block_x(1));
        petsc_util::petsc_save_ascii_vec(br_system_.get_block_x(1), "x_block_gamma1_pc.txt");
        //  2. [x_t] =  [x_t] + G*κ_1
        MatMultAdd(pc_G_.mat, kappa_1, br_system_.get_block_x(1), br_system_.get_block_x(1));
        petsc_util::petsc_save_ascii_vec(br_system_.get_block_x(1), "x_block_kappa1_pc.txt");

        // [x_o] = [x_o] + κ_2 
        petsc_util::petsc_save_ascii_vec(br_system_.get_block_x(0), "x_block_kappa2_0_pc.txt");
        VecAXPY(br_system_.get_block_x(0), 1.0, kappa_2);    // x_o = x_o + kappa_2
        petsc_util::petsc_save_ascii_vec(br_system_.get_block_x(0), "x_block_kappa2_pc.txt");

        br_system_.assemble_block_x();

        // single symmetric Gauss-Seidel sweep to the global system
        MatSOR(lhs, rhs, 1.0, SOR_SYMMETRIC_SWEEP, 0.0, 1, 1, x);

        Logger::info("T-Omega with preconditioner, total interations: "+std::to_string(total_iterations)+".");

        petsc_util::petsc_save_ascii_vec(br_system_.get_x(), "x_vec_pc.txt");

        KSPDestroy(&ksp);
        VecDestroy(&tmp);
        VecDestroy(&rho_1);
        VecDestroy(&zeta_1);
        VecDestroy(&zeta_2);
        VecDestroy(&gamma_1);
        VecDestroy(&kappa_1);
        VecDestroy(&kappa_2);
        MatDestroy(&pc_P_.mat);
        MatDestroy(&pc_G_.mat);
        MatDestroy(&pc_L_.mat);
        MatDestroy(&pc_Q_.mat);
        MatDestroy(&pc_Q_Omega_.mat);

        VecDestroy(&pc_P_.vec);
        VecDestroy(&pc_G_.vec);
        VecDestroy(&pc_L_.vec);
        VecDestroy(&pc_Q_.vec);
        VecDestroy(&pc_Q_Omega_.vec);
        
        

    }
    

#else
    Logger::error("[T_Omega] - this solver require petsc support!");
#endif

    return successful_flag;
}



bool T_Omega::solve_pc_system()
{
    const G_Matrix lhs = br_system_.get_lhs();
    const G_Vector rhs = br_system_.get_rhs();
    const G_Vector x   = br_system_.get_x();

    bool successful_flag = false;

    int total_iterations = 0;
    int iteration_buffer = 0;

    // single symmetric Gauss-Seidel sweep to the global system
    MatSOR(lhs, rhs, 1.0, SOR_SYMMETRIC_SWEEP, 0.0, 1, 1, x);

    // switch to block-wise operations
    br_system_.extract_block_system();

    /**
     *  block system view:
     * 
     *      [   Omega    ] [ coupling_o ]  *  [x_o]  =  [r_o]
     *      [ coupling_t ] [     T      ]     [x_t]     [r_t]    
     * 
     */

    Vec tmp;
    // ρ_1 = P^T * ( [r_t] - [T]*[x_t] - [coupling_t]*[x_o]) 
    // ζ_1 = G^T * ( [r_t] - [T]*[x_t] - [coupling_t]*[x_o])
    MatCreateVecs(br_system_.get_block_lhs(1,1), NULL, &tmp); 
    MatMult(br_system_.get_block_lhs(1,1), br_system_.get_block_x(1), tmp);         // tmp=[T] * [x_t]
    MatMultAdd(br_system_.get_block_lhs(1,0), br_system_.get_block_x(0), tmp, tmp); // tmp=tmp + [coupling_t]*[x_o]
    VecAYPX(tmp, -1.0, br_system_.get_block_rhs(1));                                // tmp=[r_t] - tmp

    Vec rho_1;
    MatCreateVecs(pc_P_.mat, &rho_1, NULL); 
    MatMultTranspose(pc_P_.mat, tmp, rho_1);      // ρ_1 = P^T * tmp


    Vec zeta_1;
    MatCreateVecs(pc_G_T_.mat, &zeta_1, NULL); 
    MatMultTranspose(pc_G_T_.mat, tmp, zeta_1);     // ζ_1 = G^T * tmp


    // ζ_2 = I^T * ([r_o] - [Omega]*[x_o] - [coupling_o]*[x_t])
    MatCreateVecs(br_system_.get_block_lhs(0,0), NULL, &tmp); 
    MatMult(br_system_.get_block_lhs(0,0), br_system_.get_block_x(0), tmp);         // tmp=[T] * [x_t]
    MatMultAdd(br_system_.get_block_lhs(0,1), br_system_.get_block_x(1), tmp, tmp); // tmp=tmp + [coupling_t]*[x_o]
    VecAYPX(tmp, -1.0, br_system_.get_block_rhs(0));                                // tmp=[r_t] - tmp

    Vec zeta_2;
    MatCreateVecs(pc_I_.mat, &zeta_2, NULL); 
    MatMultTranspose(pc_I_.mat, tmp, zeta_2);     // ζ_2 = I^T * tmp

    // ζ_2 = I^T * ([r_o] - [Omega]*[x_o] - [coupling_o]*[x_t])  +  G^T * ( [r_t] - [T]*[x_t] - [coupling_t]*[x_o])
    VecAXPY(zeta_2, 1.0, zeta_1);

    // ready to solve
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    //KSPSetType(ksp, KSPCG);            // need SPD
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCHYPRE);
    PCHYPRESetType(pc, "boomeramg");     // enable hypre AMG
    KSPSetFromOptions(ksp);


    // L*γ_1  = ρ_1
    KSPSetOperators(ksp, pc_L_.mat, pc_L_.mat);

    Vec gamma_1;
    MatCreateVecs(pc_L_.mat, &gamma_1,  NULL);
    
    pc_bc_T_1_v_.apply_to_system(pc_L_.mat, rho_1, gamma_1);
    
    KSPSolve(ksp, rho_1, gamma_1);       // L*γ_1=ρ_1
    petsc_util::petsc_ksp_convergence(ksp, nullptr, &iteration_buffer, "L*γ_1  = ρ_1");
    total_iterations += iteration_buffer;

    petsc_util::petsc_save_ascii_vec(gamma_1, "gamma_1.txt");



    // global 
    KSPReset(ksp);
    // Q_global * κ_global  = ζ_2_global
    KSPSetOperators(ksp, pc_global_Q_.mat, pc_global_Q_.mat);

    Vec kappa_2;
    MatCreateVecs(pc_global_Q_.mat, &kappa_2, NULL);
    
    pc_bc_global_.apply_to_system(pc_global_Q_.mat, zeta_2, kappa_2);

    KSPSolve(ksp, zeta_2, kappa_2);       // Q_global * κ_global  = ζ_2_global



    petsc_util::petsc_save_ascii_mat(pc_Q_Omega_.mat, "pc_Q_Omega_.txt");


    // [x_t] = [x_t] + P*γ_1 + G*κ_1
    // [x_0] = [x_o] + I*κ_2 
    // 
    // => [x_t] = [x_t] + P*γ_1 + G*kappa_global
    MatMultAdd(pc_P_.mat, gamma_1, br_system_.get_block_x(1), br_system_.get_block_x(1));
    MatMultAdd(pc_G_T_.mat, kappa_2, br_system_.get_block_x(1), br_system_.get_block_x(1));

    //    [x_o] = [x_o] + I*kappa_global
    MatMultAdd(pc_I_.mat, kappa_2, br_system_.get_block_x(0), br_system_.get_block_x(0));
    //    [x_global] = [x_global] + I*κ_2   


    br_system_.assemble_block_x();

    // single symmetric Gauss-Seidel sweep to the global system
    MatSOR(lhs, rhs, 1.0, SOR_SYMMETRIC_SWEEP, 0.0, 1, 1, x);

    Logger::info("T-Omega with preconditioner, total interations: "+std::to_string(total_iterations)+".");

    petsc_util::petsc_save_ascii_vec(br_system_.get_x(), "x_vec_pc.txt");

    KSPDestroy(&ksp);
    VecDestroy(&tmp);
    VecDestroy(&rho_1);
    VecDestroy(&zeta_1);
    VecDestroy(&zeta_2);
    VecDestroy(&gamma_1);
    VecDestroy(&kappa_2);
    return true;
}



scalar_t T_Omega::compute_L2_error()
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
    

    Logger::info("[T_Omega] - compute L2 error.");
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


void T_Omega::finalize()
{
    la_kernel::destroy_mat(dof_Omega_.mat);
    la_kernel::destroy_mat(dof_T_1_.mat);
    la_kernel::destroy_mat(dof_coupling_1_.mat);
    la_kernel::destroy_mat(dof_coupling_tp_1_.mat);

    la_kernel::destroy_mat(pc_P_.mat);
    la_kernel::destroy_mat(pc_G_.mat);
    la_kernel::destroy_mat(pc_L_.mat);
    la_kernel::destroy_mat(pc_Q_.mat);

    la_kernel::destroy_mat(pc_Q_Omega_.mat);
    la_kernel::destroy_mat(pc_G_T_.mat);
    la_kernel::destroy_mat(pc_I_.mat);
    la_kernel::destroy_mat(pc_global_Q_.mat);

    la_kernel::destroy_vec(dof_Omega_.vec);
    la_kernel::destroy_vec(dof_T_1_.vec);
    la_kernel::destroy_vec(dof_coupling_1_.vec);
    la_kernel::destroy_vec(dof_coupling_tp_1_.vec);

    la_kernel::destroy_vec(pc_P_.vec);
    la_kernel::destroy_vec(pc_G_.vec);
    la_kernel::destroy_vec(pc_L_.vec);
    la_kernel::destroy_vec(pc_Q_.vec);

    la_kernel::destroy_vec(pc_Q_Omega_.vec);
    la_kernel::destroy_vec(pc_G_T_.vec);
    la_kernel::destroy_vec(pc_I_.vec);
    la_kernel::destroy_vec(pc_global_Q_.vec);
    
    br_system_.finalize();
}

T_Omega::~T_Omega() 
{
    la_kernel::destroy_mat(dof_Omega_.mat);
    la_kernel::destroy_mat(dof_T_1_.mat);
    la_kernel::destroy_mat(dof_coupling_1_.mat);
    la_kernel::destroy_mat(dof_coupling_tp_1_.mat);

    la_kernel::destroy_mat(pc_P_.mat);
    la_kernel::destroy_mat(pc_G_.mat);
    la_kernel::destroy_mat(pc_L_.mat);
    la_kernel::destroy_mat(pc_Q_.mat);

    la_kernel::destroy_mat(pc_Q_Omega_.mat);
    la_kernel::destroy_mat(pc_G_T_.mat);
    la_kernel::destroy_mat(pc_I_.mat);
    la_kernel::destroy_mat(pc_global_Q_.mat);

    la_kernel::destroy_vec(dof_Omega_.vec);
    la_kernel::destroy_vec(dof_T_1_.vec);
    la_kernel::destroy_vec(dof_coupling_1_.vec);
    la_kernel::destroy_vec(dof_coupling_tp_1_.vec);

    la_kernel::destroy_vec(pc_P_.vec);
    la_kernel::destroy_vec(pc_G_.vec);
    la_kernel::destroy_vec(pc_L_.vec);
    la_kernel::destroy_vec(pc_Q_.vec);

    la_kernel::destroy_vec(pc_Q_Omega_.vec);
    la_kernel::destroy_vec(pc_G_T_.vec);
    la_kernel::destroy_vec(pc_I_.vec);
    la_kernel::destroy_vec(pc_global_Q_.vec);

    br_system_.finalize();

    //if(dof_Omega_.vec) VecDestroy(&dof_Omega_.vec);
    //if(dof_T_1_.vec) VecDestroy(&dof_T_1_.vec);
    //if(dof_coupling_1_.vec) VecDestroy(&dof_coupling_1_.vec);
    //if(dof_coupling_tp_1_.vec) VecDestroy(&dof_coupling_tp_1_.vec);

}
