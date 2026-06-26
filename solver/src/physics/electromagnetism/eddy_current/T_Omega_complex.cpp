#include "physics/electromagnetism/eddy_current/T_Omega_complex.h"

using namespace simu;


T_Omega_c::T_Omega_c(Mesh& mesh, int pc_alg) : mesh_(mesh), pc_alg_(pc_alg), fe_system_(mesh)
{

    dim_ = mesh.get_mesh_dimension();



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


}


bool T_Omega_c::initialize_system()
{
    O_space_    = H1_Space(mesh_.get_mesh_dimension(), 1);
    T_space_    = Hcurl_Space(mesh_.get_mesh_dimension(),1);

    Logger::info("[T_Omega] - register [Ki_r].");
    Ki_r_ = fe_system_.register_FE_space(O_space_, key_Omega_);
    Logger::block_info(Ki_r_.id, Ki_r_.row_offset, Ki_r_.col_offset, Ki_r_.row_size, Ki_r_.col_size);

    Logger::info("[T_Omega] - register [Ki_c].");
    Ki_c_ = fe_system_.register_FE_space(O_space_, key_Omega_, &Ki_r_);
    Logger::block_info(Ki_c_.id, Ki_c_.row_offset, Ki_c_.col_offset, Ki_c_.row_size, Ki_c_.col_size);

    
    Logger::info("[T_Omega] - register Dirichlet BC to Omega-field.");
    bc_O_o_r_ = fe_system_.register_Dirichlet_BC(Ki_r_, key_true_boundary_, Dirichlet_Type::HOMOGENEOUS);
    bc_O_i_r_ = fe_system_.register_Dirichlet_BC(Ki_r_, key_Omega_inner_boundary_1_, Dirichlet_Type::HOMOGENEOUS);
    bc_O_o_c_ = fe_system_.register_Dirichlet_BC(Ki_c_, key_true_boundary_, Dirichlet_Type::HOMOGENEOUS);
    bc_O_i_c_ = fe_system_.register_Dirichlet_BC(Ki_c_, key_Omega_inner_boundary_1_, Dirichlet_Type::HOMOGENEOUS);
    
    


    Logger::info("[T_Omega] - register [Kc_r].");
    Kc_r_ = fe_system_.register_FE_space(T_space_, key_conductor_1_);
    Logger::block_info(Kc_r_.id, Kc_r_.row_offset, Kc_r_.col_offset, Kc_r_.row_size, Kc_r_.col_size);  

    Logger::info("[T_Omega] - register [Kc_c].");
    Kc_c_ = fe_system_.register_FE_space(T_space_, key_conductor_1_, &Kc_r_);
    Logger::block_info(Kc_c_.id, Kc_c_.row_offset, Kc_c_.col_offset, Kc_c_.row_size, Kc_c_.col_size); 

    Logger::info("[T_Omega] - register [Mc_r].");
    Mc_r_ = fe_system_.register_FE_space(T_space_, key_conductor_1_);
    Logger::block_info(Mc_r_.id, Mc_r_.row_offset, Mc_r_.col_offset, Mc_r_.row_size, Mc_r_.col_size);  

    Logger::info("[T_Omega] - register [Mc_c].");
    Mc_c_ = fe_system_.register_FE_space(T_space_, key_conductor_1_, &Mc_r_);
    Logger::block_info(Mc_c_.id, Mc_c_.row_offset, Mc_c_.col_offset, Mc_c_.row_size, Mc_c_.col_size);  

    
    Logger::info("[T_Omega] - register Dirichlet BC to T");
    bc_T_1_r_ = fe_system_.register_Dirichlet_BC(Mc_r_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);
    bc_T_1_c_ = fe_system_.register_Dirichlet_BC(Mc_c_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);

    
    
    


    Logger::info("[T_Omega] - register [X__r].");
    X__r_ = fe_system_.register_FE_space_coupling(Kc_r_, Ki_r_, key_conductor_interface_layer_1_);
    Logger::block_info(X__r_.id, X__r_.row_offset, X__r_.col_offset, X__r_.row_size, X__r_.col_size);

    Logger::info("[T_Omega] - register [X__c].");
    X__c_ = fe_system_.register_FE_space_coupling(Kc_c_, Ki_c_, key_conductor_interface_layer_1_);
    Logger::block_info(X__c_.id, X__c_.row_offset, X__c_.col_offset, X__c_.row_size, X__c_.col_size);
    
    //Block dof_coupling = fe_system_.register_FE_space_coupling(dof_T_[i], dof_Omega, key_conductor_interface_layer[i]);


    Logger::info("[T_Omega] - register [Xt_r].");
    Xt_r_ = fe_system_.transpose_block(X__r_);
    Logger::block_info(Xt_r_.id,  Xt_r_.row_offset,  Xt_r_.col_offset,  Xt_r_.row_size,  Xt_r_.col_size);

    Logger::info("[T_Omega] - register [Xt_c].");
    Xt_c_ = fe_system_.transpose_block(X__c_);
    Logger::block_info(Xt_c_.id,  Xt_c_.row_offset,  Xt_c_.col_offset,  Xt_c_.row_size,  Xt_c_.col_size);


    Logger::info("[T_Omega] - register source term [Sc_r].");
    Sc_r_ = fe_system_.register_FE_space(T_space_, key_conductor_1_, &Mc_r_);   // Re(Sc)
    Logger::block_info(Sc_r_.id,  Sc_r_.row_offset,  Sc_r_.col_offset,  Sc_r_.row_size,  Sc_r_.col_size);

    Logger::info("[T_Omega] - register source term [Sc_c].");
    Sc_c_ = fe_system_.register_FE_space(T_space_, key_conductor_1_, &Mc_c_);   // Im(Sc)
    Logger::block_info(Sc_c_.id,  Sc_c_.row_offset,  Sc_c_.col_offset,  Sc_c_.row_size,  Sc_c_.col_size);
    
    Logger::info("[T_Omega] - register source term [Si_r].");
    Si_r_ = fe_system_.register_FE_space(O_space_, key_Omega_, &Ki_r_);   // Re(Si)
    Logger::block_info(Si_r_.id,  Si_r_.row_offset,  Si_r_.col_offset,  Si_r_.row_size,  Si_r_.col_size);

    Logger::info("[T_Omega] - register source term [Si_c].");
    Si_c_ = fe_system_.register_FE_space(O_space_, key_Omega_, &Ki_c_);   // Im(Si)
    Logger::block_info(Si_c_.id,  Si_c_.row_offset,  Si_c_.col_offset,  Si_c_.row_size,  Si_c_.col_size);


    if(pc_alg_!=0) initialize_pc();

    Logger::info("[T_Omega] - delete temporary block hash.");
    fe_system_.delete_block_hash();


    br_system_ = fe_system_.initialize_block_rack(4, 4);
    
    br_system_.insert_block(Mc_r_,         0, 0);
    br_system_.insert_block(X__r_,         0, 1);
    br_system_.insert_block(Kc_r_,         0, 2);
    br_system_.insert_block(Xt_r_,         1, 0);
    br_system_.insert_block(Ki_r_,         1, 1);
    br_system_.insert_block(Kc_c_,         2, 0);
    br_system_.insert_block(Mc_c_,         2, 2);
    br_system_.insert_block(X__c_,         2, 3);
    br_system_.insert_block(Xt_c_,         3, 2);
    br_system_.insert_block(Ki_c_,         3, 3);

    br_system_.insert_vec(Sc_r_,         0);
    br_system_.insert_vec(Si_r_,         1);
    br_system_.insert_vec(Sc_c_,         2);
    br_system_.insert_vec(Si_c_,         3);

    br_system_.compute_block_offset();
    Logger::info("[T_Omega] - initialize system block rack: \n"+br_system_.print_block_rack());

    

    return true;
}



bool T_Omega_c::initialize_pc()
{
    T_space_1_H1_v_ = H1_Space(mesh_.get_mesh_dimension(), 1, true, 0);
    T_space_1_H1_s_ = H1_Space(mesh_.get_mesh_dimension(), 1);

    global_H1_  = H1_Space(mesh_.get_mesh_dimension(), 1);

    Logger::info("[T_Omega preconditioner] - register [P]. ");
    pc_P_ = fe_system_.register_dual_FE_space(T_space_, T_space_1_H1_v_, key_conductor_1_, &Mc_r_, nullptr);
    Logger::block_info(pc_P_.id, pc_P_.row_offset, pc_P_.col_offset, pc_P_.row_size, pc_P_.col_size);


    Logger::info("[T_Omega preconditioner] - register [L+M] ");
    pc_LpM_ = fe_system_.register_dual_FE_space(T_space_1_H1_v_, T_space_1_H1_v_, key_conductor_1_, &pc_P_, &pc_P_);
    Logger::block_info(pc_LpM_.id, pc_LpM_.row_offset, pc_LpM_.col_offset, pc_LpM_.row_size,  pc_LpM_.col_size);

    if(pc_alg_==1){
        // decoupled pc
        Logger::info("[T_Omega preconditioner - alg1] - register [Qc] ");
        pc_Qc_ = fe_system_.register_dual_FE_space(T_space_1_H1_s_, T_space_1_H1_s_, key_conductor_1_, nullptr, nullptr);
        Logger::block_info(pc_Qc_.id, pc_Qc_.row_offset, pc_Qc_.col_offset, pc_Qc_.row_size, pc_Qc_.col_size);

        Logger::info("[T_Omega preconditioner - alg1] - register [Qc] ");
        pc_Qi_ = fe_system_.register_dual_FE_space(O_space_, O_space_, key_Omega_, nullptr, nullptr);
        Logger::block_info(pc_Qc_.id, pc_Qc_.row_offset, pc_Qc_.col_offset, pc_Qc_.row_size, pc_Qc_.col_size);
    }

    if(pc_alg_==2){
        // global gradient pc
        Logger::info("[T_Omega preconditioner - alg2] - create temporary global H1 dof table. ");
        tmp_global_ = fe_system_.register_FE_space(global_H1_, key_global_);
        Logger::info("[T_Omega preconditioner - alg2] - register [G].");
        pc_G_ = fe_system_.register_FE_space_coupling(Mc_r_, tmp_global_, key_conductor_1_);
        Logger::block_info(pc_G_.id, pc_G_.row_offset, pc_G_.col_offset, pc_G_.row_size, pc_G_.col_size);

        Logger::info("[T_Omega preconditioner - alg2] - register [I].");
        pc_I_ = fe_system_.register_FE_space_coupling(Ki_r_, tmp_global_, key_insulator_);
        Logger::block_info(pc_I_.id, pc_I_.row_offset, pc_I_.col_offset, pc_I_.row_size, pc_I_.col_size);
    }

    

    if(pc_alg_==1 || pc_alg_==3 || pc_alg_==4){
        Logger::info("[T_Omega preconditioner - alg1/3/4] - register [G].");
        pc_G_ = fe_system_.register_FE_space_coupling(Mc_r_, pc_Qc_, key_conductor_1_);
        Logger::block_info(pc_G_.id, pc_G_.row_offset, pc_G_.col_offset, pc_G_.row_size, pc_G_.col_size);

        Logger::info("[T_Omega preconditioner - alg1/3/4] - register [I].");
        pc_I_ = fe_system_.register_FE_space_coupling(Ki_r_, pc_Qi_, key_insulator_);
        Logger::block_info(pc_I_.id, pc_I_.row_offset, pc_I_.col_offset, pc_I_.row_size, pc_I_.col_size);
    }


    Logger::info("[T_Omega] - register boundary T in H1.");
    bc_T_1_H1_s_ = fe_system_.register_Dirichlet_BC(pc_Qc_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);

    Logger::info("[T_Omega] - register boundary T in (H1)^3.");
    bc_T_1_H1_v_ = fe_system_.register_Dirichlet_BC(pc_LpM_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);


    return true;
}




bool T_Omega_c::assemble_system()
{
    //coefficient:
    double mu_conductor_ = 1.;
    double mu_insulator_ = 1.;
    double sigma_conductor_ = 1.;
    double omega = 1.;

    Logger::info("[T_Omega] - assemble [Kc].");
    assemble_mat(fe_system_.assemble_mat_data(Kc_r_), [&](auto& e_data, auto& mat) {
        double sigma = 0;
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR) sigma = sigma_conductor_;
        else if(property_id == Domain::CONDUCTOR_OUTER_LAYER) sigma = sigma_conductor_;

        Integrator__s_curl_V__curl_V::assemble_element_matrix(1/sigma, e_data, mat);
    });

    Logger::info("[T_Omega] - set [_Kc = -Kc].");
    la_kernel::duplicate_mat(Kc_r_.mat, Kc_c_.mat);
    la_kernel::scale_mat(-1.0, Kc_c_.mat);

    Logger::info("[T_Omega] - assemble [Mc].");
    assemble_mat(fe_system_.assemble_mat_data(Mc_r_), [&](auto& e_data, auto& mat) {
        double mu = 0;
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR) mu = mu_conductor_;
        else if(property_id == Domain::CONDUCTOR_OUTER_LAYER) mu = mu_conductor_;

        Integrator__s_V__V::assemble_element_matrix(mu, e_data, mat);
    });

    Mc_c_.mat = Mc_r_.mat;

    Logger::info("[T_Omega] - assemble [Ki].");
    assemble_mat(fe_system_.assemble_mat_data(Ki_r_), [&](auto& e_data, auto& mat) {
        double mu = 0;
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR_OUTER_LAYER) mu = mu_conductor_;
        else if(property_id == Domain::EMPTY) mu = mu_insulator_;

        Integrator__s_grad_S__grad_S::assemble_element_matrix(mu, e_data, mat);
    });

    Ki_c_.mat = Ki_r_.mat;


    Logger::info("[T_Omega] - assemble [X__r_].");
    assemble_mat(fe_system_.assemble_mat_data(X__r_), [&](auto& e_data, auto& mat) {
        double mu = 0;
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR_OUTER_LAYER) mu = mu_conductor_;

        Integrator__s_grad_S__V::assemble_element_matrix(-mu, e_data, mat);
    });

    Logger::info("[T_Omega] - assemble [X__c_].");
    X__c_.mat = X__r_.mat;

    Logger::info("[T_Omega] - assemble [Xt_r_].");
    la_kernel::copy_transpose(X__r_.mat, Xt_r_.mat);

    Logger::info("[T_Omega] - assemble [Xt_c_].");
    Xt_c_.mat = Xt_r_.mat;



    // rhs linear form integrator
    const double pi = CONST::PI;


    // ar = [1, 2, -1]   ai = [-2, 1, 1]

    V_Field_function ball_field_Sc_r(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        double rr = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
        double q = 0.25 - 2*rr;
        double phi = (0.25-rr)*(0.25-rr);
        v(0) = -8*(( -2*x(0)+x(1)+x(2) )*x(0) - 2*q) - phi - 2*x(0);
        v(1) = -8*(( -2*x(0)+x(1)+x(2) )*x(1) + q) - 2*phi - 2*x(1);
        v(2) = -8*(( -2*x(0)+x(1)+x(2) )*x(2) + q) + phi - 2*x(2);
    });

    V_Field_function ball_field_Sc_r_neg(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        double rr = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
        double q = 0.25 - 2*rr;
        double phi = (0.25-rr)*(0.25-rr);
        v(0) = -8*(( -2*x(0)+x(1)+x(2) )*x(0) - 2*q) - phi - 2*x(0);
        v(1) = -8*(( -2*x(0)+x(1)+x(2) )*x(1) + q) - 2*phi - 2*x(1);
        v(2) = -8*(( -2*x(0)+x(1)+x(2) )*x(2) + q) + phi - 2*x(2);

        v(0) = -v(0);
        v(1) = -v(1);
        v(2) = -v(2);
    });

    V_Field_function ball_field_Sc_c(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        double rr = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
        double q = 0.25 - 2*rr;
        double phi = (0.25-rr)*(0.25-rr);
        v(0) = 8*((x(0)+2*x(1)-x(2))*x(0) + q) + 2*phi + 9./4. - rr - 2*x(0)*x(0);
        v(1) = 8*((x(0)+2*x(1)-x(2))*x(1) + 2*q) - phi - 2*x(0)*x(1);
        v(2) = 8*((x(0)+2*x(1)-x(2))*x(2) - q) - phi - 2*x(0)*x(2);
    });

    V_Field_function ball_field_Sc_c_neg(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        double rr = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
        double q = 0.25 - 2*rr;
        double phi = (0.25-rr)*(0.25-rr);
        v(0) = 8*((x(0)+2*x(1)-x(2))*x(0) + q) + 2*phi + 9./4. - rr - 2*x(0)*x(0);
        v(1) = 8*((x(0)+2*x(1)-x(2))*x(1) + 2*q) - phi - 2*x(0)*x(1);
        v(2) = 8*((x(0)+2*x(1)-x(2))*x(2) - q) - phi - 2*x(0)*x(2);

        v(0) = -v(0);
        v(1) = -v(1);
        v(2) = -v(2);
    });

    V_Field_function ball_field_Si_r(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        v(0) = -2*x(0);
        v(1) = -2*x(1);
        v(2) = -2*x(2);
    });

    V_Field_function ball_field_Si_c(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        double rr = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
        v(0) = 9./4. - rr - 2*x(0)*x(0);
        v(1) = -2*x(0)*x(1);
        v(2) = -2*x(0)*x(2);
    });




    Logger::info("[T_Omega] - assemble RHS vector (real part in conductor) [Sc_r].");
    assemble_vec(fe_system_.assemble_vec_data(Sc_r_), [&](auto& e_data, auto& vec) {
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR){
            Integrator__v__V::assemble_element_vector(ball_field_Sc_r_neg, e_data, vec);  

        }else if(property_id == Domain::CONDUCTOR_OUTER_LAYER){
            Integrator__v__V::assemble_element_vector(ball_field_Sc_r_neg, e_data, vec);  
        }       
    });


    Logger::info("[T_Omega] - assemble RHS vector (real part in omega) [Si_r].");
    assemble_vec(fe_system_.assemble_vec_data(Si_r_), [&](auto& e_data, auto& vec) {
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR_OUTER_LAYER){
            Integrator__v__grad_S::assemble_element_vector(ball_field_Sc_r, e_data, vec);   

        }else if(property_id == Domain::EMPTY){
            Integrator__v__grad_S::assemble_element_vector(ball_field_Si_r, e_data, vec);   
        }       
    });

    Logger::info("[T_Omega] - assemble RHS vector (complex part in conductor) [Sc_r].");
    assemble_vec(fe_system_.assemble_vec_data(Sc_c_), [&](auto& e_data, auto& vec) {
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR){
            Integrator__v__V::assemble_element_vector(ball_field_Sc_c_neg, e_data, vec);  

        }else if(property_id == Domain::CONDUCTOR_OUTER_LAYER){
            Integrator__v__V::assemble_element_vector(ball_field_Sc_c_neg, e_data, vec);   
        }       
    });


    Logger::info("[T_Omega] - assemble RHS vector (complex part in omega) [Si_r].");
    assemble_vec(fe_system_.assemble_vec_data(Si_c_), [&](auto& e_data, auto& vec) {
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR_OUTER_LAYER){
            Integrator__v__grad_S::assemble_element_vector(ball_field_Sc_c, e_data, vec);   //  -Hs = T - ∇Ω 

        }else if(property_id == Domain::EMPTY){
            Integrator__v__grad_S::assemble_element_vector(ball_field_Si_c, e_data, vec);   //   -Hs = -∇Ω
        }       
    });

 

    // scale with frequency  (Mc_r_ and Mc_c_ share the same matrix, same for other _r and _c)
    la_kernel::scale_mat(omega_, Mc_r_.mat);
    la_kernel::scale_mat(omega_, X__r_.mat);
    la_kernel::scale_mat(omega_, Ki_r_.mat);
    la_kernel::scale_mat(omega_, Xt_r_.mat);

    la_kernel::scale_vec(omega_, Sc_r_.vec);
    la_kernel::scale_vec(omega_, Sc_c_.vec);
    la_kernel::scale_vec(omega_, Si_r_.vec);
    la_kernel::scale_vec(omega_, Si_c_.vec);


    Logger::info("[T_Omega] - build linear system.");
    br_system_.build_linear_system();

    Logger::info("[T_Omega] - apply Dirichlet BC to the linear system.");
    bc_O_o_r_.apply_to_system(br_system_);
    bc_O_i_r_.apply_to_system(br_system_);
    bc_T_1_r_.apply_to_system(br_system_);
    bc_O_o_c_.apply_to_system(br_system_);
    bc_O_i_c_.apply_to_system(br_system_);
    bc_T_1_c_.apply_to_system(br_system_);


    return true;
}


bool T_Omega_c::assemble_pc()
{
    Logger::info("[T_Omega - preconditioner] - assemble edge interpolation matrix [P].");
    assemble_mat(fe_system_.assemble_mat_data(pc_P_), [&](auto& e_data, auto& mat) {
        Interpolator__H1_to_Hcurl::interpolate_element(e_data, mat);
    });


    // K + wM => 
    Logger::info("[T_Omega - preconditioner] - assemble preconditioner in (H1)^3 [LpM].");
    assemble_mat(fe_system_.assemble_mat_data(pc_LpM_), [&](auto& e_data, auto& mat) {
        double mu = 0;
        double sigma = 0;
        size_t property_id = e_data.e->get_property_id();
        if(property_id == Domain::CONDUCTOR) { mu = mu_conductor_; sigma = sigma_conductor_;}
        else if(property_id == Domain::CONDUCTOR_OUTER_LAYER) { mu = mu_conductor_; sigma = sigma_conductor_;}

        Integrator__s_grad_V__grad_V::assemble_element_matrix(1/sigma, e_data, mat);
        Integrator_H1__s_V__V::assemble_element_matrix(omega_*mu, e_data, mat);
    });


    Logger::info("[T_Omega - preconditioner] - assemble discrete gradient matrix [G].");
    assemble_mat(fe_system_.assemble_mat_data(pc_G_), [&](auto& e_data, auto& mat) {
        Interpolator__grad_H1_to_Hcurl::interpolate_element(e_data, mat);
    });

    Logger::info("[T_Omega - preconditioner] - assemble dof mapping from global H1 to omega H1 field [I].");
    assemble_mat(fe_system_.assemble_mat_data(pc_I_), [&](auto& e_data, auto& mat) {

        Identity_Mapping::direct_mapping(e_data, mat);
    });


    if(pc_alg_==1){
        Logger::info("[T_Omega - preconditioner] - assemble preconditioner [Qc].");
        assemble_mat(fe_system_.assemble_mat_data(pc_Qc_), [&](auto& e_data, auto& mat) {
            double mu = 0;
            size_t property_id = e_data.e->get_property_id();
            if(property_id == Domain::CONDUCTOR) mu = mu_conductor_; 
            else if(property_id == Domain::CONDUCTOR_OUTER_LAYER) mu = mu_conductor_; 
            Integrator__s_grad_S__grad_S::assemble_element_matrix(omega_*mu, e_data, mat);
        });


        Logger::info("[T_Omega - preconditioner] - assemble preconditioner [Qi].");
        assemble_mat(fe_system_.assemble_mat_data(pc_Qi_), [&](auto& e_data, auto& mat) {
            double mu = 0;
            size_t property_id = e_data.e->get_property_id();
            if(property_id == Domain::CONDUCTOR_OUTER_LAYER) mu = mu_conductor_; 
            else if(property_id == Domain::EMPTY) mu = mu_insulator_;
            Integrator__s_grad_S__grad_S::assemble_element_matrix(mu, e_data, mat);
        });
    }else if(pc_alg_==2){

        Mat M1, M2, M3, M4;
        Mat tM, oM;
        PetscCall(MatPtAP(Mc_c_.mat,pc_G_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M1));
        PetscCall(MatPtAP(Ki_c_.mat,pc_I_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M2));

        PetscCall(MatTransposeMatMult(pc_G_.mat,X__c_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&tM));
        PetscCall(MatTransposeMatMult(pc_I_.mat,Xt_c_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&oM));
        PetscCall(MatMatMult(tM, pc_I_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M3));
        PetscCall(MatMatMult(oM, pc_G_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M4));

        PetscCall(MatDuplicate(M1, MAT_COPY_VALUES,&pc_H_.mat));
        PetscCall(MatAXPY(pc_H_.mat, 1.0, M2, DIFFERENT_NONZERO_PATTERN));
        PetscCall(MatAXPY(pc_H_.mat, 1.0, M3, DIFFERENT_NONZERO_PATTERN));
        PetscCall(MatAXPY(pc_H_.mat, 1.0, M4, DIFFERENT_NONZERO_PATTERN));

        PetscCall(MatDestroy(&M1));
        PetscCall(MatDestroy(&M2));
        PetscCall(MatDestroy(&M3));
        PetscCall(MatDestroy(&M4));
        PetscCall(MatDestroy(&tM));
        PetscCall(MatDestroy(&oM));

    }else if(pc_alg_==3){
        Mat M1, M2, M3, M4;
        Mat tM, oM;
        PetscCall(MatPtAP(Mc_c_.mat,pc_G_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M1));
        PetscCall(MatPtAP(Ki_c_.mat,pc_I_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M2));

        PetscCall(MatTransposeMatMult(pc_G_.mat,X__c_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&tM));
        PetscCall(MatTransposeMatMult(pc_I_.mat,Xt_c_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&oM));
        PetscCall(MatMatMult(tM, pc_I_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M3));
        PetscCall(MatMatMult(oM, pc_G_.mat,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&M4));

        std::vector<Mat> blocks = {M1, M3, M4, M2};

        la_kernel::create_nest_mat(2, 2, blocks, pc_J_.mat);

        PetscCall(MatDestroy(&M1));
        PetscCall(MatDestroy(&M2));
        PetscCall(MatDestroy(&M3));
        PetscCall(MatDestroy(&M4));
        PetscCall(MatDestroy(&tM));
        PetscCall(MatDestroy(&oM));
    }



    return true;
}



bool T_Omega_c::solve_system()
{
    const G_Matrix lhs = br_system_.get_lhs();
    const G_Vector rhs = br_system_.get_rhs();
    const G_Vector x   = br_system_.get_x();

    petsc_util::petsc_save_ascii_mat(lhs, "lhs_Ac.txt");
    petsc_util::petsc_save_ascii_vec(rhs, "rhs_Ac.txt");


    bool successful_flag = false;
#ifdef LOAD_PETSC
    if(pc_alg_==0){
        // direct solver
        // Krilov subspace solver
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);

        // Set the operator (lhs matrix)
        KSPSetOperators(ksp, lhs, lhs);

        // Set solver type
        //KSPSetType(ksp, KSPMINRES);
        KSPSetType(ksp, KSPGMRES);
        
        //PetscCall(KSPGMRESSetRestart(ksp, 10000)); 

        // Configure the preconditioner (e.g., Jacobi)
        PC pc;
        KSPGetPC(ksp, &pc); 
        PCSetType(pc, PCILU);  
        PCFactorSetLevels(pc, 3); 
        //PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO); 
        //PCSetType(pc, PCHYPRE);
        //PCHYPRESetType(pc, "boomeramg");

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


        // Clean up
        KSPDestroy(&ksp);
    
    }else{


    }
#else
    Logger::error("[T_Omega] - this solver require petsc support!");
#endif

    return true;
}


bool T_Omega_c::solve_pc_system()
{

    return true;
}


scalar_t T_Omega_c::compute_L2_error()
{

    const double pi = CONST::PI;


    // manufactured solution
    // 3D vector field.  u = T - ∇Ω

    // ar = [1, 2, -1]   ai = [-2, 1, 1]
    V_Field_function T_r(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        double rr = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
        double phi = (0.25-rr)*(0.25-rr);

        v(0) = phi;
        v(1) = 2*phi;
        v(2) = -phi;

    });

    V_Field_function T_c(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        double rr = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
        double phi = (0.25-rr)*(0.25-rr);

        v(0) = -2*phi;
        v(1) = phi;
        v(2) = phi;

    });

    V_Field_function grad_O_r(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        v(0) = -2*x(0);
        v(1) = -2*x(1);
        v(2) = -2*x(2);

    });

    V_Field_function grad_O_c(mesh_, [&](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        
        double rr = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);

        v(0) = 9./4. - rr - 2*x(0)*x(0);
        v(1) = -2*x(0)*x(1);
        v(2) = -2*x(0)*x(2);

    });



    std::string debug_file_str = (pc_alg_!=0) ? "/pc_c_debug.txt" : "/c_debug.txt";

    std::string filepath = std::string(DEBUG_DATA_OUTPUT_DIR) + debug_file_str;
    std::ofstream file(filepath, std::ios::trunc);  
    

    Logger::info("[T_Omega] - compute L2 error.");
    scalar_t l2_error = integrate_element(br_system_, fe_system_, [&](Element_Data<3, 3>& e_data, scalar_t& result) {
        scalar_t local_integral = 0.;

        //if(e_data.e->get_property_id()!=5) return;
        //if(e_data.e->get_property_id()==1) return;

        const std::vector<const FEM_Space*>& space_list = *e_data.space_list;
        const std::vector<std::vector<scalar_t>>& dof_value_list = *e_data.dof_value_list; 
        const std::vector<size_t>& block_id_list = *e_data.block_id_list; 


        // compute T - grad Omega
        const std::vector<Integration_Point>& i_p_list = Integration::get_rule(e_data.b_shape, 3);

        Matrix<4, 3> H1_grad_basis;
        Matrix<4, 3> H1_phy_grad_basis;

        Matrix<6, 3> Hcurl_basis;
        Matrix<6, 3> Hcurl_phy_basis;

        Matrix<6, 6> dof_transform;

        Vector<3> solved_field_r = Vector<3>::Zero();
        Vector<3> solution_field_r = Vector<3>::Zero();

        Vector<3> solved_field_c = Vector<3>::Zero();
        Vector<3> solution_field_c = Vector<3>::Zero();

        Vector<3> last_solve_f_r; // for test
        Vector<3> last_exact_f_r; // for test
        Vector<3> last_solve_f_c; // for test
        Vector<3> last_exact_f_c; // for test

        Vector<3> temp_r;

        Vector<3> temp_c;

        for(const Integration_Point& i_p : i_p_list)
        {
            solved_field_r.setZero();
            solved_field_c.setZero();
            solution_field_r.setZero();
            solution_field_c.setZero();
            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));
            const Matrix<3, 3>& J_inv = e_data.get_inv_J(i_p.coord); 

            for(int i=0; i<space_list.size(); ++i)
            {
                int block_id = block_id_list[i];
                const FEM_Space* space = space_list[i];
                const std::vector<scalar_t>& dof_value = dof_value_list[i];
                if(block_id == Mc_r_.id){
                    // real part of T field
                    space->dof_transformation(e_data.e->get_node_idx(), dof_transform);
                    space->get_basis_v(i_p.coord, Hcurl_basis);
                    Hcurl_phy_basis = dof_transform*Hcurl_basis * J_inv;
                    for (int j = 0; j < dof_value.size(); ++j) {
                        solved_field_r += dof_value[j] * Hcurl_phy_basis.row(j).transpose();
                    }
                }

                if(block_id == Ki_r_.id){
                    // real part of Omega field
                    space->get_ED_basis_v(i_p.coord, H1_grad_basis);
                    H1_phy_grad_basis = H1_grad_basis * J_inv;
                    for (int j = 0; j < dof_value.size(); ++j) {
                        solved_field_r -= dof_value[j] * H1_phy_grad_basis.row(j).transpose();
                    }
                }

                if(block_id == Mc_c_.id){
                    // complex part of T field
                    space->dof_transformation(e_data.e->get_node_idx(), dof_transform);
                    space->get_basis_v(i_p.coord, Hcurl_basis);
                    Hcurl_phy_basis = dof_transform*Hcurl_basis * J_inv;
                    for (int j = 0; j < dof_value.size(); ++j) {
                        solved_field_c += dof_value[j] * Hcurl_phy_basis.row(j).transpose();
                    }
                }

                if(block_id == Ki_c_.id){
                    // complex part of Omega field
                    space->get_ED_basis_v(i_p.coord, H1_grad_basis);
                    H1_phy_grad_basis = H1_grad_basis * J_inv;
                    for (int j = 0; j < dof_value.size(); ++j) {
                        solved_field_c -= dof_value[j] * H1_phy_grad_basis.row(j).transpose();
                    }
                }

            }

            size_t property_id = e_data.e->get_property_id();
            if(property_id == Domain::CONDUCTOR){
                T_r.eval(i_p.coord, *e_data.e, temp_r);
                T_c.eval(i_p.coord, *e_data.e, temp_c);
                solution_field_r += temp_r;
                solution_field_c += temp_c;
                grad_O_r.eval(i_p.coord, *e_data.e, temp_r);
                grad_O_c.eval(i_p.coord, *e_data.e, temp_c);
                solution_field_r -= temp_r;
                solution_field_c -= temp_c;
            }else if(property_id == Domain::CONDUCTOR_OUTER_LAYER){
                T_r.eval(i_p.coord, *e_data.e, temp_r);
                T_c.eval(i_p.coord, *e_data.e, temp_c);
                solution_field_r += temp_r;
                solution_field_c += temp_c;
                grad_O_r.eval(i_p.coord, *e_data.e, temp_r);
                grad_O_c.eval(i_p.coord, *e_data.e, temp_c);
                solution_field_r -= temp_r;
                solution_field_c -= temp_c;
            }else if(property_id == Domain::EMPTY){
                grad_O_r.eval(i_p.coord, *e_data.e, temp_r);
                grad_O_c.eval(i_p.coord, *e_data.e, temp_c);
                solution_field_r -= temp_r;
                solution_field_c -= temp_c;
            }   

            last_solve_f_r = solved_field_r;
            last_exact_f_r = solution_field_r;
            last_solve_f_c = solved_field_c;
            last_exact_f_c = solution_field_c;

            double diff_sq = (solved_field_r - solution_field_r).squaredNorm() + (solved_field_c - solution_field_c).squaredNorm();
            local_integral += i_p.weight * abs_det_J * diff_sq;

        }

        result += local_integral;

        Vector<3> node_phys = e_data.physical_point(i_p_list.back().coord);
        size_t property_id = e_data.e->get_property_id();
        file <<"e id: "<<e_data.e->get_id()<<";   p id: "<<property_id<<std::endl;
        file << "phy coord: " << node_phys(0) << ", " << node_phys(1) << ", " << node_phys(2) << ", " << std::endl;
        file << "my result_r: " <<last_solve_f_r.transpose()  << std::endl;
        file << "solution_r: "  <<last_exact_f_r.transpose()<< std::endl;
        file << "my result_c: " <<last_solve_f_c.transpose()  << std::endl;
        file << "solution_c: "  <<last_exact_f_c.transpose()<< std::endl;
        file << "==============================" << std::endl;

    });

    file.close();

    return l2_error;

}


void T_Omega_c::finalize()
{
    la_kernel::destroy_mat(Kc_r_.mat);  
    la_kernel::destroy_mat(Mc_r_.mat);  
    la_kernel::destroy_mat(Ki_r_.mat);  
    la_kernel::destroy_mat(X__r_.mat);  
    la_kernel::destroy_mat(Xt_r_.mat);  

    la_kernel::destroy_mat(Kc_c_.mat);  
    Mc_c_.mat = nullptr;
    Ki_c_.mat = nullptr;
    X__c_.mat = nullptr;
    Xt_c_.mat = nullptr;

    la_kernel::destroy_vec(Sc_r_.vec);
    la_kernel::destroy_vec(Si_r_.vec);
    la_kernel::destroy_vec(Sc_c_.vec);  
    la_kernel::destroy_vec(Si_c_.vec);  

    la_kernel::destroy_mat(tmp_global_.mat);  
    la_kernel::destroy_mat(pc_Qc_.mat);  
    la_kernel::destroy_mat(pc_Qi_.mat);  
    la_kernel::destroy_mat(pc_H_.mat);  
    la_kernel::destroy_mat(pc_J_.mat); 

}



T_Omega_c::~T_Omega_c() 
{


}



struct Option 
{
    std::string_view long_name; 
    char             short_name; 
    bool             takes_value;
    std::function<void(std::string_view)> handler;
};


int main(int argc, char** argv) 
{
    std::string mesh_file    = "test_bball_2.geo";
    bool h_refine = false;
    bool l2_test = false;
    int pc_alg = 0;  // 0: no pc,  1: decouped_pc,  2: global_pc,  3: coupled_pc


    const std::vector<Option> options = {
        {"mesh",     'm',   true,  [&](std::string_view v) { mesh_file = std::string(v);          }},
        {"h-refine", 'h',   false, [&](std::string_view)   { h_refine = true;                     }},
        {"l2-test",  '\0',  false, [&](std::string_view)   { l2_test = true;                      }},
        {"pc",       '\0',  true,  [&](std::string_view v) { pc_alg = std::stoi(std::string(v));  }}
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



    if(h_refine)
    {
        const std::vector<std::pair<std::string, double>> mesh_sweep = {
            //{"test_cc_0.geo", 0.500000000},

            {"test_bball_2.geo", 0.250000000},
            {"test_bball_3.geo", 0.176776695},
            //{"test_bball_4.geo", 0.125000000},
            //{"test_cc_5.geo", 0.088388348},
            //{"test_cc_6.geo", 0.062500000},
            //{"test_cc_7.geo", 0.044200000},
            //{"test_cc_8.geo", 0.031250000},


        };

        std::string file_str = (pc_alg!=0) ? "/T_Omega_c_pcAlg"+std::to_string(pc_alg)+".dat" : "/T_Omega_c.dat";

        const std::string dat_path = DATA_OUTPUT_DIR + file_str;
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
                T_Omega_c T_O(mesh, pc_alg);
                T_O.initialize_system();
                if(pc_alg!=0) T_O.initialize_pc();
                Logger::stop_timer("Initialize T-Omega solver");

                Logger::start_timer("Assemble T-Omega matrix system");
                T_O.assemble_system();
                assemble_time += Logger::stop_timer("Assemble T-Omega matrix system");

                if(pc_alg!=0){
                    Logger::start_timer("Assemble T-Omega preconditioner");
                    T_O.assemble_pc();
                    assemble_time += Logger::stop_timer("Assemble T-Omega preconditioner");
                }

                for(int i=0; i<1; ++i)
                {

                    Logger::start_timer("Solve T-Omega matrix system");
                    if(pc_alg!=0){
                        T_O.solve_pc_system();
                    }else{
                        T_O.solve_system();
                    }
                    solving_time = Logger::stop_timer("Solve T-Omega matrix system");

                    n_iteration = T_O.get_n_iteration();
                    condition_numbder = T_O.get_system_condition();

                    if(l2_test){
                        Logger::start_timer("Compute L2 error.");
                        l2_error = T_O.compute_L2_error();
                        Logger::stop_timer("Compute L2 error.");
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
                    l2_convergence.flush();  

                }

                
                

                T_O.finalize();
            }

            

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
        T_Omega_c T_O(mesh, pc_alg);
        T_O.initialize_system();
        Logger::stop_timer("Initialize T-Omega solver");

        Logger::start_timer("Assemble T-Omega matrix system");
        T_O.assemble_system();
        Logger::stop_timer("Assemble T-Omega matrix system");

        if(pc_alg!=0){
            Logger::start_timer("Assemble T-Omega preconditioner");
            T_O.assemble_pc();
            Logger::stop_timer("Assemble T-Omega preconditioner");
        }

        Logger::start_timer("Solve T-Omega matrix system");
        if(pc_alg!=0){
            T_O.solve_pc_system();
        }else{
            T_O.solve_system();
        }
        Logger::stop_timer("Solve T-Omega matrix system");

        if(l2_test){
            Logger::start_timer("Compute L2 error.");
            scalar_t l2_error = T_O.compute_L2_error();
            Logger::stop_timer("Compute L2 error.");

            std::ostringstream ss;
            ss << std::scientific << std::setprecision(15) << l2_error;
            Logger::info("[T-Omega] - test case L2 error: " + ss.str());
        }

        

        T_O.finalize();
    }

    //PetscLogView(PETSC_VIEWER_STDOUT_WORLD);

    la_kernel::finalize();
    
    Logger::export_to_file("simuEM.log");
    return 0;
}