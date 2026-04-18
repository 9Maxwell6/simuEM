#include "physics/electromagnetism/formulation/T_Omega.h"

using namespace simu;

T_Omega::T_Omega(Mesh& mesh) : mesh_(mesh), fe_system_(mesh)
{

    dim_ = mesh.get_mesh_dimension();

    Omega_space_ = H1_Space(mesh.get_mesh_dimension(), 1);
    T_space_1_ = Hcurl_Space(mesh.get_mesh_dimension(),1);

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


        if(util::a_contains_b(description, {{"Source, Current"}}))
        {   
            key_source_ = mesh_key;
            mesh_.set_group_property_id(mesh_key, 1);
            Logger::debug("T_Omega - Found source current: " + description);
        }
        else if(util::a_contains_b(description, {{"Insulator"}, {"Insulating"}}))
        {
            key_insulator_=mesh_key;
            mesh_.set_group_property_id(mesh_key, 2);
            Logger::debug("T_Omega - Found insulating region: " + description);
        }
        else if(util::a_contains_b(description, {{"Conductor", "1"}, {"Conducting", "1"}}))
        {
            key_conductor_1_ = mesh_key;
            mesh_.set_group_property_id(mesh_key, 3);
            Logger::debug("T_Omega - Found conductor: " + description);
        }
    }

    // mark all elements in conductors' boundary. (in 3D mesh, they are the 3D elements at boundary layer inside conductor)
    
    std::string description = mesh_.get_group_description(key_conductor_1_);
    std::cout<<key_conductor_interface_1_.dim<<std::endl;
    std::cout<<key_conductor_interface_1_.id<<std::endl;
    key_conductor_interface_layer_1_ =  mesh_.mark_elements(conductor_interface_layer_filter(key_conductor_interface_1_), key_conductor_1_, "Conductor Outerlayer | "+description);

    std::string new_description_conductor_interface_layer = mesh_.get_group_description(key_conductor_interface_layer_1_);
    
    Logger::info("T_Omega - Create conductor interface layer, marked as: " + new_description_conductor_interface_layer + ", #element = " + std::to_string(mesh_.get_element_group(key_conductor_interface_layer_1_).size()));
    Logger::mesh_entity(key_conductor_interface_layer_1_.dim, -1, key_conductor_interface_layer_1_.id, 
                                                                    mesh.get_element_group(key_conductor_interface_layer_1_).size(), 
                                                                    mesh.get_element_geometry_size_group(key_conductor_interface_layer_1_).size(),
                                                                    mesh.get_element_size_group(key_conductor_interface_layer_1_),
                                                                    mesh.get_group_description(key_conductor_interface_layer_1_));


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
    //   i_r_list      - integration rules per order
    //
    //
    // Template parameters:
    //   phy_dim - physical space dimension (1, 2, or 3)
    //   ref_dim - reference element dimension (1, 2, or 3), ref_dim <= phy_dim
    //
    // Note: accessed via auto& in user callbacks. Use e_data.e, e_data.J, etc.

    
    Logger::info("[T_Omega] - assemble Omega-field block matrix.");
    assemble_mat(fe_system_.assemble_mat_data(dof_Omega_), [&](auto& e_data, auto& mat) {
        double sigma = 1.;
        if(e_data.e->get_property_id()==3) sigma = 1.;
        //std::cout<<e_data.e->get_property_id()<<std::endl;

        Integrator__s_S__S::assemble_element_matrix(sigma, e_data, mat);

        Integrator__s_grad_S__grad_S::assemble_element_matrix(sigma, e_data, mat);

    });

    //*
    Logger::info("[T_Omega] - assemble T-field block matrix.");
    assemble_mat(fe_system_.assemble_mat_data(dof_T_1_), [&](auto& e_data, auto& mat) {
        double sigma = 1.;
        if(e_data.e->get_property_id()==3) sigma = 1.;

        Integrator__s_curl_V__curl_V::assemble_element_matrix(sigma, e_data, mat);
        Integrator__s_V__V::assemble_element_matrix(sigma, e_data, mat);

    });

    
    //*/

    //*
    Logger::info("[T_Omega] - assemble coupling block matrix.");
    assemble_mat(fe_system_.assemble_mat_data(dof_coupling_1_), [&](auto& e_data, auto& mat) {
        double sigma = 1.;
        if(e_data.e->get_property_id()==3) sigma = 1.;

        
        //Integrator__s_curl_V__curl_V::assemble_element_matrix(sigma, e_data, mat);
        Integrator__s_grad_S__V::assemble_element_matrix(sigma, e_data, mat);

    });
    
    dof_coupling_tp_1_.block_transpose(dof_coupling_1_);
    
    //*/

    


    V_Field_function f([](Eigen::Ref<const VectorXd> x, const Field_Data& fd, Eigen::Ref<VectorXd> v) {
        v(0) = std::sin(x(0));
        v(1) = std::cos(x(1));
        v(2) = std::cos(x(2));
    });

    Logger::info("[T_Omega] - assemble RHS vector for Omega block.");
    assemble_vec(fe_system_.assemble_vec_data(dof_Omega_), [&](auto& e_data, auto& vec) {
        double mu = 1.;
        if(e_data.e->get_property_id()==3) mu = 1.;

        Integrator__v__grad_S::assemble_element_vector(f, e_data, vec);
        

    });



    Logger::info("[T_Omega] - assemble RHS vector for T-field block.");
    assemble_vec(fe_system_.assemble_vec_data(dof_T_1_), [&](auto& e_data, auto& vec) {


        Integrator__v__V::assemble_element_vector(f, e_data, vec);

    });

    


    Logger::info("[T_Omega] - build linear system.");
    br_system_.build_linear_system();

    Logger::info("[T_Omega] - apply Dirichlet BC to the linear system.");
    bc_Omega_out_.apply_to_system(br_system_);
    bc_Omega_in_.apply_to_system(br_system_);
    bc_T_1_.apply_to_system(br_system_);

    br_system_.solve();

    
    return true;
}


