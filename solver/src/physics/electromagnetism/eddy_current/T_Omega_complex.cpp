#include "physics/electromagnetism/eddy_current/T_Omega_complex.h"

using namespace simu;


T_Omega_c::T_Omega_c(Mesh& mesh, bool enable_preconditioner) : mesh_(mesh), enable_pc_(enable_preconditioner), fe_system_(mesh)
{

    dim_ = mesh.get_mesh_dimension();

    O_space_    = H1_Space(mesh.get_mesh_dimension(), 1);
    T_space_    = Hcurl_Space(mesh.get_mesh_dimension(),1);

    T_space_H1_ = H1_Space(mesh.get_mesh_dimension(), 1, true, 0);
    //T_H1_s_     = H1_Space(mesh.get_mesh_dimension(), 1);





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
    Ki_ = fe_system_.register_FE_space(O_space_, key_Omega_);

    
    Logger::info("[T_Omega] - register Dirichlet BC to Omega-field.");
    bc_O_o_ = fe_system_.register_Dirichlet_BC(Ki_, key_true_boundary_, Dirichlet_Type::HOMOGENEOUS);
    bc_O_i_ = fe_system_.register_Dirichlet_BC(Ki_, key_Omega_inner_boundary_1_, Dirichlet_Type::HOMOGENEOUS);
    
    Logger::block_info(Ki_.id, Ki_.row_offset, Ki_.col_offset, Ki_.row_size, Ki_.col_size);


    Logger::info("[T_Omega] - Assign space Hcurl to T-field region.");
    Kc_ = fe_system_.register_FE_space(T_space_, key_conductor_1_);
    Mc_ = fe_system_.register_FE_space(T_space_, key_conductor_1_);
    
    Logger::info("[T_Omega] - register Dirichlet BC to T-field.");
    bc_T_1_ = fe_system_.register_Dirichlet_BC(Kc_, key_conductor_interface_1_, Dirichlet_Type::HOMOGENEOUS);

    Logger::block_info(Kc_.id, Kc_.row_offset, Kc_.col_offset, Kc_.row_size, Kc_.col_size);  
    Logger::block_info(Mc_.id, Mc_.row_offset, Mc_.col_offset, Mc_.row_size, Mc_.col_size);  


    Logger::info("[T_Omega] - initialize coupling between T-field and Omega-field.");
    X__ = fe_system_.register_FE_space_coupling(Ki_, Kc_, key_conductor_interface_layer_1_);
    Logger::block_info(X__.id, X__.row_offset, X__.col_offset, X__.row_size, X__.col_size);
    
    //Block dof_coupling = fe_system_.register_FE_space_coupling(dof_T_[i], dof_Omega, key_conductor_interface_layer[i]);



    Xt_ = fe_system_.transpose_block(X__);
    Logger::block_info(Xt_.id,  Xt_.row_offset,  Xt_.col_offset,  Xt_.row_size,  Xt_.col_size);

    
}
