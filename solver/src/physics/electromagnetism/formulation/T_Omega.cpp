#include "physics/electromagnetism/formulation/T_Omega.h"

using namespace simu;

T_Omega::T_Omega(Mesh& mesh) : mesh_(mesh), fe_system(mesh)
{

    dim_ = mesh.get_mesh_dimension();

    for(const Key& mesh_key : mesh_.get_keys_true_boundary())
    {
        std::string description = mesh.get_group_description(mesh_key);
        key_true_boundary.push_back(mesh_key);
        key_Omega_field_boundary.push_back(mesh_key);
        Logger::info("T_Omega - Found simulation boundary: " + description);
    }

    for(const Key& mesh_key : mesh_.get_keys_internal_surface())
    {
        std::string description = mesh.get_group_description(mesh_key);
        if(util::a_contains_b(description, {{"Conductor", "Boundary"}, {"Conducting", "Boundary"}})){
            key_conductor_interface.push_back(mesh_key);
            Logger::debug("T_Omega - Found conductor boundary: " + description);
        }
    }
    
    for(const Key& mesh_key : mesh_.get_keys_domain())
    {
        std::string description = mesh.get_group_description(mesh_key);
        if(util::a_contains_b(description, {{"Source, Current"}}))
        {
            key_source.push_back(mesh_key);
            Logger::debug("T_Omega - Found source current: " + description);
        }
        else if(util::a_contains_b(description, {{"Insulator"}, {"Insulating"}}))
        {
            key_insulator.push_back(mesh_key);
            Logger::debug("T_Omega - Found insulating region: " + description);
        }
        else if(util::a_contains_b(description, {{"Conductor"}, {"Conducting"}}))
        {
            key_conductor.push_back(mesh_key);
            Logger::debug("T_Omega - Found conductor: " + description);
        }
    }

    // mark all elements in conductors' boundary. (in 3D mesh, they are the 3D elements at boundary layer inside conductor)
    for(Key& key : key_conductor)
    {
        std::string description = mesh_.get_group_description(key);

        Key new_key_conductor_interface_layer =  mesh_.mark_elements(conductor_interface_layer_filter(), key, "Conductor Outerlayer | "+description);

        std::string new_description_conductor_interface_layer = mesh_.get_group_description(new_key_conductor_interface_layer);
        key_conductor_interface_layer.push_back(new_key_conductor_interface_layer);
        Logger::info("T_Omega - Create conductor interface layer, marked as: " + new_description_conductor_interface_layer + ", #element = " + std::to_string(mesh_.get_element_group(new_key_conductor_interface_layer).size()));
        Logger::mesh_entity(new_key_conductor_interface_layer.dim, -1, new_key_conductor_interface_layer.id, 
                                                                       mesh.get_element_group(new_key_conductor_interface_layer).size(), 
                                                                       mesh.get_element_geometry_size_group(new_key_conductor_interface_layer).size(),
                                                                       mesh.get_element_size_group(new_key_conductor_interface_layer),
                                                                       mesh.get_group_description(new_key_conductor_interface_layer));

        Key new_key_Omega_field_inner_boundary =  mesh_.mark_new_elements(scalar_field_Omega_inner_boundary_filter(), dim_-1, new_key_conductor_interface_layer, "Omega field inner boundary | "+description);
        
        std::string new_description_Omega_field_inner_boundary = mesh_.get_group_description(new_key_Omega_field_inner_boundary);
        key_Omega_field_boundary.push_back(new_key_Omega_field_inner_boundary);
        Logger::info("T_Omega - Create Omega-field inner boundary, marked as: " + new_description_Omega_field_inner_boundary + ", #element = " + std::to_string(mesh_.get_element_group(new_key_Omega_field_inner_boundary).size()));
        Logger::mesh_entity(new_key_Omega_field_inner_boundary.dim, -1, new_key_Omega_field_inner_boundary.id, 
                                                                        mesh.get_element_group(new_key_Omega_field_inner_boundary).size(), 
                                                                        mesh.get_element_geometry_size_group(new_key_Omega_field_inner_boundary).size(),
                                                                        mesh.get_element_size_group(new_key_Omega_field_inner_boundary),
                                                                        mesh.get_group_description(new_key_Omega_field_inner_boundary));
    }

    Logger::info("T_Omega - Create Omega-field group.");
    Key new_key_Omega_field = mesh_.group_union(key_insulator[0], key_conductor_interface_layer[0], "Omega field");
    
    Logger::mesh_entity(new_key_Omega_field.dim, -1, new_key_Omega_field.id, 
                                                     mesh.get_element_group(new_key_Omega_field).size(), 
                                                     mesh.get_element_geometry_size_group(new_key_Omega_field).size(),
                                                     mesh.get_element_size_group(new_key_Omega_field),
                                                     mesh.get_group_description(new_key_Omega_field));

    Logger::info("[T_Omega] - Create function space Hcurl");
    Hcurl_Space field_T(dim_,1);

    Logger::info("[T_Omega] - Assign space Hcurl to T-field region");
    Block dof_T = fe_system.register_FE_space(field_T, key_conductor[0]);
    Logger::block_info(dof_T.id, dof_T.row_offset, dof_T.col_offset, dof_T.row_size, dof_T.col_size);                                                 

    Logger::info("[T_Omega] - Create function space H1");
    H1_Space field_Omega(dim_,1);

    Logger::info("[T_Omega] - Assign space H1 to Omega-field region");
    Block dof_Omega = fe_system.register_FE_space(field_Omega, new_key_Omega_field);
    Logger::block_info(dof_Omega.id, dof_Omega.row_offset, dof_Omega.col_offset, dof_Omega.row_size, dof_Omega.col_size);


    Logger::info("[T_Omega] - initialize coupling between T-field and Omega-field");
    Block dof_coupling = fe_system.register_FE_space_coupling(dof_T, dof_Omega, key_conductor_interface_layer[0]);
    Logger::block_info(dof_coupling.id, dof_coupling.row_offset, dof_coupling.col_offset, dof_coupling.row_size, dof_coupling.col_size);


    Block dof_coupling_tp = fe_system.transpose_block(dof_coupling);
    Logger::block_info(dof_coupling_tp.id, dof_coupling_tp.row_offset, dof_coupling_tp.col_offset, dof_coupling_tp.row_size, dof_coupling_tp.col_size);


    Logger::info("[T_Omega] - delete temporary block hash.");
    fe_system.delete_block_hash();

    Block_Rack br_l = fe_system.initialize_block_rack(2, 2);
    br_l.insert_block(dof_T,           0, 0);
    br_l.insert_block(dof_Omega,       1, 1);
    br_l.insert_block(dof_coupling,    0, 1);
    br_l.insert_block(dof_coupling_tp, 1, 0);
    br_l.compute_block_offset();
    Logger::info("[T_Omega] - initialize block rack: \n"+br_l.print_block_rack());

    
    fe_system.assemble_data(dof_T);
};


