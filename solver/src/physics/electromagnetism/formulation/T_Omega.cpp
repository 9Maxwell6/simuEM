#include "physics/electromagnetism/formulation/T_Omega.h"

using namespace simu;

T_Omega::T_Omega(Mesh& mesh) : mesh_(mesh), fe_system_(mesh), Omega_space_(mesh.get_mesh_dimension(), 1), T_space_(mesh.get_mesh_dimension(),1)
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
            mesh_.set_group_property_id(mesh_key, 1);
            Logger::debug("T_Omega - Found source current: " + description);
        }
        else if(util::a_contains_b(description, {{"Insulator"}, {"Insulating"}}))
        {
            key_insulator.push_back(mesh_key);
            mesh_.set_group_property_id(mesh_key, 2);
            Logger::debug("T_Omega - Found insulating region: " + description);
        }
        else if(util::a_contains_b(description, {{"Conductor"}, {"Conducting"}}))
        {
            key_conductor.push_back(mesh_key);
            mesh_.set_group_property_id(mesh_key, 3);
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




    Logger::info("[T_Omega] - Assign space H1 to Omega-field region");
    //Block dof_Omega = fe_system_.register_FE_space(field_Omega, new_key_Omega_field);
    dof_Omega_ = fe_system_.register_FE_space(Omega_space_, new_key_Omega_field);
    Logger::block_info(dof_Omega_.id, dof_Omega_.row_offset, dof_Omega_.col_offset, dof_Omega_.row_size, dof_Omega_.col_size);



    Logger::info("[T_Omega] - Assign space Hcurl to T-field region");
    for(size_t i = 0; i < key_conductor.size(); ++i)
    {
        dof_T_.push_back(fe_system_.register_FE_space(T_space_, key_conductor[i]));
        Logger::block_info(dof_T_[i].id, 
                           dof_T_[i].row_offset, 
                           dof_T_[i].col_offset, 
                           dof_T_[i].row_size, 
                           dof_T_[i].col_size);  
    }
    //Block dof_T = fe_system_.register_FE_space(field_T, key_conductor[0]);
    
    
    



    Logger::info("[T_Omega] - initialize coupling between T-field and Omega-field");
    for(size_t i = 0; i < key_conductor_interface_layer.size(); ++i)
    {
        dof_coupling_.push_back(fe_system_.register_FE_space_coupling(dof_Omega_, dof_T_[i], key_conductor_interface_layer[i]));
        Logger::block_info(dof_coupling_[i].id, 
                           dof_coupling_[i].row_offset, 
                           dof_coupling_[i].col_offset, 
                           dof_coupling_[i].row_size, 
                           dof_coupling_[i].col_size);
    }
    //Block dof_coupling = fe_system_.register_FE_space_coupling(dof_T_[i], dof_Omega, key_conductor_interface_layer[i]);
    



    for(size_t i = 0; i < dof_coupling_.size(); ++i)
    {
        dof_coupling_tp_.push_back(fe_system_.transpose_block(dof_coupling_[i]));

        Logger::block_info(dof_coupling_tp_[i].id, 
                           dof_coupling_tp_[i].row_offset, 
                           dof_coupling_tp_[i].col_offset, 
                           dof_coupling_tp_[i].row_size, 
                           dof_coupling_tp_[i].col_size);
    }
    //Block dof_coupling_tp = fe_system_.transpose_block(dof_coupling);


    Logger::info("[T_Omega] - delete temporary block hash.");
    fe_system_.delete_block_hash();

    Block_Rack br_l = fe_system_.initialize_block_rack(2, 2);
    
    br_l.insert_block(dof_Omega_,          0, 0);
    br_l.insert_block(dof_T_[0],           1, 1);
    br_l.insert_block(dof_coupling_[0],    0, 1);
    br_l.insert_block(dof_coupling_tp_[0], 1, 0);
    br_l.compute_block_offset();
    Logger::info("[T_Omega] - initialize block rack: \n"+br_l.print_block_rack());

    
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
    assemble_block(fe_system_.assemble_data(dof_Omega_), [&](auto& e_data, auto& mat) {
        double sigma = 0.;
        if(e_data.e->get_property_id()==3) sigma = 1.;
        //std::cout<<e_data.e->get_property_id()<<std::endl;

        Integrator__s_S__S::assemble_element_matrix(sigma, e_data, mat);

    });

    /*
    Logger::info("[T_Omega] - assemble T-field block matrix.");
    for(Block& dof_T : dof_T_)
    {
        assemble_block(fe_system_.assemble_data(dof_T), [&](auto& e_data, auto& mat) {
            double sigma = 0.;
            if(e_data.e->get_property_id()==1) sigma = 1.;

            Integrator__s_S__S::assemble_element_matrix(sigma, e_data, mat);

        });

    }
    */


    

}


