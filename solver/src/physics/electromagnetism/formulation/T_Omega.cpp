#include "physics/electromagnetism/formulation/T_Omega.h"

using namespace simu;

T_Omega::T_Omega(Mesh& mesh) : mesh_(mesh)
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


        Key new_key_Omega_field_inner_boundary =  mesh_.mark_new_elements(scalar_field_Omega_inner_boundary_filter(), dim_-1, new_key_conductor_interface_layer, "Omega field inner boundary | "+description);
        
        std::string new_description_Omega_field_inner_boundary = mesh_.get_group_description(new_key_Omega_field_inner_boundary);
        key_Omega_field_boundary.push_back(new_key_Omega_field_inner_boundary);
        Logger::info("T_Omega - Create Omega field inner boundary, marked as: " + new_description_Omega_field_inner_boundary + ", #element = " + std::to_string(mesh_.get_element_group(new_key_Omega_field_inner_boundary).size()));
    }
    
        
};


