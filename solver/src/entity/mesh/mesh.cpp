#include "entity/mesh/mesh.h"



const std::vector<Element *>& Mesh::get_group(Key mesh_key) const
{
    auto it = element_group.find(mesh_key);
    if (it != element_group.end())
        return it->second;
    Logger::error("Mesh::get_group - failed: key not found in mesh, return empty vector<Element *>.");
    static const std::vector<Element*> empty;
    return empty;
}


const std::string Mesh::get_group_description(Key mesh_key) const
{ 
    auto it = element_group_description.find(mesh_key);
    if (it != element_group_description.end())
        return it->second;
    Logger::error("Mesh::get_group_description - failed: key not found in mesh, return empty description.");
    return "";
}

