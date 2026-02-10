#pragma once

#include "mesh.h"
#include "config.h"
#include "utils/logger.h"

#ifdef LOAD_GMSH
  #include "gmsh.h"
#endif

#include <string>


enum class Mesh_Format
{
    GMSH
};

class Mesh_Parser
{
private:
    Mesh_Format format_;
    static std::vector<size_t> build_tag_to_index_map(const std::vector<size_t>& node_tags);

    static Element* create_element(int gmsh_type, const size_t* node_indices, size_t property_id);

public:
    Mesh_Parser(Mesh_Format format);

    Mesh load_mesh(const std::string& filename);

    Mesh load_gmsh(const std::string& filename);


};