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

    std::tuple<Type, int, int> convert_Type(int gmsh_type);

    static Element * create_element(Type element_type, std::vector<std::size_t> node_idx, size_t property_id=0, int o=1);


public:
    Mesh_Parser(Mesh_Format format);

    Mesh load_mesh(const std::string& filename);

    Mesh load_gmsh(const std::string& filename);


};