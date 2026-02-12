#pragma once

#include "mesh.h"
#include "config.h"
#include "utils/logger.h"
#include "utils/string_utils.h"

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

    size_t total_node = 0;
    size_t total_edge = 0;
    size_t total_triangle = 0;
    size_t total_quadrilateral = 0;
    size_t total_tetrahedron = 0;
    size_t total_hexahedron = 0;
    size_t total_prism = 0;
    size_t total_pyramid = 0;

    size_t counter_node = 0;
    size_t counter_edge = 0;
    size_t counter_triangle = 0;
    size_t counter_quadrilateral = 0;
    size_t counter_tetrahedron = 0;
    size_t counter_hexahedron = 0;
    size_t counter_prism = 0;
    size_t counter_pyramid = 0;

    Mesh_Format format_;

    std::tuple<Type, int, int> convert_Type(int gmsh_type);

    void create_mesh_element(Element *& element_pointer, Mesh& mesh, Type element_type, std::vector<std::size_t> node_idx, size_t element_id, size_t property_id=0, int o=1);
    Element * create_element(Type element_type, std::vector<std::size_t> node_idx, size_t element_id, size_t property_id=0, int o=1);


    void count_element_gmsh();


    // initialize vector size of mesh
    void initialize_mesh(Mesh& mesh);

public:
    Mesh_Parser(Mesh_Format format);

    Mesh load_mesh(const std::string& filename);

    Mesh load_gmsh(const std::string& filename);
    


};