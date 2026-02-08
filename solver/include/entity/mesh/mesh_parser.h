#pragma once

#include <string>
#include "mesh.h"

class Mesh_Parser {
public:
    static Mesh load_mesh(const std::string& filename);

private:
    static std::vector<size_t> build_tag_to_index_map(
        const std::vector<size_t>& node_tags);

    static Element* create_element(int gmsh_type,
                                   const size_t* node_indices,
                                   size_t property_id);
};