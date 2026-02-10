#include "entity/mesh/mesh_parser.h"
#include "entity/mesh/e_tetrahedron.h"

#include <map>
#include <algorithm>
#include <filesystem>



Mesh_Parser::Mesh_Parser(Mesh_Format format):format_(format)
{

}

Mesh Mesh_Parser::load_mesh(const std::string& filename)
{
    std::string ext = std::filesystem::path(filename).extension().string();
    switch (format_) {
        case Mesh_Format::GMSH:
            if (ext == ".msh" || ext == ".geo")
            {
                return load_gmsh(filename);
            } else {
                throw std::runtime_error("Error: Gmsh does not support " + ext + " format for file: " + filename);
            }
        default:
            throw std::runtime_error("Mesh format " + ext + " is not supported yet.");;
    }
}


std::vector<size_t> Mesh_Parser::build_tag_to_index_map(const std::vector<size_t>& node_tags)
{
    size_t max_tag = 0;
    for (auto t : node_tags)
        if (t > max_tag) max_tag = t;

    std::vector<size_t> map(max_tag + 1, static_cast<size_t>(-1));
    for (size_t i = 0; i < node_tags.size(); ++i)
        map[node_tags[i]] = i;

    return map;
}

Element* Mesh_Parser::create_element(int gmsh_type,
                                    const size_t* node_indices,
                                    size_t property_id)
{
    switch (gmsh_type) {
        case 4:  return new Tetrahedron(node_indices, property_id);
        // case 2:  return new Triangle(node_indices, property_id);
        // case 1:  return new Edge(node_indices, property_id);
        // case 15: return new NodeElement(node_indices, property_id);
        default: return nullptr;
    }
}


Mesh Mesh_Parser::load_gmsh(const std::string& filename)
{
#ifdef LOAD_GMSH
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::open(filename);

    Mesh mesh;
    mesh.dim_ = gmsh::model::getDimension();

    // ---- Nodes ----
    std::vector<size_t> node_tags;
    std::vector<double> coords;
    std::vector<double> parametric; // unused

    //gmsh::model::mesh::generate(mesh.dim_);
    gmsh::model::mesh::getNodes(node_tags, coords, parametric);
    mesh.n_node = node_tags.size();

    //auto tag_to_idx = build_tag_to_index_map(node_tags);

    // Store node coordinates indexed by 0-based index.
    // Gmsh returns nodes in node_tags order, which may not be contiguous,
    // so we use the tag_to_idx mapping.
    mesh.nodes.resize(mesh.n_node);
    for (size_t i = 0; i < node_tags.size(); ++i) {
        //size_t idx = tag_to_idx[node_tags[i]];
        mesh.nodes[node_tags[i]-1].x = coords[3 * i];
        mesh.nodes[node_tags[i]-1].y = coords[3 * i + 1];
        mesh.nodes[node_tags[i]-1].z = coords[3 * i + 2];
    }

    // ---- Build entity -> physical tag mapping ----
    // For each entity (dim, entity_tag), find which physical groups it belongs to.
    // An entity may have no physical group or one physical group.
    // We store only (dim, phys_tag) per entity.
    gmsh::vectorpair phys_dim_tags;
    gmsh::model::getPhysicalGroups(phys_dim_tags);



    gmsh::finalize();
    //std::cout<<mesh.n_element<<std::endl;
    return mesh;
#else
    throw std::runtime_error("Gmsh support not compiled!");
#endif
}