#include "entity/mesh/mesh_parser.h"

#include <map>
#include <algorithm>

#include "gmsh.h"
#include "entity/mesh/e_tetrahedron.h"
// #include "Triangle.h"
// #include "Edge.h"
// #include "NodeElement.h"

std::vector<size_t> Mesh_Parser::build_tag_to_index_map(
    const std::vector<size_t>& node_tags)
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


Mesh Mesh_Parser::load_mesh(const std::string& filename) {
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::open(filename);

    Mesh mesh;
    mesh.dim_ = gmsh::model::getDimension();

    // ---- Nodes ----
    std::vector<size_t> node_tags;
    std::vector<double> coords;
    std::vector<double> parametric; // unused
    gmsh::model::mesh::getNodes(node_tags, coords, parametric);
    mesh.n_node = node_tags.size();
    // TODO: store coords into mesh node storage

    auto tag_to_idx = build_tag_to_index_map(node_tags);

    // Store node coordinates indexed by 0-based index.
    // Gmsh returns nodes in node_tags order, which may not be contiguous,
    // so we use the tag_to_idx mapping.
    mesh.nodes.resize(mesh.n_node);
    for (size_t i = 0; i < node_tags.size(); ++i) {
        size_t idx = tag_to_idx[node_tags[i]];
        mesh.nodes[idx].x = coords[3 * i];
        mesh.nodes[idx].y = coords[3 * i + 1];
        mesh.nodes[idx].z = coords[3 * i + 2];
    }

    // ---- Build entity -> physical tag mapping ----
    // For each entity (dim, entity_tag), find which physical groups it belongs to.
    // An entity may have no physical group or one physical group.
    // We store only (dim, phys_tag) per entity.
    gmsh::vectorpair phys_dim_tags;
    gmsh::model::getPhysicalGroups(phys_dim_tags);

    // entity (dim, entity_tag) -> physical group tag
    std::map<std::pair<int,int>, int> entity_to_phys;

    // physical (dim, phys_tag) -> name
    std::map<std::pair<int,int>, std::string> phys_names;

    for (const auto& [dim, phys_tag] : phys_dim_tags) {
        std::string name;
        gmsh::model::getPhysicalName(dim, phys_tag, name);
        phys_names[{dim, phys_tag}] = name;

        std::vector<int> entity_tags;
        gmsh::model::getEntitiesForPhysicalGroup(dim, phys_tag, entity_tags);
        for (int etag : entity_tags) {
            entity_to_phys[{dim, etag}] = phys_tag;
        }
    }

    // ---- Map each (dim, phys_tag) to a Key ----
    // When we encounter a new physical group, increment dim_keys[dim].id
    std::map<std::pair<int,int>, Key> phys_to_key;

    // ---- Load elements entity-by-entity ----
    gmsh::vectorpair all_entities;
    gmsh::model::getEntities(all_entities);

    for (const auto& [ent_dim, ent_tag] : all_entities) {
        std::vector<int> elem_types;
        std::vector<std::vector<size_t>> elem_tags_vec;
        std::vector<std::vector<size_t>> node_tags_vec;

        gmsh::model::mesh::getElements(
            elem_types, elem_tags_vec, node_tags_vec, ent_dim, ent_tag);

        // Look up physical group for this entity
        int phys_tag = 0;
        auto it = entity_to_phys.find({ent_dim, ent_tag});
        if (it != entity_to_phys.end())
            phys_tag = it->second;

        // Get or create the Key for this physical group
        Key key{ent_dim, 0}; // key.id = 0 means no physical group
        if (phys_tag != 0) {
            auto kt = phys_to_key.find({ent_dim, phys_tag});
            if (kt != phys_to_key.end()) {
                key = kt->second;
            } else {
                // New physical group: allocate a key
                mesh.dim_keys[ent_dim].id++;
                key = mesh.dim_keys[ent_dim];
                phys_to_key[{ent_dim, phys_tag}] = key;

                // Initialize group
                mesh.element_group[key] = {};
                mesh.group_description[key] = phys_names[{ent_dim, phys_tag}];
            }
        }

        // Process elements
        for (size_t t = 0; t < elem_types.size(); ++t) {
            int gmsh_type = elem_types[t];
            size_t num_elems = elem_tags_vec[t].size();
            if (num_elems == 0) continue;

            size_t nodes_per_elem = node_tags_vec[t].size() / num_elems;

            for (size_t i = 0; i < num_elems; ++i) {
                std::vector<size_t> idx(nodes_per_elem);
                for (size_t j = 0; j < nodes_per_elem; ++j) {
                    idx[j] = tag_to_idx[node_tags_vec[t][i * nodes_per_elem + j]];
                }

                size_t property_id = static_cast<size_t>(key.id);
                Element* elem = create_element(gmsh_type, idx.data(), property_id);
                if (!elem) continue;

                // Highest-dimension elements go into mesh.elements
                if (ent_dim == mesh.dim_) {
                    mesh.elements.push_back(elem);
                }

                // If it belongs to a physical group, add to element_group
                if (phys_tag != 0) {
                    mesh.element_group[key].push_back(elem);
                }
            }
        }
    }

    mesh.n_element = mesh.elements.size();

    gmsh::finalize();
    std::cout<<mesh.n_element<<std::endl;
    return mesh;
}