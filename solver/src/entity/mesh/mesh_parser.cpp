#include "entity/mesh/mesh_parser.h"

#include <map>
#include <algorithm>
#include <filesystem>

using namespace simu;

Mesh_Parser::Mesh_Parser(Mesh_Format format):format_(format)
{

}



/**
 * @brief Converts a Gmsh element type ID into internal mesh properties.
 * * This function maps the standard Gmsh element identifiers (found in the MSH 
 * file format or returned by gmsh::model::mesh::getElements) to the internal 
 * Type enum, the polynomial order, and the total number of nodes per element.
 * 
 * 
 * @param gmsh_type The integer ID provided by the Gmsh API (e.g., 2 for a 3-node triangle).
 * @return std::tuple<Type, int, int> A tuple containing:
 * - [0]: Internal Type enum (e.g., Type::TRIANGLE).
 * - [1]: Polynomial order (1 for linear, 2 for quadratic, etc.).
 * - [2]: Number of nodes per element (used as a stride for the node_tags vector).
 * * @throws std::invalid_argument If the gmsh_type does not map to a supported internal shape.
 * * @note Current support is limited to 1st-order linear elements:
 * - Type 1: 2-node line
 * - Type 2: 3-node triangle
 * - Type 4: 4-node tetrahedron
 * - Type 15: 1-node point
 */
std::tuple<Type, int, int> Mesh_Parser::convert_Type(int gmsh_type) {
    switch (gmsh_type) {
        case 1:  return {Type::EDGE,        1, 2}; // 2-node line
        case 2:  return {Type::TRIANGLE,    1, 3}; // 3-node triangle
        case 4:  return {Type::TETRAHEDRON, 1, 4}; // 4-node tetrahedron
        case 15: return {Type::NODE,        1, 1}; // 1-node point
        default:
            // Handle unsupported shapes (Quads, Hexahedra, etc.)
            Logger::error("Mesh_Parser::convert_Type faliure - Unsupported Gmsh element type: "+ std::to_string(gmsh_type));
            throw std::invalid_argument("Unsupported Gmsh element type: " + std::to_string(gmsh_type));
    }
}



void Mesh_Parser::create_mesh_element(Element*& element_pointer, Mesh& mesh, Type element_type, std::vector<std::size_t> node_idx, size_t element_id, size_t property_id, int o)
{
    try {
        switch (element_type) 
        {
            case Type::EDGE:
            {
                mesh.elements_edge_.at(counter_edge) = Edge(node_idx, element_id, property_id, 1);
                element_pointer = &mesh.elements_edge_[counter_edge];
                counter_edge++;
                break;
            }
            case Type::TRIANGLE: 
            {    
                mesh.elements_triangle_.at(counter_triangle) = Triangle(node_idx, element_id, property_id, 1);
                element_pointer = &mesh.elements_triangle_[counter_triangle];
                counter_triangle++;
                break;
            }
            case Type::TETRAHEDRON:
            {
                mesh.elements_tetrahedron_.at(counter_tetrahedron) = Tetrahedron(node_idx, element_id, property_id, 1);
                element_pointer = &mesh.elements_tetrahedron_[counter_tetrahedron];
                counter_tetrahedron++;
                break;
            }
            default: Logger::warning("Mesh_Parser::create_mesh_element - failed: element_type not available");
        }
    } catch (const std::out_of_range& e) {
        Logger::error(std::string("Mesh_Parser::create_mesh_element - failed: ")+e.what());
    }
}

Element * Mesh_Parser::create_element(Type element_type, std::vector<std::size_t> node_idx, size_t element_id, size_t property_id, int o)
{
    switch (element_type) 
    {
        case Type::EDGE: return new Edge(node_idx, element_id, property_id, o);
        case Type::TRIANGLE: return new Triangle(node_idx, element_id, property_id, o);
        case Type::TETRAHEDRON: return new Tetrahedron(node_idx, element_id, property_id, o);
        default: 
            Logger::warning("Mesh_Parser::create_element - failed: element_type not available, return nullptr");
            return nullptr;
    }

}

void Mesh_Parser::count_element_gmsh()
{
    total_node = 0;
    total_edge = 0;
    total_triangle = 0;
    total_quadrilateral = 0;
    total_tetrahedron = 0;
    total_hexahedron = 0;
    total_prism = 0;
    total_pyramid = 0;

    for (int dim = 0; dim <= 3; dim++) 
    {   
        std::vector<int> element_types;
        std::vector<std::vector<size_t>> element_tags;
        std::vector<std::vector<size_t>> node_tags;
        gmsh::model::mesh::getElements(element_types, element_tags, node_tags, dim);

        for (size_t i = 0; i < element_types.size(); i++) 
        {
            std::string name;
            std::vector<double> ref_coords;
            int dim, order, n_node, n_vertex;
            gmsh::model::mesh::getElementProperties(element_types[i], name, dim, order, n_node, ref_coords, n_vertex);

            size_t count = element_tags[i].size();

            if (dim == 0) {
                total_node += count;
            } else if (dim == 1 && n_vertex == 2) {
                total_edge += count;
            } else if (dim == 2 && n_vertex == 3) {
                total_triangle += count;
            } else if (dim == 2 && n_vertex == 4) {
                total_quadrilateral += count;
            } else if (dim == 3 && n_vertex == 4) {
                total_tetrahedron += count;
            } else if (dim == 3 && n_vertex == 5) {
                total_pyramid += count;
            } else if (dim == 3 && n_vertex == 6) {
                total_prism += count;
            } else if (dim == 3 && n_vertex == 8) {
                total_hexahedron += count;
            }
        }
    }
}



void Mesh_Parser::initialize_mesh(Mesh& mesh)
{
    switch (mesh.dim_) 
    {
        case 0: mesh.nodes_.resize(total_node); 
        case 1: 
        {
            mesh.elements_edge_.resize(total_edge);
            mesh.elements_.reserve(total_edge);
            break;
        }
        case 2: 
        {
            mesh.elements_triangle_.resize(total_triangle);
            mesh.elements_.reserve(total_triangle);
            break;
        }
        case 3: 
        {
            mesh.elements_tetrahedron_.resize(total_tetrahedron);
            mesh.elements_.reserve(total_tetrahedron); // total_tetrahedron + total_hexahedron
            break;
        }
        default:
            Logger::error("Mesh_Parser::initialize_mesh - impossible mesh dimension: return bad mesh.");
    }
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
                Logger::error("Mesh_Parser::load_mesh - failed: Gmsh does not support " + ext + " format for file: " + filename);
                throw std::runtime_error("Error: Gmsh does not support " + ext + " format for file: " + filename);
            }
        default:
            Logger::error("Mesh format " + ext + " is not supported yet.");
            throw std::runtime_error("Mesh format " + ext + " is not supported yet.");
    }
}




Mesh Mesh_Parser::load_gmsh(const std::string& filename)
{
#ifdef LOAD_GMSH
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1);
    gmsh::open(filename);

    Logger::info("Load mesh: " + filename);

    Mesh mesh;
    mesh.dim_ = gmsh::model::getDimension();

    // initialize element counter
    count_element_gmsh();
    initialize_mesh(mesh);

    // ---- Nodes ----
    std::vector<size_t> node_tags;
    std::vector<double> coords;
    std::vector<double> parametric; // unused

    gmsh::model::mesh::generate(mesh.dim_);
    gmsh::model::mesh::getNodes(node_tags, coords, parametric);
    mesh.n_node = node_tags.size();


    // Store node coordinates indexed by 0-based index.
    // Gmsh returns nodes in node_tags order, which may not be contiguous,
    // so we use the tag_to_idx mapping.
    mesh.nodes_.resize(mesh.n_node);
    for (size_t i = 0; i < node_tags.size(); ++i) 
    {
        //size_t idx = tag_to_idx[node_tags[i]];
        mesh.nodes_[node_tags[i]-1].x = coords[3 * i];
        mesh.nodes_[node_tags[i]-1].y = coords[3 * i + 1];
        mesh.nodes_[node_tags[i]-1].z = coords[3 * i + 2];
    }



    // Store curve elements into mesh.curve
    

    // Store surface elements into mesh.surface  (only if mesh is 2D!)


    // Store volume elements into mesh.volume   (only if mesh is 3D!)


    // Store all elements with the same dimension of mesh
    std::unordered_map<size_t, Element*> element_map_md;
    
    std::vector<int> element_types_mesh_md;
    std::vector<std::vector<size_t>> element_tags_mesh_md;
    std::vector<std::vector<size_t>> node_tags_mesh_md;
    // element_types and element_tags should have same size
    gmsh::model::mesh::getElements(element_types_mesh_md, element_tags_mesh_md, node_tags_mesh_md, mesh.dim_, -1);
    //element_tags.size() is number of different shapes/order
    for (size_t i = 0; i < element_tags_mesh_md.size(); ++i) 
    {   
        auto [type, order, n_node] = convert_Type(element_types_mesh_md[i]);

        size_t n_elements = element_tags_mesh_md[i].size();
        for (size_t j=0; j<n_elements; ++j)
        {   
            size_t element_id = element_tags_mesh_md[i][j];
            auto startIt = node_tags_mesh_md[i].begin() + (j * n_node);
            auto endIt   = startIt + n_node;
            std::vector<size_t> element_node_tags(startIt, endIt);
            //Element * element = create_element(type, element_node_tags, 0, order);

            Element * element;
            create_mesh_element(element, mesh, type, element_node_tags, element_id, 0, order);
            
            mesh.elements_.push_back(element);
        }
    }

    for (Element* el : mesh.elements_) {
        element_map_md[el->get_Id()] = el; 
    }
    


    // Store volume elements into mesh.volume   (only if mesh is 3D!)
 
    


    // create key of element groups generated by gmsh file
    gmsh::vectorpair                   physicalGroups_dim_tag;
    std::map<std::pair<int, int>, Key> physicalGroups_to_key;
    gmsh::model::getPhysicalGroups(physicalGroups_dim_tag);

    Logger::info("Node number: " + std::to_string(mesh.n_node));


    Logger::info("Physical groups: Load elements...");
    // get all physical tags in dimension 0/1/2/3
    for(auto& [dim, tag] : physicalGroups_dim_tag)
    {
        std::string name;                                      // discription in gmsh
        gmsh::model::getPhysicalName(dim, tag, name);


        // generate next key for mesh.element_group
        mesh.dim_keys[dim].id++;
        Key gmsh_key = mesh.dim_keys[dim];
        


        // store key
        if(dim == mesh.dim_){
            mesh.key_domain.push_back(gmsh_key);
        }else if (dim == mesh.dim_-1){
            if(util::a_contains_b(name, {{"True", "Boundary"}})){
                mesh.key_true_boundary.push_back(gmsh_key);
            }else{
                mesh.key_internal_surface.push_back(gmsh_key); 
            }
        }else{
            mesh.key_others.push_back(gmsh_key); 
        }
        

        // get entity tags from physical groups
        std::vector<int> entity_tags;
        gmsh::model::getEntitiesForPhysicalGroup(dim, tag, entity_tags);

        auto& group = mesh.element_group[gmsh_key];
        mesh.element_group_description[gmsh_key] = name;

        // get all elements from entity
        for (int entity_tag : entity_tags)
        {   
            std::vector<int> element_types;
            std::vector<std::vector<size_t>> element_tags;
            std::vector<std::vector<size_t>> node_tags;
            // element_types and element_tags should have same size
            gmsh::model::mesh::getElements(element_types, element_tags, node_tags, dim, entity_tag);
            
            // one entity can be made of different type of elements, element_tags.size() is number of different shapes
            for (size_t i = 0; i < element_tags.size(); ++i) 
            {   
                auto [type, order, n_node] = convert_Type(element_types[i]);

                size_t n_elements = element_tags[i].size();
                group.reserve(group.size() + n_elements);

                // we pre-store all elements in vector for 3d elements
                if(mesh.dim_== dim){
                    for (size_t j = 0; j < n_elements; ++j)
                    {
                        auto it = element_map_md.find(element_tags[i][j]);  // element id = element_tags[i][j]
                        if (it != element_map_md.end()) {
                            group.push_back(it->second);
                        } else {
                            Logger::error("Mesh_Parser::load_gmsh - fatal"+std::to_string(dim)+"D elements not found in mesh.elements, return bad mesh.");
                        }
                    }
                }else{
                    for (size_t j = 0; j < n_elements; ++j)
                    {
                        auto startIt = node_tags[i].begin() + (j * n_node);
                        std::vector<size_t> element_node_tags(startIt, startIt + n_node);
                        group.push_back(create_element(type, element_node_tags, element_tags[i][j], gmsh_key.id, order));
                    }
                }
                
            }
        }


        Logger::mesh_entity(dim, tag, mesh.dim_keys[dim].id, mesh.element_group[gmsh_key].size(), name);
    }
    


    gmsh::finalize();
    //std::cout<<mesh.n_element<<std::endl;
    return mesh;
#else
    Logger::error("Mesh_Parser::load_gmsh - Gmsh support not compiled!");
    throw std::runtime_error("Gmsh support not compiled!");
#endif
}