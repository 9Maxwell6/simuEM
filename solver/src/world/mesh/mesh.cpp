#include "world/mesh/mesh.h"

using namespace simu;

Mesh::~Mesh() 
{   
    // Delete all element_group of dimensions other than mesh dimension,
    // since elements with mesh dimension are directly stored in vector.
    for (auto& [key, elements] : element_group)
    {
        if (key.dim != dim_) {
            for (auto* e : elements) delete e;
        }
    }
}

const std::vector<Element *>& Mesh::get_mesh_elements() const { return elements_; }
const std::map<Geometry, size_t>& Mesh::get_mesh_element_geometry_size() const { return geometry_size_; }


const std::vector<Element *>& Mesh::get_element_group(const Key& mesh_key) const
{
    auto it = element_group.find(mesh_key);
    if (it != element_group.end()) return it->second;

    Logger::error("Mesh::get_element_group - failed: key not found in mesh, return reference to empty vector<Element *>.");
    static const std::vector<Element*> empty;
    return empty;
}

const std::map<Geometry, size_t>& Mesh::get_element_geometry_size_group(const Key& mesh_key) const
{
    auto it = element_geometry_size_group.find(mesh_key);
    if (it != element_geometry_size_group.end()) return it->second;

    Logger::error("Mesh::get_element_geometry_size_group - failed: key not found in mesh, return reference to empty std::set<Geometry>.");
    static const std::map<Geometry, size_t> empty;
    return empty;
}

const std::array<size_t, 4>& Mesh::get_element_size_group(const Key& mesh_key) const
{
    auto it = element_size_group.find(mesh_key);
    if (it != element_size_group.end()) return it->second;

    Logger::error("Mesh::get_element_size_group - failed: key not found in mesh, return reference to empty std::array<size_t, 4>.");
    static const std::array<size_t, 4> empty = {0,0,0,0};
    return empty;
}


/**
 * @brief Get all nodal indices of elements in group with given mesh_key.
 *
 * @param mesh_key (lambda function parameter) pointer to the target element.
 * @return (lambda function return) true if at least one element in conductor_interface is covered by the target element,
 * otherwise return false.
 */
const std::set<size_t>& Mesh::get_node_group(const Key& mesh_key)
{
    // TODO: std::set has poor locality, change to other array structure.
    auto node_it = node_group.find(mesh_key);
    if (node_it != node_group.end()) {
        return node_it->second;
    }

    auto& group = node_group[mesh_key];

    auto element_it = element_group.find(mesh_key);
    if (element_it != element_group.end()){
        for(Element * e : element_it->second)
        {
            const size_t * node_ids =  e->get_node_idx();
            int nodes_size = e->get_node_num();
            for (size_t i = 0; i < nodes_size; ++i) 
            {
                group.insert(node_ids[i]);
            }

        }
        return group;
    }

    Logger::error("Mesh::get_node_group - failed: key not found in mesh, return reference to empty vector<size_t>.");
    static const std::set<size_t> empty;
    return empty;
}


const std::string& Mesh::get_group_description(const Key& mesh_key) const
{ 
    auto it = element_group_description.find(mesh_key);
    if (it != element_group_description.end())
        return it->second;
    Logger::error("Mesh::get_group_description - failed: key not found in mesh, return reference to empty description.");
    static const std::string empty = "";
    return empty;
}


void Mesh::set_group_property_id(const Key& mesh_key, size_t property_id)
{
    if(property_id == 0) Logger::warning("Mesh::set_group_property_id - set to default id=0, not recommanded.");

    auto it = element_group.find(mesh_key);
    if (it != element_group.end()){
        for(Element* e : it->second)
        {
            size_t old_id = e->get_property_id();
            if(old_id !=0 ) Logger::warning("Mesh::set_group_property_id - overwrite element property id: old id: "+std::to_string(old_id)+"; new id: "+std::to_string(property_id)+".");
            e->set_property_id(property_id);
        }
    }else{
        Logger::error("Mesh::set_group_property_id - failed: key not found in mesh.");
    }

    
}


Element * Mesh::create_element(Geometry element_type, std::vector<std::size_t> node_idx, size_t element_id, size_t property_id, int o)
{
    switch (element_type) 
    {
        case Geometry::EDGE: return new Edge(node_idx, element_id, property_id, o);
        case Geometry::TRIANGLE: return new Triangle(node_idx, element_id, property_id, o);
        case Geometry::TETRAHEDRON: return new Tetrahedron(node_idx, element_id, property_id, o);
        default: 
            Logger::warning("Mesh_Parser::create_element - failed: element_type not available, return nullptr.");
            return nullptr;
    }

}

Key Mesh::group_union(const Key& group_1_key, const Key& group_2_key, const std::string& description)
{
    if(group_1_key.dim != group_2_key.dim) { Logger::error("Key Mesh::group_union - failed: group dimension not match, return bad key."); return {0,0}; }

    auto it_1 = element_group.find(group_1_key);
    if (it_1 == element_group.end()){ Logger::warning("Mesh::group_union - group_1_key not found, return bad key."); return {0,0};}

    auto it_2 = element_group.find(group_2_key);
    if (it_2 == element_group.end()){ Logger::warning("Mesh::group_union - group_2_key not found, return bad key."); return {0,0};}

    const std::vector<Element*>& group_1 =it_1->second;
    const std::vector<Element*>& group_2 =it_2->second;

    dim_keys[group_1_key.dim].id++;
    Key new_key = dim_keys[group_1_key.dim];
    element_group_description[new_key] = description;

    std::vector<Element*>& group_1_2 = element_group[new_key];
    std::map<Geometry,size_t>& g_group_1_2 = element_geometry_size_group[new_key];
    
    group_1_2.insert(group_1_2.end(), group_1.begin(), group_1.end());

    util::Block_Hash_D bh(32*1024);

    for (auto* e : group_1)
    {
        const size_t* idx  = e->get_node_idx();
        int n = e->get_node_num();
        bh.get_id(idx, n, 0);
        g_group_1_2[e->get_geometry()]++;
    }

    for (auto* e : group_2)
    {
        const size_t* idx  = e->get_node_idx();
        int n = e->get_node_num();
        if(!bh.if_exist(idx, n, 0))
        {
            group_1_2.push_back(e);
            g_group_1_2[e->get_geometry()]++;
        }

    }

    element_size_group[new_key] = count_node_edge_face_volume(element_group[new_key]);
    

    return new_key;
}

/**
 * @brief Count unique topological entities in a mesh.
 *
 * Traverses all elements and counts unique nodes, edges, faces, and volumes
 * without double-counting shared entities between adjacent elements.
 * Deduplication is based on sorted corner node indices, so curved/high-order
 * elements are handled correctly (mid-edge nodes are ignored).
 *
 * @param elements List of mesh elements.
 * @return std::array<size_t, 4> where:
 *         [0] = number of unique nodes
 *         [1] = number of unique edges
 *         [2] = number of unique faces
 *         [3] = number of volumes
 */
std::array<size_t, 4> Mesh::count_node_edge_face_volume(const std::vector<Element*>& elements) {

    util::Block_Hash bh(1024, 1024, 1024, 1024);
    size_t unique_node_count = 0;
    size_t unique_edge_count = 0;
    size_t unique_face_count = 0;
    size_t unique_volume_count = 0;

    //std::set<size_t>              unique_nodes;
    //std::set<std::vector<size_t>> unique_edges;
    //std::set<std::vector<size_t>> unique_faces;
    //size_t                     num_volumes = 0;

    //auto makeKey = [](std::initializer_list<size_t> list) {
    //    std::vector<size_t> v(list);
    //    std::sort(v.begin(), v.end());
    //    return v;
    //};

    for (auto* e : elements) {
        const size_t* idx  = e->get_node_idx();
        int n = e->get_node_num();

        auto g = e->get_geometry();

        //for (int i = 0; i < n; ++i) unique_nodes.insert(idx[i]);

        for (int i = 0; i < n; ++i) 
            if(!bh.if_exist(idx[i], 0)) { unique_node_count++; bh.get_id(idx[i], 0); }

        switch (g) {
            case Geometry::EDGE:
                if(!bh.if_exist(idx[0], idx[1], 0)) { unique_node_count++; bh.get_id(idx[0], idx[1], 0); }

                //unique_edges.insert(makeKey({idx[0], idx[1]}));
                break;

            case Geometry::TRIANGLE:
                if(!bh.if_exist(idx[0], idx[1], 0)) { unique_edge_count++; bh.get_id(idx[0], idx[1], 0); }
                if(!bh.if_exist(idx[1], idx[2], 0)) { unique_edge_count++; bh.get_id(idx[1], idx[2], 0); }
                if(!bh.if_exist(idx[0], idx[2], 0)) { unique_edge_count++; bh.get_id(idx[0], idx[2], 0); }
                if(!bh.if_exist(idx[0], idx[1], idx[2], 0)) { unique_face_count++; bh.get_id(idx[0], idx[1], idx[2], 0); }
                //unique_edges.insert(makeKey({idx[0], idx[1]}));
                //unique_edges.insert(makeKey({idx[1], idx[2]}));
                //unique_edges.insert(makeKey({idx[0], idx[2]}));
                //unique_faces.insert(makeKey({idx[0], idx[1], idx[2]}));
                break;

            case Geometry::TETRAHEDRON:
                if(!bh.if_exist(idx[0], idx[1], 0)) { unique_edge_count++; bh.get_id(idx[0], idx[1], 0); }
                if(!bh.if_exist(idx[0], idx[2], 0)) { unique_edge_count++; bh.get_id(idx[0], idx[2], 0); }
                if(!bh.if_exist(idx[0], idx[3], 0)) { unique_edge_count++; bh.get_id(idx[0], idx[3], 0); }
                if(!bh.if_exist(idx[1], idx[2], 0)) { unique_edge_count++; bh.get_id(idx[1], idx[2], 0); }
                if(!bh.if_exist(idx[1], idx[3], 0)) { unique_edge_count++; bh.get_id(idx[1], idx[3], 0); }
                if(!bh.if_exist(idx[2], idx[3], 0)) { unique_edge_count++; bh.get_id(idx[2], idx[3], 0); }
                if(!bh.if_exist(idx[0], idx[1], idx[2], 0)) { unique_face_count++; bh.get_id(idx[0], idx[1], idx[2], 0); }
                if(!bh.if_exist(idx[0], idx[1], idx[3], 0)) { unique_face_count++; bh.get_id(idx[0], idx[1], idx[3], 0); }
                if(!bh.if_exist(idx[0], idx[2], idx[3], 0)) { unique_face_count++; bh.get_id(idx[0], idx[2], idx[3], 0); }
                if(!bh.if_exist(idx[1], idx[2], idx[3], 0)) { unique_face_count++; bh.get_id(idx[1], idx[2], idx[3], 0); }
                unique_volume_count++;

                //unique_edges.insert(makeKey({idx[0], idx[1]}));
                //unique_edges.insert(makeKey({idx[0], idx[2]}));
                //unique_edges.insert(makeKey({idx[0], idx[3]}));
                //unique_edges.insert(makeKey({idx[1], idx[2]}));
                //unique_edges.insert(makeKey({idx[1], idx[3]}));
                //unique_edges.insert(makeKey({idx[2], idx[3]}));
                //unique_faces.insert(makeKey({idx[0], idx[1], idx[2]}));
                //unique_faces.insert(makeKey({idx[0], idx[1], idx[3]}));
                //unique_faces.insert(makeKey({idx[0], idx[2], idx[3]}));
                //unique_faces.insert(makeKey({idx[1], idx[2], idx[3]}));
                //num_volumes++;
                break;


            default: break;
        }
    }

    return {unique_node_count, unique_edge_count, unique_face_count, unique_volume_count};
}



/**
 * @brief Create sub-element for the given element.
 * 
 * e.g., Tetrahedron element consists of 4 Triangle elements, 6 Edge elements.
 * 
 * Based on the exclude_id and dim parameter, the function will return the sub-elements
 * of the given dimension which do not contain any id listed in exclude_id.
 *
 * @param e Pointer to the target element.
 * @param exclude_id List of ids that should not appear in sub-elements.
 * @param dim Desired dimension of sub-elements.
 * @return List of sub-elements.
 */
std::vector<Element *> Mesh::create_sub_element(Element * e, std::vector<size_t>& exclude_ids, int dim)
{
    int element_id = e->get_id();
    int property_id = e->get_property_id();
    int order = e->get_geometry_order();

    const size_t * node_ids =  e->get_node_idx();
    int nodes_size = e->get_node_num();

    std::vector<Element *> new_sub_elements;
    
    switch (e->get_geometry()) 
    {
        case Geometry::EDGE: 
        {
            Logger::warning("No need to create sub-elements for Geometry::EDGE element, return empty vector.");
            return {};
        }
        case Geometry::TRIANGLE: 
        {
            // sub-element -> edges
            if(dim==1){
                size_t n0 = node_ids[0];
                size_t n1 = node_ids[1]; 
                size_t n2 = node_ids[2];

                // by convension, node index from small to large
                if(util::a_not_in_b({n0, n1}, exclude_ids))
                    new_sub_elements.push_back(new Edge(util::sort_a_b(n0, n1), element_id, property_id, order));
                if(util::a_not_in_b({n0, n2}, exclude_ids))
                    new_sub_elements.push_back(new Edge(util::sort_a_b(n0, n2), element_id, property_id, order));
                if(util::a_not_in_b({n1, n2}, exclude_ids))
                    new_sub_elements.push_back(new Edge(util::sort_a_b(n1, n2), element_id, property_id, order));
                return new_sub_elements;
            }

            Logger::warning("No need to create sub-elements of dimension-"+std::to_string(dim)+" for Geometry::TRIANGLE element, return empty vector.");
            return {};
        }
        case Geometry::TETRAHEDRON: 
        {
            size_t n0 = node_ids[0];
            size_t n1 = node_ids[1]; 
            size_t n2 = node_ids[2];
            size_t n3 = node_ids[3];
            // sub-element -> edges
            if(dim==1){
                if(util::a_not_in_b({n0, n1}, exclude_ids))
                    new_sub_elements.push_back(new Edge(util::sort_a_b(n0, n1), element_id, property_id, order));
                if(util::a_not_in_b({n0, n2}, exclude_ids))
                    new_sub_elements.push_back(new Edge(util::sort_a_b(n0, n2), element_id, property_id, order));
                if(util::a_not_in_b({n0, n3}, exclude_ids))
                    new_sub_elements.push_back(new Edge(util::sort_a_b(n0, n3), element_id, property_id, order));
                if(util::a_not_in_b({n1, n2}, exclude_ids))
                    new_sub_elements.push_back(new Edge(util::sort_a_b(n1, n2), element_id, property_id, order));
                if(util::a_not_in_b({n1, n3}, exclude_ids))
                    new_sub_elements.push_back(new Edge(util::sort_a_b(n1, n3), element_id, property_id, order));
                if(util::a_not_in_b({n2, n3}, exclude_ids))
                    new_sub_elements.push_back(new Edge(util::sort_a_b(n2, n3), element_id, property_id, order));
                return new_sub_elements;
            }
            
            // sub-element -> triangles
            if(dim==2){
                const Node& p0 = nodes_[n0];
                const Node& p1 = nodes_[n1];
                const Node& p2 = nodes_[n2];
                const Node& p3 = nodes_[n3];
                if(util::a_not_in_b({n0, n1, n2}, exclude_ids)) 
                    new_sub_elements.push_back(new Triangle(util::sort_outward_triangle(n0, n1, n2, p0, p1, p2, p3), element_id, property_id, order));
                if(util::a_not_in_b({n0, n1, n3}, exclude_ids)) 
                    new_sub_elements.push_back(new Triangle(util::sort_outward_triangle(n0, n1, n3, p0, p1, p3, p2), element_id, property_id, order));
                if(util::a_not_in_b({n0, n2, n3}, exclude_ids)) 
                    new_sub_elements.push_back(new Triangle(util::sort_outward_triangle(n0, n2, n3, p0, p2, p3, p1), element_id, property_id, order));
                if(util::a_not_in_b({n1, n2, n3}, exclude_ids)) 
                    new_sub_elements.push_back(new Triangle(util::sort_outward_triangle(n1, n2, n3, p1, p2, p3, p0), element_id, property_id, order));
                return new_sub_elements;
            }

            Logger::warning("No need to create sub-elements of dimension-"+std::to_string(dim)+" for Geometry::TETRAHEDRON element, return empty vector.");
            return {};
            //return new Tetrahedron(node_idx, element_id, property_id, order);
        }
        default: 
            Logger::warning("Mesh_Parser::create_element - failed: element_type not available, return empty vector");
            return {};
    }
}

