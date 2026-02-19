#include "entity/mesh/mesh.h"

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

const std::vector<Element *>& Mesh::get_element_group(Key mesh_key) const
{
    auto it = element_group.find(mesh_key);
    if (it != element_group.end())
        return it->second;
    Logger::error("Mesh::get_element_group - failed: key not found in mesh, return reference to empty vector<Element *>.");
    static const std::vector<Element*> empty;
    return empty;
}


/**
 * @brief Get all nodal indices of elements in group with given mesh_key.
 *
 * @param mesh_key (lambda function parameter) pointer to the target element.
 * @return (lambda function return) true if at least one element in conductor_interface is covered by the target element,
 * otherwise return false.
 */
const std::set<size_t>& Mesh::get_node_group(Key mesh_key)
{
    auto node_it = node_group.find(mesh_key);
    if (node_it != node_group.end()) {
        return node_it->second;
    }

    auto& group = node_group[mesh_key];

    auto element_it = element_group.find(mesh_key);
    if (element_it != element_group.end()){
        for(Element * e : element_it->second)
        {
            const size_t * node_ids =  e->get_nodeIdx();
            int nodes_size = e->get_nodeNum();
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


const std::string& Mesh::get_group_description(Key mesh_key) const
{ 
    auto it = element_group_description.find(mesh_key);
    if (it != element_group_description.end())
        return it->second;
    Logger::error("Mesh::get_group_description - failed: key not found in mesh, return reference to empty description.");
    static const std::string empty = "";
    return empty;
}


Element * Mesh::create_element(Type element_type, std::vector<std::size_t> node_idx, size_t element_id, size_t property_id, int o)
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
    int element_id = e->get_Id();
    int property_id = e->get_propertyId();
    int order = e->get_geometry_order();

    const size_t * node_ids =  e->get_nodeIdx();
    int nodes_size = e->get_nodeNum();

    std::vector<Element *> new_sub_elements;
    
    switch (e->get_Type()) 
    {
        case Type::EDGE: 
        {
            Logger::warning("No need to create sub-elements for Type::EDGE element, return empty vector.");
            return {};
        }
        case Type::TRIANGLE: 
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

            Logger::warning("No need to create sub-elements of dimension-"+std::to_string(dim)+" for Type::TRIANGLE element, return empty vector.");
            return {};
        }
        case Type::TETRAHEDRON: 
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

            Logger::warning("No need to create sub-elements of dimension-"+std::to_string(dim)+" for Type::TETRAHEDRON element, return empty vector.");
            return {};
            //return new Tetrahedron(node_idx, element_id, property_id, order);
        }
        default: 
            Logger::warning("Mesh_Parser::create_element - failed: element_type not available, return empty vector");
            return {};
    }
}

