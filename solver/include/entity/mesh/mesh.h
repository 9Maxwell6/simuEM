#pragma once
#include "element.h"
#include "node.h"

#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <functional>
#include <iostream>


struct Key 
{
    int dim;
    uint32_t id; 

    bool operator==(const Key& k) const {
        return dim == k.dim && id == k.id;
    }
};


struct Key_Hash {
    size_t operator()(const Key& k) const {
        return (static_cast<size_t>(k.dim) << 32) | k.id;
    }
};


class Mesh
{
    friend class Mesh_Parser;
protected:
    int dim_; // space dimension
    size_t n_node    , n_edge    , n_face    ;  // actual number of node/edge/face in mesh.
    //size_t N_dof_node, N_dof_edge, N_dof_face; 
    size_t n_element;

    size_t n_exterior_boundary_node; // true boundary of simulation domain.
    size_t n_interior_boundary_node;

    //std::vector<Element> elements;  // all elements.
    std::vector<Element *> elements;
    std::vector<Node> nodes;

    // 
    std::set<size_t> exterior_boundary_node;
    std::set<size_t> interior_boundary_node;
    
    // in 3D, the boundary is 2D Element
    // in 2D, the boundary is 1D Element
    std::vector<Element *> exterior_boundary_elements;  
    std::vector<Element *> interior_boundary_elements;
    

    int key_positive_idx = 0;
    int key_negative_idx = 0;
    
    // whenever distribute a key, key.id++.  
    // id=0 means bad key.
    Key dim0_key {0,0};  // key to find node group.   
    Key dim1_key {1,0};  // key to find edge group.
    Key dim2_key {2,0};  // key to find face group    (e.g., triangle).
    Key dim3_key {3,0};  // key to find volume group  (e.g., tetrahedron).
    Key dim_keys[4] = {dim0_key, dim1_key, dim2_key, dim3_key};
    std::unordered_map<Key, std::vector<Element *>, Key_Hash> element_group;
    std::unordered_map<Key, std::string, Key_Hash> group_description;
    


public:

    /**
     * @brief Creates a new named group of elements based on the Filter function.
     * 
     * Iterates over a pool of elements and applies a filter function to each.
     * Elements that pass the filter are collected into a new group, stored in
     * element_group with an auto-incremented key for the given dimension.
     * 
     * @tparam Filter  Callable with signature bool(const Element*).
     *                 Can be a lambda, function pointer, or functor.
     *                 Lambdas may capture the Mesh via [this] for access to mesh data.
     * 
     * @param filter       The filter function to apply to each element.
     * @param dim          Dimension of the new group (0=node, 1=edge, 2=triangle, 3=tet).
     * @param search_key   Key of an existing group to search within.
     *                     Default {0,0} searches all elements.
     * @param description  A label for the new group. Default "None".
     * 
     * @return Key of the newly created group, or {dim, 0} if:
     *         - the search pool is empty,
     *         - dim is out of range,
     *         - no elements matched the filter.
     */
    template <typename Filter>
    Key mark_elements(Filter&& filter, int dim, Key search_key = {0,0}, const std::string& description = "None")
    {
        // Determine which group of elements to search
        const std::vector<Element*>& search_pool = (search_key.id == 0) ? elements : element_group[search_key];

        if (search_pool.empty())
        {
            std::cerr << "Warning: search pool is empty for key {dim="<<search_key.dim<<", id=" << search_key.id << "}, return bad key.\n";
            return {dim, 0};
        }

        if (dim < 0 || dim > 3)
        {
            std::cerr << "Warning: impossible dimension " << dim << ", return bad key.\n";
            return {dim, 0};
        }

        std::vector<Element*> group;
        for (Element* e : search_pool)
        {
            if (filter(e))
            {
                group.push_back(e);
            }
        }
        
        if (!group.empty())
        {
            dim_keys[dim].id++;
            Key new_key = dim_keys[dim];
            element_group[new_key] = std::move(group);
            group_description[new_key] = description;
            return new_key;
        }

        std::cerr << "Warning: no elements matched the filter, return bad key.\n";
        return {dim, 0};
    }





};