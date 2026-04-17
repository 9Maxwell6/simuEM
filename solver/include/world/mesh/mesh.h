#pragma once

#include "world/mesh/e_collection.h"
#include "world/mesh/group_key.h"

#include "utils/logger.h"
#include "utils/util_math.h"
#include "utils/util_hash.h"



#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <functional>
#include <iostream>

namespace simu {




class Mesh
{
    friend class Mesh_Parser;

protected:
    int dim_; // space dimension

    // total number of node/edge/face in mesh.
    size_t n_node_;
    size_t n_edge_;
    size_t n_face_; 
    size_t n_volume_;

    size_t n_exterior_boundary_node; // true boundary of simulation domain.
    size_t n_interior_boundary_node;
    

    std::vector<Node> nodes_;       // 0d element: contain pure 3d coordinates, node index is its index in this vector.
    std::vector<Edge> elements_edge_;
    std::vector<Triangle> elements_triangle_;
    //std::vector<Quadrilateral> elements_quadrilateral;  // not implemented.
    std::vector<Tetrahedron> elements_tetrahedron_;

    // TODO: For higher order elements, indices of extra nodes other than vertices will stored here.
    // not implemented.
    // 

    
    //std::vector<Element*> curve;    // 1d element: contain pointers to 2d elements. e.g. in std::vector<Edge>
    //std::vector<Element*> surface;  // 2d element: contain pointers to 2d elements. e.g. in std::vector<Triangle>, std::vector<Quadrilateral>
    //std::vector<Element*> volume;   // 3d element: contain index to nodes

    // all elements with the same dimension of mesh
    std::vector<Element*> elements_;
    // all geometry types of the elements, and number of elements per type
    std::map<Geometry, size_t> geometry_size_;

    // count size of faces for each different number of vertices, used for constructing dof structure.
    //  map<#vertices, #faces>
    //std::map<size_t, size_t> face_size_;



    

    int key_positive_idx = 0;
    int key_negative_idx = 0;
    
    // whenever distribute a key, key.id++.  
    // id=0 means bad key.
    //Key dim0_key {0,0};  // key to find node group.   
    //Key dim1_key {1,0};  // key to find edge group.
    //Key dim2_key {2,0};  // key to find face group    (e.g., triangle).
    //Key dim3_key {3,0};  // key to find volume group  (e.g., tetrahedron).
    //Key dim_keys[4] = {dim0_key, dim1_key, dim2_key, dim3_key};

    // id=0 means bad key.
    // usage: whenever distribute a key, key.id++.  
    Key dim_keys[4] = {{0,0},    // key to find node group.  
                       {1,0},    // key to find curve group.
                       {2,0},    // key to find surface group    (e.g., triangle).
                       {3,0}};   // key to find volume group     (e.g., tetrahedron).



    // caution !!!
    //      for d-dimensional mesh, the element group of dimension < d
    //      can have duplicate elements, this is due to the perpose of keeping info related to its parent element. 
    //      (e.g. surface of triangle has normal pointing outward from tetrahedron).
    //      
    //      hence for element_geometry_size_group, the size of each geometry can double count elements,
    //      if the dimension of this group < dimension of the mesh.
    std::unordered_map<Key, std::vector<Element *>,      Key::Hash> element_group;                // group of elements
    std::unordered_map<Key, std::map<Geometry, size_t>,  Key::Hash> element_geometry_size_group;  // element geometry from the group, and number of elements per type
    std::unordered_map<Key, std::set<size_t>,            Key::Hash> node_group;                   // all nodes from the group
    std::unordered_map<Key, std::string,                 Key::Hash> element_group_description;    // description of the group

    std::unordered_map<Key, std::array<size_t, 4>,       Key::Hash> element_size_group;           // total number of node/edge/face/volume of the group
    std::unordered_map<Key, std::map<size_t, size_t>,    Key::Hash> element_face_size_group;      // total number of of faces for each different number of vertices of the group

    std::vector<Key> key_true_boundary;
    std::vector<Key> key_internal_surface; 
    std::vector<Key> key_domain;
    std::vector<Key> key_others;


    std::array<size_t, 4> count_node_edge_face_volume(const std::vector<Element*>& elements);
    


public:

    int get_mesh_dimension() const {return dim_; }

    const Node& get_node(size_t idx) const { return nodes_[idx]; }

    const std::vector<Element *>& get_mesh_elements() const;
    const std::map<Geometry, size_t>& get_mesh_element_geometry_size() const;


    const std::vector<Element *>&     get_element_group(const Key& mesh_key) const;
    const std::map<Geometry, size_t>& get_element_geometry_size_group(const Key& mesh_key) const;
    const std::array<size_t, 4>&      get_element_size_group(const Key& mesh_key) const;
    const std::set<size_t>&           get_node_group(const Key& mesh_key);
    const std::string&                get_group_description(const Key& mesh_key) const;


    const std::vector<Key>& get_keys_true_boundary() const { return key_true_boundary; }
    const std::vector<Key>& get_keys_internal_surface() const { return key_internal_surface; }
    const std::vector<Key>& get_keys_domain() const { return key_domain; }
    const std::vector<Key>& get_keys_others() const { return key_others; }

    void set_group_property_id(const Key& mesh_key, size_t property_id);

    Element * create_element(Geometry element_type, std::vector<std::size_t> node_idx, size_t element_id, size_t property_id, int o);


    std::vector<Element *> create_sub_element(Element * e, std::vector<size_t>& exclude_ids, int dim);

    // TODO：
    Key group_intersection(const Key& group_1, const Key& group_2, const std::string& description);
    Key group_union(const Key& group_1, const Key& group_2, const std::string& description);
    void delete_group(const Key& group);


    ~Mesh();


    /**
     * @brief Creates a new named group of elements based on the Filter function.
     * 
     * Iterates over a pool of elements and applies a filter function to each.
     * Elements that pass the filter are collected into a new group, stored in
     * element_group with an auto-incremented key.
     * 
     * The dimension of the new element group must be same with the element group 
     * provided by search_key.
     * 
     * @tparam Filter  Callable with signature bool(const Element*).
     *                 Can be a lambda (recommand) or function pointer.
     *                 Lambdas may capture the Mesh via [this] for access to mesh data.
     * 
     * @param filter       The filter function to apply to each element.
     * @param search_key   Key of an existing group to search within.
     *                     Default {0,0} searches all elements.
     * @param description  A label for the new group. Default "None".
     * @param property_id  Element property id, used for applying coefficients.
     * 
     * @return Key of the newly created group, or {dim, 0} if:
     *         - the search pool is empty,
     *         - no elements matched the filter.
     */
    template <typename Filter>
    Key mark_elements(Filter&& filter, Key search_key = {0,0}, const std::string& description = "None", size_t property_id=0)
    {
        auto it = element_group.find(search_key);
        if (it == element_group.end()) Logger::info("Mesh::mark_elements - search_key not found: search from all elements with the same dimension of mesh.");

        // Determine which group of elements to search
        const std::vector<Element*>& search_pool = (it != element_group.end()) ? it->second : elements_;
        
        //const std::vector<Element*>& search_pool; = (search_key.id == 0) ? elements : element_group[search_key];

        if (search_pool.empty())
        {
            Logger::error("Mesh::mark_elements - failed: search pool is empty for key {dim="+ std::to_string(search_key.dim)+", id=" +std::to_string(search_key.id) + "}, return bad key.\n");
            return {static_cast<uint32_t>(search_key.dim), 0};
        }


        std::vector<Element*> e_group;
        std::map<Geometry, size_t> g_group;

        for (Element* e : search_pool)
        {
            if (filter(e))
            {
                if(property_id!=0) e->set_property_id(property_id);
                e_group.push_back(e);
                g_group[e->get_geometry()]++;
            }
        }
        
        if (!e_group.empty())
        {
            dim_keys[search_key.dim].id++;
            Key new_key = dim_keys[search_key.dim];
            
            element_group[new_key] = std::move(e_group);

            element_geometry_size_group[new_key] = std::move(g_group);

            element_size_group[new_key] = count_node_edge_face_volume(element_group[new_key]);

            element_group_description[new_key] = description;

            return new_key;
        }

        Logger::error("Mesh::mark_elements - failed: no elements matched the filter, return bad key.\n");
        return {static_cast<uint32_t>(search_key.dim), 0};
    }



    /**
     * @brief Creates a new named group of elements based on the Filter function,
     * these marked elements are previously not contained in any group, thus they
     * are newly created within this function.
     * 
     * e.g., marking certain faces of tetrahedron, which previously not initialized 
     * during mesh loading phase. 
     * 
     * Iterates over a pool of elements and applies a filter function to each.
     * 
     * The surface/edge elements of the that pass the filter are collected into a new group, stored in
     * element_group with an auto-incremented key for the given dimension.
     * 
     * @tparam Filter  Callable with signature Element*(const Element*).
     *                 Can be a lambda (recommand) or function pointer.
     *                 Lambdas may capture the Mesh via [this] for access to mesh data.
     * 
     * @param filter       The filter function to apply to each element.
     * @param dim          Dimension of the new group
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
    Key mark_new_elements(Filter&& filter, int dim, Key search_key = {0,0}, const std::string& description = "None", size_t property_id=0)
    {

        if (dim < 0 || dim > dim_)
        {
            Logger::error("Mesh::mark_elements - failed: impossible dimension "+std::to_string(dim) + ", return bad key.\n");
            return {static_cast<uint32_t>(dim), 0};
        }


        auto it = element_group.find(search_key);
        if (it == element_group.end()) Logger::info("Mesh::mark_elements - search_key not found: search from all elements with the same dimension of mesh.");

        // Determine which group of elements to search
        const std::vector<Element*>& search_pool = (it != element_group.end()) ? it->second : elements_;
        
        //const std::vector<Element*>& search_pool; = (search_key.id == 0) ? elements : element_group[search_key];

        if (search_pool.empty())
        {
            Logger::error("Mesh::mark_elements - failed: search pool is empty for key {dim="+ std::to_string(search_key.dim)+", id=" +std::to_string(search_key.id) + "}, return bad key.\n");
            return {static_cast<uint32_t>(dim), 0};
        }


        std::vector<Element*> e_group;
        std::map<Geometry, size_t> g_group;

        for (Element* e : search_pool)
        {
            std::vector<Element *> new_es = filter(e);
            for(Element * new_e : new_es){
                // TODO: check if new element is duplicate ?  (probobly not, element contains exterior normal info)
                if(property_id!=0) new_e->set_property_id(property_id);
                e_group.push_back(new_e);
                g_group[new_e->get_geometry()]++;
                const size_t * e_ids = new_e->get_node_idx();
            }

        }
        
        if (!e_group.empty())
        {
            dim_keys[dim].id++;
            Key new_key = dim_keys[dim];
            
            element_group[new_key] = std::move(e_group);

            element_geometry_size_group[new_key] = std::move(g_group);

            element_size_group[new_key] = count_node_edge_face_volume(element_group[new_key]);

            element_group_description[new_key] = description;

            return new_key;
        }

        Logger::error("Mesh::mark_elements - failed: no elements matched the filter, return bad key.\n");
        return {static_cast<uint32_t>(dim), 0};
    }


};



}