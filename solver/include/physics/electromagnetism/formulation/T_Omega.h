#pragma once

#include "entity/mesh/mesh.h"
#include "math/fem/fem_space.h"
#include "math/fem/fem_system.h"
#include "math/fem/assemble.h"

#include "math/field/field_function.h"

#include "utils/util_string.h"


#include <functional>

namespace simu {


class T_Omega
{

private:
    Mesh& mesh_;
    FEM_System fe_system_;
    Block_Rack br_l;

    H1_Space           Omega_space_;
    Block              dof_Omega_;

    Hcurl_Space        T_space_;
    std::vector<Block> dof_T_;

    std::vector<Block> dof_coupling_;
    std::vector<Block> dof_coupling_tp_;




    int dim_;

    std::vector<Key> key_true_boundary;              // key to 1D/2D true boundary element groups
    std::vector<Key> key_Omega_field_boundary;       // key to 1D/2D omega field boundary element groups
    std::vector<Key> key_conductor_interface;        // key to 1D/2D conductor boundary element groups
    std::vector<Key> key_source;                     // key to 2D/3D source element groups
    std::vector<Key> key_insulator;                  // key to 2D/3D insulating region element groups
    std::vector<Key> key_conductor;                  // key to 2D/3D conducting region element groups
    std::vector<Key> key_conductor_interface_layer;  // key to 2D/3D conducting region element groups in contact with conductor boundarys.
    


    /**
     * @brief Wapper of lambda filter function used for Mesh::mark_elements, it will checks if 
     * at least one element in the conductor_interface group is covered by the target element
     *
     * element A is considered "covered" by element B if there is at least one A's vertices
     * contained in B's vertices.
     *
     * @param e (lambda function parameter) pointer to the target element.
     * @return (lambda function return) true  If at least one element in conductor_interface is covered by the target element,
     * otherwise return false.
     */
    std::function<bool(Element*)> conductor_interface_layer_filter() {
        return [this](Element* e) -> bool {
            const size_t * e_ids = e->get_node_idx();
            int e_size = e->get_node_num();
            //std::cout<<e_ids[0]<<" "<< e_ids[1]<<" " << e_ids[2]<<" "<< e_ids[3]<<" " <<std::endl;
            for(Key& key: this->key_conductor_interface)
            {    
                for (Element* ie : this->mesh_.get_element_group(key))
                {
                    const size_t * ie_ids = ie->get_node_idx();
                    int ie_size = ie->get_node_num();

                    for (size_t i = 0; i < ie_size; ++i) 
                    {
                        for (size_t j = 0; j < e_size; ++j) 
                        {
                            if (ie_ids[i] == e_ids[j]) { return true; }
                        }
                    }
                }
            }
            return false;
        };
    }


    /**
     * @brief Wapper of lambda filter function used for Mesh::mark_elements, it will help to mark the
     * inner boundary elements for Omega field (scalar filed). Omega field covers the
     * insulating region and the outer layer of conducting region, the inner boundary 
     * elements correspond to the boundary of outer conducting layers inside the conductor.
     *
     *
     * @param e (lambda function return) pointer to the target element.
     * @return (lambda function return) pointer to the surface element if exist, otherwise return nullptr.
     */
    std::function<std::vector<Element *>(Element*)> scalar_field_Omega_inner_boundary_filter() {
        return [this](Element* e) -> std::vector<Element *> {
            const size_t * e_ids = e->get_node_idx();
            int e_size = e->get_node_num();

            for(Key& key: this->key_conductor_interface)
            {    
                std::vector<size_t> exclude_ids;
                // get all node index on conductor interface
                for (size_t id : this->mesh_.get_node_group(key))
                {
                    for (size_t i = 0; i < e_size; ++i) 
                    {      
                        if (id == e_ids[i]) { exclude_ids.push_back(id); }                 
                    }
                }
                if(exclude_ids.size()!=0) {
                    std::vector<Element *> new_element = this->mesh_.create_sub_element(e, exclude_ids, dim_-1);
                    if(new_element.size()!=0) return new_element;
                }
            }
            return {};
        };
    }



public:
    T_Omega(Mesh& mesh);

    bool assemble_system();

};


}