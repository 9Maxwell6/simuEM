#pragma once

#include "world/mesh/mesh.h"
#include "world/structure/structure.h"

#include "math/fem/fem_space.h"
#include "math/fem/fem_system.h"
#include "math/fem/assemble_mat.h"
#include "math/fem/assemble_vec.h"

#include "math/fem/bc_dirichlet.h"

#include "math/fem/post_processing.h"


#include "math/field/field_function.h"

#include "utils/util_string.h"
#include "utils/util_constant.h"



#include <functional>
#include <unordered_map>

namespace simu {

enum Domain 
{
    EMPTY = 1,
    CONDUCTOR = 2,
    INSULATOR = 3,
    SOURCE = 4,
    CONDUCTOR_OUTER_LAYER = 5
};

class T_Omega
{

private:
    Mesh& mesh_;
    FEM_System fe_system_;

    H1_Space           Omega_space_;
    Block              dof_Omega_;

    Hcurl_Space        T_space_1_;
    Block              dof_T_1_;

    Block              dof_coupling_1_;
    Block              dof_coupling_tp_1_;

    H1_Space           T_H1_v_; //  vector H1 field in T
    H1_Space           T_H1_s_; //  scalar H1 field in T
    Block              pc_P_;   //  edge interpolation matrix for preconditioner.
    Block              pc_G_;   //  discrete gradient matrix  for preconditioner.
    Block              pc_L_;   //  L = integral_{ σ^-1  ∇u : ∇v  + μ u.v  dx }    (H1)^3 x (H1)^3
    Block              pc_Q_;   //  Q = integral_{ μ grad u. grad v  dx }           H1   x   H1


    Dirichlet_BC bc_Omega_out_;
    Dirichlet_BC bc_Omega_in_;

    Dirichlet_BC bc_T_1_;


    Block_Rack br_system_;


    int dim_;

    Key key_true_boundary_;                // key to 1D/2D true boundary element groups

    Key key_Omega_inner_boundary_1_;       // key to 1D/2D omega field boundary element groups


    Key key_conductor_interface_1_;        // key to 1D/2D conductor boundary element groups
    Key key_conductor_interface_layer_1_;  // key to 2D/3D conducting region element groups in contact with conductor boundarys.

    Key key_source_;                       // key to 2D/3D source element groups
    Key key_insulator_;                    // key to 2D/3D insulating region element groups
    Key key_conductor_1_;                  // key to 2D/3D conducting region element groups
    

    Key key_Omega_;

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
    std::function<bool(Element*)> conductor_interface_layer_filter(Key& key_conductor_interface) {
        return [this, &key_conductor_interface](Element* e) -> bool {
            const size_t * e_ids = e->get_node_idx();
            int e_size = e->get_node_num();
            //std::cout<<e_ids[0]<<" "<< e_ids[1]<<" " << e_ids[2]<<" "<< e_ids[3]<<" " <<std::endl;
            for (Element* ie : this->mesh_.get_element_group(key_conductor_interface))
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
    std::function<std::vector<Element *>(Element*)> scalar_field_Omega_inner_boundary_filter(Key& key_conductor_interface) {
        return [this, &key_conductor_interface](Element* e) -> std::vector<Element *> {
            const size_t * e_ids = e->get_node_idx();
            int e_size = e->get_node_num();
            std::vector<size_t> exclude_ids;
            // get all node index on conductor interface
            for (size_t id : this->mesh_.get_node_group(key_conductor_interface))
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
            
            return {};
        };
    }



public:
    T_Omega(Mesh& mesh);

    bool assemble_system();

    bool assemble_preconditioner();

    bool solve_system();

    scalar_t compute_L2_error();

};


}