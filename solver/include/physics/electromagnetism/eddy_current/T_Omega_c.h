#pragma once

#include "world/mesh/mesh.h"
#include "world/mesh/mesh_parser.h"

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

#include "physics/electromagnetism/eddy_current/auxiliary_solver_c.h"




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

class T_Omega_c
{

private:

    // coefficient:
    //coefficient:
    double mu_conductor_ = 1.;
    double mu_insulator_ = 1.;
    double sigma_conductor_ = 1.;
    double omega_ = 1.;

    int pc_alg_ = 0;  // preconditioner flag

    int n_iteration_ = 0;
    double system_condition_ = 0.;

    Mesh& mesh_;
    FEM_System fe_system_;

    H1_Space           O_space_;
    Hcurl_Space        T_space_;
    H1_Space        T_space_1_H1_s_;
    H1_Space        T_space_1_H1_v_;

    H1_Space           global_H1_;

    // notice that block should be unique inside block_rack, because each block has their own row_offset and col_offset,
    // using same block twice inside block_rack will destory the offset record.
    // if the same matrix is appeared multiple times inside block rack, using multiple blocks with same matrix pointer instead!

    // for real part
    Block              Kc_r_;
    Block              Mc_r_;    
    Block              Ki_r_;    
    Block              X__r_;   
    Block              Xt_r_;    

    // for imaginary part
    Block              Kc_c_;
    Block              Mc_c_;    
    Block              Ki_c_;    
    Block              X__c_;    
    Block              Xt_c_;    

    // rhs vector
    Block              Sc_r_;    
    Block              Si_r_;    
    Block              Sc_c_;    
    Block              Si_c_;  


    // preconditioning operator

    Block              pc_P_;
    Block              pc_LpM_; // -L + w*M

    Block              pc_PtP_;


    Block              pc_G_;  // alg1/3/4 use the same G, alg2 use a different one.
    Block              pc_I_;  // alg1/3/4 use the same I, alg2 use a different one.

    // global mapping (for alg2)
    Block              tmp_global_;


    
    // preconditioning algorithm 1  (decoupled)
    Block              pc_Qc_;
    Block              pc_Qi_;

    // preconditioning algorithm 2  (global)
    Block              pc_H1_;
    Block              pc_H2_;

    // preconditioning algorithm 3  (coupled)
    Block              pc_J1_;
    Block              pc_J2_;





    Dirichlet_BC bc_O_o_r_;
    Dirichlet_BC bc_O_i_r_;
    Dirichlet_BC bc_T_1_r_;

    Dirichlet_BC bc_O_o_c_;
    Dirichlet_BC bc_O_i_c_;
    Dirichlet_BC bc_T_1_c_;


    Dirichlet_BC bc_T_1_H1_s_;
    Dirichlet_BC bc_T_1_H1_v_;
    Dirichlet_BC bc_T_1_H1_v2_;

    Dirichlet_BC bc_global_H1_;



    
    Block_Rack br_system_;
    //Block_Rack pc_br_system_;


    int dim_;

    Key key_true_boundary_;                // key to 1D/2D true boundary element groups

    Key key_Omega_inner_boundary_1_;       // key to 1D/2D omega field boundary element groups


    Key key_conductor_interface_1_;        // key to 1D/2D conductor boundary element groups
    Key key_conductor_interface_layer_1_;  // key to 2D/3D conducting region element groups in contact with conductor boundarys.

    Key key_global_;
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
    T_Omega_c(Mesh& mesh, int pc_alg=0 );
    ~T_Omega_c();

    bool initialize_system();
    bool initialize_pc();

    bool assemble_system();

    bool assemble_pc();

    bool solve_system();

    bool solve_pc_system();


    scalar_t compute_L2_error();

    int get_n_iteration() const { return n_iteration_; }
    double get_system_condition() const { return system_condition_; }

    void finalize();


    void set_omega(double omega) {omega_ = omega; }

};


}