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
    INTERIOR = 1
};

class CurlCurl
{

private:
    double system_condition_ = 0.;

    Mesh& mesh_;
    FEM_System fe_system_;

    Hcurl_Space            space_;
    Block              dof_field_;

    Dirichlet_BC bc_;


    // preconditioner
    H1_Space           T_H1_v_;
    H1_Space           T_H1_s_;
    Block              pc_P_;
    Block              pc_G_;
    Block              pc_L_;   
    Block              pc_Q_;  

    Dirichlet_BC bc_v_;
    Dirichlet_BC bc_s_;


    Block_Rack br_system_;


    int dim_;

    Key key_true_boundary_;                // key to 1D/2D true boundary element groups
    Key key_interior_;



public:
    CurlCurl(Mesh& mesh);

    bool assemble_system();

    bool assemble_pc_system();

    int solve_system();

    int solve_pc_system();

    scalar_t compute_L2_error();

    double get_condition_number() {return system_condition_; }

};


}