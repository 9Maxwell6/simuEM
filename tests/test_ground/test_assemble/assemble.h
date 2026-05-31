#pragma once

#include "world/mesh/mesh.h"
#include "world/structure/structure.h"

#include "math/fem/fem_space.h"
#include "math/fem/fem_system.h"
#include "math/fem/assemble_mat.h"
#include "math/fem/assemble_vec.h"

#include "math/fem/bc_dirichlet.h"


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

class Assemble
{

private:
    double system_condition_ = 0.;

    Mesh& mesh_;
    FEM_System fe_system_;

    Block              dof_field_;

    Dirichlet_BC bc_;


    // preconditioner
    H1_Space           H1_v_;
    H1_Space           H1_s_;
    Hcurl_Space        Hcurl_;
 




    int dim_;

    Key key_true_boundary_;                // key to 1D/2D true boundary element groups
    Key key_interior_;

    Block_Rack br_system_;


public:
    Assemble(Mesh& mesh);

    bool assemble_system();


};


}