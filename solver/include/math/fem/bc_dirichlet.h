#pragma once

#include "math/fem/space_collection.h"
#include "entity/mesh/e_collection.h"

#include "math/fem/block.h"
#include "math/data_format.h"

#include "utils/logger.h"
#include "utils/util_la.h"

#include <vector>

namespace simu {

enum class Dirichlet_Type
{
    HOMOGENEOUS,  // default
    CONSTANT,
    FIELD
};

struct Dirichlet_BC
{
    Dirichlet_Type bc_type; 

    const Block* block;
    const FEM_Space* fe_space;
    const std::vector<Element*>* bdr_elements;

    std::vector<dof_idx> bc_dofs;

    // size 0 for HOMOGENEOUS, 1 for CONSTANT, dofs.size() for FIELD
    std::vector<scalar_t> bc_values;

    // store all dof per elements: [element_1 dofs, element_2 dofs, ...], only needed for DirichletKind::FIELD
    // can be cleared after all values are computed.
    std::vector<dof_idx> bc_element_dof;

    double get_bc(size_t i) const 
    {
        switch (bc_type) {
            case Dirichlet_Type::HOMOGENEOUS: return 0.0;
            case Dirichlet_Type::CONSTANT:    return bc_values[0];
            case Dirichlet_Type::FIELD:       return bc_values[i];
        }
        return 0.;  // HOMOGENEOUS by default
    }

    // TODO: for now we only support constant BC,
    // e.g. for Field type in Hcurl, need to compute the surface linearform integral of (g⋅t) * N
    // where (g⋅t) is tangential component of user defined vector field, N is Hcurl basis function at surface.
    // void evaluate();

};


}
