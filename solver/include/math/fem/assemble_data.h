#pragma once
#include "world/mesh/e_collection.h"
#include "math/fem/space_collection.h"
#include "world/mesh/mesh.h"
#include "math/operator/integrator.h"
#include "math/fem/block.h"


#include <mutex>

namespace simu {


struct Assemble_Data
{
    int mesh_dim;
    int element_dim;

    size_t row_size;
    size_t col_size;

    mutable size_t row_dof_offset;
    mutable size_t col_dof_offset;

    const Mesh* mesh;

    const FEM_Space* space_1;  // test  space
    const FEM_Space* space_2;  // trial space

    const std::vector<Element*>* elements;

    const std::vector<dof_idx>* row_dof;
    const std::vector<dof_idx>* col_dof;

    mutable G_Matrix block_matrix;   // shared pointer to global block matrix
    mutable G_Vector block_vector;   // shared pointer to global block vector

    std::unordered_map<Basis_Shape , std::vector<const std::vector<Integration_Point>*>, Shape_Hash>& integration_rule;

    mutable std::array<std::once_flag, Integrator::SIZE> integrator_check_flags;


    
};




}
