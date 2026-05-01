#pragma once
#include "world/mesh/e_collection.h"
#include "math/fem/space_collection.h"
#include "world/mesh/mesh.h"
#include "math/operator/operator_collection.h"
#include "math/fem/block.h"




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

    // to check if the constraints of each operator applied is satisfied.
    // e.g. FEM_Space space_1 must be Hcurl space.
    mutable Check_Flag check;


    
};




}
