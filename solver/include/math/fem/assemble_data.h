#pragma once
#include "entity/mesh/e_collection.h"
#include "math/fem/space_collection.h"

namespace simu {


struct Assemble_Data {
    const FEM_Space* space_1;
    const FEM_Space* space_2;

    const std::vector<Element*>* elements;

    size_t row_size;
    size_t col_size;

    const std::vector<size_t>* row_dof;
    const std::vector<size_t>* col_dof;
    
    // ... whatever assemble needs

    
};

}
