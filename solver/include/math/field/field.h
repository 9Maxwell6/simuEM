#pragma once

#include "entity/mesh/mesh.h"
#include "entity/mesh/e_collection.h"
#include "math/fem/assemble_data.h"
#include "math/data_format.h"


namespace simu {


enum class Field_Type 
{ 
    CONSTANT,
    FUNCTION,
    FEM_SPACE
};

template<int phy_dim>
class Field
{

protected:
    const Mesh* mesh_ = nullptr;

    double t_ = 0;  // time       [s ]
    double f_ = 0;  // frequency  [Hz]

    Field(const Mesh& mesh);

public:
    

};


}