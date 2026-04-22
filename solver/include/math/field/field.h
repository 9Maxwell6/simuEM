#pragma once

#include "world/mesh/mesh.h"
#include "world/mesh/e_collection.h"
#include "math/fem/assemble_data.h"
#include "math/fem/element_data.h"

#include "math/data_format.h"
#include "math/ref_coord.h"

#include "math/fem/assemble_data.h"


#include "math/field/field_util.h"



namespace simu {


enum class Field_Type 
{ 
    CONSTANT,
    FUNCTION,
    FEM_SPACE
};

struct Field_Data
{
    const Mesh* mesh = nullptr;
    mutable double t = 0;  // time       [s ]
    mutable double f = 0;  // frequency  [Hz]
    size_t element_counter = 0;
    size_t dof_counter = 0;

};

class Field
{

protected:
    Field_Data f_data;


    Field() = default;
    Field(const Mesh& mesh) { f_data.mesh = &mesh; }

public:
    void set_mesh(const Mesh& m)   { f_data.mesh = &m; }
    void set_time(double t)        { f_data.t = t; }
    void set_frequency(double f)   { f_data.f = f; }

};


}