#pragma once

#include "world/mesh/mesh.h"
#include "world/mesh/e_collection.h"

#include "math/data_format.h"
#include "math/ref_coord.h"



#include "math/field/field_util.h"
#include <mutex>



namespace simu {


enum class Field_Type 
{ 
    CONSTANT,
    FUNCTION,
    FEM_SPACE,
    NONE
};

struct Field_Data
{
    Field_Type f_type = Field_Type::NONE;
    const Mesh* mesh = nullptr;

    // use [dof_offset + 0, ... , dof_offset + dof_size-1] to get dof_idx from S/V_Field_fespace::dof_,
    // then use the dof_idx to fetch the dof value from S/V_Field_fespace::value_.
    size_t dof_size = 0;
    size_t dof_offset = 0;

    mutable double t = 0;  // time       [s ]
    mutable double f = 0;  // frequency  [Hz]
    size_t element_counter = 0;
    size_t dof_counter = 0;

};


class Field
{

protected:
    
    Field() = default;
    Field(const Mesh& mesh);

public:
    Field_Data f_data;

    void set_mesh(const Mesh& m)   { f_data.mesh = &m; }
    void set_time(double t)        { f_data.t = t; }
    void set_frequency(double f)   { f_data.f = f; }


    virtual Field_Type get_field_type() const = 0;

};


}