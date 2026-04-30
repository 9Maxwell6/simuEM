#pragma once

#include "math/field/field.h"

namespace simu {

class Field_constant : public Field
{

private:


    Field_constant(const Mesh& mesh);

public:

    Field_Type get_field_type() const override { return Field_Type::CONSTANT; }


};


}