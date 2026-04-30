#include "math/field/field.h"


using namespace simu;


Field::Field(const Mesh& mesh) 
{ 
    f_data.mesh = &mesh; 
}