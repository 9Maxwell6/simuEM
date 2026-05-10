#pragma once

#include "math/operator/operator.h"
#include "math/operator/integration.h"
#include "math/data_format.h"



#include <Eigen/Dense>



namespace simu {

/**
 * s: scalar coefficient
 * v: vector coefficient
 * S: scalar function
 * V: vector function
 * 
 * @note each integrator has unique id, start from 0 to SIZE-1, 
 * and each block will contain a flag array, 
 * 
 */
class Integrator : public Operator
{

protected:
    Integrator();

    // TODO： support more integration rules: point/edge/squares/cube/ etc.s

    //virtual  void assemble_element_matrix();
public:
    static constexpr int SIZE = 12;

};


}