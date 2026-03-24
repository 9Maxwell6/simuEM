#pragma once



#include "entity/mesh/e_collection.h"
#include "math/fem/space_collection.h"
#include "math/matrix.h"

#include <Eigen/Dense>


namespace simu {

/*
/**
 * contains all information for local Galerkin matrix assemble.
 */

struct Element_Info
{
private:
    friend class FEM_System;

    const Element* element_;

    const FEM_Space * space_1_;
    const FEM_Space * space_2_;



    Eigen::Ref<MatrixXd> J_;      
    Eigen::Ref<MatrixXd> inv_J_;
    double det_J_;
};

class Operation
{
protected:
    Operation();
};


}

