#pragma once



#include "entity/mesh/e_collection.h"
#include "math/fem/space_collection.h"
#include "math/fem/fem_util.h"

#include "math/matrix.h"



namespace simu {

template<int phy_dim, int ref_dim>
struct Element_Data;


class Operation
{
protected:
    Operation();

public:
    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static dof_transformation(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);
};


}

