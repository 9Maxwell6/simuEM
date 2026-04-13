#pragma once



#include "entity/mesh/e_collection.h"
#include "math/fem/space_collection.h"
#include "math/fem/fem_util.h"


#include "math/data_format.h"

#include "utils/util_la.h"




namespace simu {

struct Assemble_Data;

template<int phy_dim, int ref_dim>
struct Element_Data;


class Operation
{
protected:
    Operation();

public:
    template<int phy_dim, int ref_dim, typename Mat_Type>
    void static dof_transformation_mat(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix);

    template<int phy_dim, int ref_dim, typename Vec_Type>
    void static dof_transformation_vec(Element_Data<phy_dim, ref_dim>& e_data, Vec_Type& element_vector);

    template<typename Mat_Type>
    void static add_to_global_mat(const Assemble_Data& data, Mat_Type& element_matrix);

    template<typename Vec_Type>
    void static add_to_global_vec(const Assemble_Data& data, Vec_Type& element_vector);
};


}

