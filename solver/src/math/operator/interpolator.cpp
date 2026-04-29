#include "math/operator/interpolator.h"

#include "math/fem/element_data.h"


using namespace simu;

template<int phy_dim, int ref_dim, typename Mat_Type>
void Interpolator_H1_to_Hcurl::interpolate_element(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once(e_data.check->interpolator_check[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        const FEM_Space* trial_space = e_data.shape_space_2;

        
        
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Interpolator_H1_to_Hcurl, interpolate_element)