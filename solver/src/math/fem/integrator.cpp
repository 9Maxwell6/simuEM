#include "math/fem/integrator.h"

#include <stdexcept>


using namespace simu;

template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_S__S::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    if constexpr  (phy_dim == ref_dim && R == C) 
    {
        
    }
}

INSTANTIATE_INTEGRATOR_TEMPLATE(Integrator__s_S__S, assemble_element_matrix, double)


