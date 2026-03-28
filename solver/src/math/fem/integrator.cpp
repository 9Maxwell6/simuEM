#include "math/fem/integrator.h"



using namespace simu;

template<int phy_dim, int ref_dim, typename Mat_Type>
void Integrator__s_S__S::assemble_element_matrix(double coeff, Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;


    if constexpr  (phy_dim == ref_dim && R == C) 
    {
        const Element* e = e_data.e;
        const FEM_Space* trial_space = e_data.shape_space_1;
        const FEM_Space* test__space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order() + test__space->get_basis_order();

        if(e_data.i_r_list->size() < order) e_data.i_r_list->resize(order+1, nullptr);
        
        const std::vector<Integration_Point>* i_p_list = (*e_data.i_r_list)[order];
        if(!i_p_list) {
            i_p_list = &Integration::get_integration_points(e_data.b_shape, order);
            (*e_data.i_r_list)[order] = i_p_list;
        }

        
        
        for(int i=0; i<i_p_list->size(); ++i)
        {
            const Integration_Point& i_p = (*i_p_list)[i];
            Vector<R> basis;
            trial_space->get_basis_s(i_p, basis);
            
            double det_J = e_data.get_det_Jacobian(i_p);
            
            element_matrix += coeff * i_p.weight * det_J * basis * basis.transpose();

        }

    }
}

INSTANTIATE_INTEGRATOR_TEMPLATE(Integrator__s_S__S, assemble_element_matrix, double)


