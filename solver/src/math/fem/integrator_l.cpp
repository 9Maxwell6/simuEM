#include "math/fem/integrator_l.h"

#include "math/fem/assemble_data.h"



using namespace simu;


template<typename Field, int phy_dim, int ref_dim, typename Vec_Type>
void Integrator__s__S::assemble_element_vector(Field&& F, Element_Data<phy_dim, ref_dim>& e_data, Vec_Type& element_matrix)
{
    constexpr int R = Vec_Type::RowsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        const FEM_Space* trial_space = e_data.shape_space_2;
        int order = e->get_geometry_order() + trial_space->get_basis_order() + test__space->get_basis_order();

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);
        
        Vector<R> basis(e_data.rows);

        for(const Integration_Point& i_p : i_p_list)
        {
            trial_space->get_basis_s(i_p, basis);
            
            double abs_det_J = std::abs(e_data.get_det_J(i_p));

            //element_matrix += coeff * i_p.weight * abs_det_J * basis * basis.transpose();
        }
    }
}
INSTANTIATE_ELEMENT_VEC_TEMPLATE_ARGS(Integrator__s__S, assemble_element_vector, S_Field_3d&&)
INSTANTIATE_ELEMENT_VEC_TEMPLATE_ARGS(Integrator__s__S, assemble_element_vector, S_Field_2d&&)