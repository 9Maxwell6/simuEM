#include "math/fem/integrator_l.h"

#include "math/fem/assemble_data.h"



using namespace simu;


template<typename Field, int phy_dim, int ref_dim, typename Vec_Type>
void Integrator__s__S::assemble_element_vector(Field& F, Element_Data<phy_dim, ref_dim>& e_data, Vec_Type& element_vector)
{
    constexpr int R = Vec_Type::RowsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        int order = e->get_geometry_order() + 2*test__space->get_basis_order();

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);
        
        Vector<R> basis(e_data.rows);
        
        for(const Integration_Point& i_p : i_p_list)
        {
            test__space->get_basis_s(i_p.coord, basis);

            double coeff = F.eval(i_p.coord, e_data);
            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));

            element_vector += coeff * i_p.weight * abs_det_J * basis;
        }
    }
}
INSTANTIATE_ELEMENT_VEC_TEMPLATE_ARGS(Integrator__s__S, assemble_element_vector, S_Field_fespace&)
INSTANTIATE_ELEMENT_VEC_TEMPLATE_ARGS(Integrator__s__S, assemble_element_vector, S_Field_function&)





template<typename Field, int phy_dim, int ref_dim, typename Vec_Type>
void Integrator__v__grad_S::assemble_element_vector(Field& F, Element_Data<phy_dim, ref_dim>& e_data, Vec_Type& element_vector)
{
    constexpr int R = Vec_Type::RowsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        int order = e->get_geometry_order() + 2*test__space->get_basis_order();

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);
        

        Vector<phy_dim> coeff(phy_dim);
        Matrix<R,ref_dim> grad_basis(e_data.rows, ref_dim);
        Matrix<R, phy_dim> phy_grad_basis(e_data.rows, phy_dim);

        for(const Integration_Point& i_p : i_p_list)
        {            
            test__space->get_ED_basis_v(i_p.coord, grad_basis);

            F.eval(i_p.coord, e_data, coeff);

            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));
            const Matrix<ref_dim, phy_dim>& inv_J = e_data.get_inv_J(i_p.coord);

            phy_grad_basis = grad_basis*inv_J;

            element_vector += i_p.weight * abs_det_J * phy_grad_basis * coeff;
        }
    }
}
INSTANTIATE_ELEMENT_VEC_TEMPLATE_ARGS(Integrator__v__grad_S, assemble_element_vector, V_Field_fespace&)
INSTANTIATE_ELEMENT_VEC_TEMPLATE_ARGS(Integrator__v__grad_S, assemble_element_vector, V_Field_function&)





template<typename Field, int phy_dim, int ref_dim, typename Vec_Type>
void Integrator__v__V::assemble_element_vector(Field& F, Element_Data<phy_dim, ref_dim>& e_data, Vec_Type& element_vector)
{
    constexpr int R = Vec_Type::RowsAtCompileTime;

    if constexpr  (phy_dim == ref_dim) 
    {
        std::call_once((*e_data.integrator_check)[INTEGRATOR_ID], [&]{ check_precondition(e_data.space_1, e_data.space_2); }); 
        
        const Element* e = e_data.e;
        const FEM_Space* test__space = e_data.shape_space_1;
        int order = e->get_geometry_order() + 2*test__space->get_basis_order();

        const std::vector<Integration_Point>& i_p_list = Integration::integration_rule_update(*e_data.i_r_list, e_data.b_shape, order);

        Vector<phy_dim> coeff(phy_dim);
        Matrix<R, ref_dim> basis(e_data.rows, ref_dim);
        Matrix<R, phy_dim> phy_basis(e_data.rows, phy_dim);
        for(const Integration_Point& i_p : i_p_list)
        {
            test__space->get_basis_v(i_p.coord, basis);
            
            F.eval(i_p.coord, e_data, coeff);

            double abs_det_J = std::abs(e_data.get_det_J(i_p.coord));
            const Matrix<ref_dim, phy_dim>& J_inv = e_data.get_inv_J(i_p.coord);   

            phy_basis = basis * J_inv;

            element_vector += i_p.weight * abs_det_J * phy_basis * coeff;
        }
    }
}
INSTANTIATE_ELEMENT_VEC_TEMPLATE_ARGS(Integrator__v__V, assemble_element_vector, V_Field_fespace&)
INSTANTIATE_ELEMENT_VEC_TEMPLATE_ARGS(Integrator__v__V, assemble_element_vector, V_Field_function&)