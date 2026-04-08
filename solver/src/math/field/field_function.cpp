#include "math/field/field_function.h"

using namespace simu;

template<int phy_dim, int ref_dim>
double S_Field_function::eval(const Ref_Coord& ref_coord, const Element_Data<phy_dim, ref_dim>& e_data) const
{
    Vector<phy_dim> phy_point = e_data.physical_point(ref_coord);
    return func_(phy_point, f_data);
}
INSTANTIATE_FIELD_EVAL(S_Field_function, eval)


template<int phy_dim, int ref_dim>
void V_Field_function::eval(const Ref_Coord& ref_coord, const Element_Data<phy_dim, ref_dim>& e_data, Eigen::Ref<VectorXd> value) const
{
    Vector<phy_dim> phy_point = e_data.physical_point(ref_coord);
    func_(phy_point, f_data, value);
}
INSTANTIATE_FIELD_EVAL_VEC(V_Field_function, eval)


