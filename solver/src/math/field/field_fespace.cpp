#include "math/field/field_fespace.h"


using namespace simu;

template<int phy_dim, int ref_dim>
double S_Field_fespace::eval(const Ref_Coord& ref_coord, const Element_Data<phy_dim, ref_dim>& e_data) const
{
    Logger::error("S_Field_fespace::eval - not implemented.");
    return 0;
}
INSTANTIATE_FIELD_EVAL(S_Field_fespace, eval)


template<int phy_dim, int ref_dim>
void V_Field_fespace::eval(const Ref_Coord& ref_coord, const Element_Data<phy_dim, ref_dim>& e_data, Eigen::Ref<VectorXd> value) const
{
    Logger::error("S_Field_fespace::eval - not implemented.");
}
INSTANTIATE_FIELD_EVAL_VEC(V_Field_fespace, eval)
