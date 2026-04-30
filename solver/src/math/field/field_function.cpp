#include "math/field/field_function.h"

using namespace simu;

double S_Field_function::eval(const Ref_Coord& ref_coord, const Element& e) const
{
    VectorXd phy_point(f_data.mesh->get_mesh_dimension());
    e.physical_point(*(f_data.mesh), ref_coord, phy_point);
    return func_(phy_point, f_data);
}


void V_Field_function::eval(const Ref_Coord& ref_coord, const Element& e, Eigen::Ref<VectorXd> value) const
{
    VectorXd phy_point(f_data.mesh->get_mesh_dimension());
    e.physical_point(*(f_data.mesh), ref_coord, phy_point);
    func_(phy_point, f_data, value);
}


