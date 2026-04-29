#pragma once

#include "math/field/field.h"

#include "math/fem/element_data.h"

#include <variant>

namespace simu {


class S_Field_function : public Field
{
    std::function<double(const Eigen::Ref<VectorXd>, const Field_Data&)> func_;
public:

    template<typename Func>
    S_Field_function(Func&& f) : func_(std::forward<Func>(f)) {}

    template<int phy_dim, int ref_dim>
    double eval(const Ref_Coord& ref_coord, const Element_Data<phy_dim, ref_dim>& e_data) const;
};



class V_Field_function : public Field
{
    std::function<void(const Eigen::Ref<VectorXd>, const Field_Data&, Eigen::Ref<VectorXd>)> func_;
public:
    template<typename Func>
    V_Field_function(Func&& f) : func_(std::forward<Func>(f)) {}

    template<int phy_dim, int ref_dim>
    void eval(const Ref_Coord& ref_coord, const Element_Data<phy_dim, ref_dim>& e_data, Eigen::Ref<VectorXd> value) const;
};


}