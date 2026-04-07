#pragma once

#include "math/field/field.h"

#include <variant>

namespace simu {



template<int phy_dim>
class S_Field_function : public Field<phy_dim>
{
    std::function<double(const Vector<phy_dim>&)> func_;
public:

    template<typename Func>
    S_Field_function(Func&& f) : func_(std::forward<Func>(f)) {}

    template<int ref_dim>
    double eval(const Integration_Point& i_p, const Element_Data<phy_dim, ref_dim>& e_data) const {
        //return func_(e_data.physical_point(q));
    }
};



template<int phy_dim>
class V_Field_function : public Field<phy_dim>
{
    std::function<Vector<phy_dim>(const Vector<phy_dim>&)> func_;
public:
    template<typename Func>
    V_Field_function(Func&& f) : func_(std::forward<Func>(f)) {}

    template<int ref_dim>
    Vector<phy_dim> eval(const Integration_Point& i_p, const Element_Data<phy_dim, ref_dim>& e_data) const {
        //return func_(e_data.physical_point(q));
    }
};


}