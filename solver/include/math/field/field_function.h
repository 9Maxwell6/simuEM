#pragma once

#include "math/field/field.h"


namespace simu {


class S_Field_function : public Field
{
    std::function<double(const Eigen::Ref<VectorXd>, const Field_Data&)> func_;
public:

    template<typename Func>
    S_Field_function(const Mesh& mesh, Func&& f) : Field(mesh), func_(std::forward<Func>(f)) {f_data.f_type=Field_Type::FUNCTION;}

    double eval(const Ref_Coord& ref_coord, const Element& e) const;

    Field_Type get_field_type() const override { return Field_Type::FUNCTION; }
};



class V_Field_function : public Field
{
    std::function<void(const Eigen::Ref<VectorXd>, const Field_Data&, Eigen::Ref<VectorXd>)> func_;
public:
    template<typename Func>
    V_Field_function(const Mesh& mesh, Func&& f) : Field(mesh), func_(std::forward<Func>(f)) {f_data.f_type=Field_Type::FUNCTION;}

    void eval(const Ref_Coord& ref_coord, const Element& e, Eigen::Ref<VectorXd> value) const;

    Field_Type get_field_type() const override { return Field_Type::FUNCTION; }
};


}