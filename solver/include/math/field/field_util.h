#pragma once

// For scalar eval: double Class::Method(const Ref_Coord&, const Element_Data<P,R>&) const
#define INSTANTIATE_FIELD_EVAL(Class, Method)                                                               \
    template double Class::Method(const Ref_Coord&, const Element_Data<3,3>&) const;                        \
    template double Class::Method(const Ref_Coord&, const Element_Data<3,2>&) const;                        \
    template double Class::Method(const Ref_Coord&, const Element_Data<3,1>&) const;                        \
    template double Class::Method(const Ref_Coord&, const Element_Data<2,2>&) const;                        \
    template double Class::Method(const Ref_Coord&, const Element_Data<2,1>&) const;                        \
    template double Class::Method(const Ref_Coord&, const Element_Data<1,1>&) const;


// For vector eval: void Class::Method(const Ref_Coord&, const Element_Data<P,R>&, Ref<VectorXd>) const
#define INSTANTIATE_FIELD_EVAL_VEC(Class, Method)                                                           \
    template void Class::Method(const Ref_Coord&, const Element_Data<3,3>&, Eigen::Ref<VectorXd>) const;              \
    template void Class::Method(const Ref_Coord&, const Element_Data<3,2>&, Eigen::Ref<VectorXd>) const;              \
    template void Class::Method(const Ref_Coord&, const Element_Data<3,1>&, Eigen::Ref<VectorXd>) const;              \
    template void Class::Method(const Ref_Coord&, const Element_Data<2,2>&, Eigen::Ref<VectorXd>) const;              \
    template void Class::Method(const Ref_Coord&, const Element_Data<2,1>&, Eigen::Ref<VectorXd>) const;              \
    template void Class::Method(const Ref_Coord&, const Element_Data<1,1>&, Eigen::Ref<VectorXd>) const;


