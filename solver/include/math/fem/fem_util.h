#pragma once

#define INSTANTIATE_ELEMENT_MAT_TEMPLATE_BASE(Class, Method, ...)                                    \
    /* Matrix<3,3> */                                                                                \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, Matrix<3,3>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, Matrix<3,3>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, Matrix<3,3>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, Matrix<3,3>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, Matrix<3,3>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, Matrix<3,3>&);                       \
    /* Matrix<4,4> */                                                                                \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, Matrix<4,4>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, Matrix<4,4>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, Matrix<4,4>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, Matrix<4,4>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, Matrix<4,4>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, Matrix<4,4>&);                       \
    /* Matrix<4,6> */                                                                                \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, Matrix<4,6>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, Matrix<4,6>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, Matrix<4,6>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, Matrix<4,6>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, Matrix<4,6>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, Matrix<4,6>&);                       \
    /* Matrix<6,4> */                                                                                \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, Matrix<6,4>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, Matrix<6,4>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, Matrix<6,4>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, Matrix<6,4>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, Matrix<6,4>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, Matrix<6,4>&);                       \
    /* Matrix<6,6> */                                                                                \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, Matrix<6,6>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, Matrix<6,6>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, Matrix<6,6>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, Matrix<6,6>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, Matrix<6,6>&);                       \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, Matrix<6,6>&);                       \
    /* MatrixXd */                                                                                   \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, MatrixXd&);                          \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, MatrixXd&);                          \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, MatrixXd&);                          \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, MatrixXd&);                          \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, MatrixXd&);                          \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, MatrixXd&);

#define INSTANTIATE_ELEMENT_MAT_TEMPLATE(Class, Method)                                              \
    INSTANTIATE_ELEMENT_MAT_TEMPLATE_BASE(Class, Method,)

#define INSTANTIATE_ELEMENT_MAT_TEMPLATE_ARGS(Class, Method, ...)                                    \
    INSTANTIATE_ELEMENT_MAT_TEMPLATE_BASE(Class, Method, __VA_ARGS__,)

#define INSTANTIATE_MAT_TEMPLATE_BASE(Class, Method, ...)                                            \
    template void Class::Method(__VA_ARGS__ Matrix<3,3>&, G_Matrix&);                              \
    template void Class::Method(__VA_ARGS__ Matrix<4,4>&, G_Matrix&);                              \
    template void Class::Method(__VA_ARGS__ Matrix<4,6>&, G_Matrix&);                              \
    template void Class::Method(__VA_ARGS__ Matrix<6,4>&, G_Matrix&);                              \
    template void Class::Method(__VA_ARGS__ Matrix<6,6>&, G_Matrix&);                              \
    template void Class::Method(__VA_ARGS__ MatrixXd&,    G_Matrix&);                                            


#define INSTANTIATE_MAT_TEMPLATE(Class, Method)                                                      \
    INSTANTIATE_MAT_TEMPLATE_BASE(Class, Method,)

#define INSTANTIATE_MAT_TEMPLATE_ARGS(Class, Method, ...)                                            \
    INSTANTIATE_MAT_TEMPLATE_BASE(Class, Method, __VA_ARGS__,)

#include "entity/mesh/e_collection.h"
#include "math/fem/space_collection.h"

namespace simu {

inline Basis_Shape to_basis_shape(Geometry t)
{
    switch (t) {
        case Geometry::TETRAHEDRON: return Basis_Shape::TETRAHEDRON;
        default:
        {
            Logger::error("FEM_System::to_basis_geometry - type not supported yet.");
            throw std::invalid_argument("geometry not supported yet.");
        }
    }
}

inline Geometry to_element_geometry(Basis_Shape g)
{
    switch (g) {
        case Basis_Shape::TETRAHEDRON: return Geometry::TETRAHEDRON;
        default:
        {
            Logger::error("FEM_System::to_basis_geometry - geometry not supported yet.");
            throw std::invalid_argument("geometry not supported yet.");
        }
    }
}

}