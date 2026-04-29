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
    template void Class::Method(__VA_ARGS__ Matrix<3,3>&);                                           \
    template void Class::Method(__VA_ARGS__ Matrix<4,4>&);                                           \
    template void Class::Method(__VA_ARGS__ Matrix<4,6>&);                                           \
    template void Class::Method(__VA_ARGS__ Matrix<6,4>&);                                           \
    template void Class::Method(__VA_ARGS__ Matrix<6,6>&);                                           \
    template void Class::Method(__VA_ARGS__ MatrixXd&   );                                            


#define INSTANTIATE_MAT_TEMPLATE(Class, Method)                                                      \
    INSTANTIATE_MAT_TEMPLATE_BASE(Class, Method,)

#define INSTANTIATE_MAT_TEMPLATE_ARGS(Class, Method, ...)                                            \
    INSTANTIATE_MAT_TEMPLATE_BASE(Class, Method, __VA_ARGS__,)


#define INSTANTIATE_ELEMENT_VEC_TEMPLATE_BASE(Class, Method, ...)                                    \
    /* Vector<2> */                                                                                  \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, Vector<2>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, Vector<2>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, Vector<2>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, Vector<2>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, Vector<2>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, Vector<2>&);                         \
    /* Vector<3> */                                                                                  \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, Vector<3>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, Vector<3>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, Vector<3>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, Vector<3>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, Vector<3>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, Vector<3>&);                         \
    /* Vector<4> */                                                                                  \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, Vector<4>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, Vector<4>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, Vector<4>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, Vector<4>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, Vector<4>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, Vector<4>&);                         \
    /* Vector<6> */                                                                                  \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, Vector<6>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, Vector<6>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, Vector<6>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, Vector<6>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, Vector<6>&);                         \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, Vector<6>&);                         \
    /* VectorXd */                                                                                   \
    template void Class::Method(__VA_ARGS__ Element_Data<3,3>&, VectorXd&);                          \
    template void Class::Method(__VA_ARGS__ Element_Data<3,2>&, VectorXd&);                          \
    template void Class::Method(__VA_ARGS__ Element_Data<3,1>&, VectorXd&);                          \
    template void Class::Method(__VA_ARGS__ Element_Data<2,2>&, VectorXd&);                          \
    template void Class::Method(__VA_ARGS__ Element_Data<2,1>&, VectorXd&);                          \
    template void Class::Method(__VA_ARGS__ Element_Data<1,1>&, VectorXd&);

#define INSTANTIATE_ELEMENT_VEC_TEMPLATE(Class, Method)                                              \
    INSTANTIATE_ELEMENT_VEC_TEMPLATE_BASE(Class, Method,)

#define INSTANTIATE_ELEMENT_VEC_TEMPLATE_ARGS(Class, Method, ...)                                    \
    INSTANTIATE_ELEMENT_VEC_TEMPLATE_BASE(Class, Method, __VA_ARGS__,)


#define INSTANTIATE_VEC_TEMPLATE_BASE(Class, Method, ...)                                            \
    template void Class::Method(__VA_ARGS__ Vector<2>&);                                             \
    template void Class::Method(__VA_ARGS__ Vector<3>&);                                             \
    template void Class::Method(__VA_ARGS__ Vector<4>&);                                             \
    template void Class::Method(__VA_ARGS__ Vector<6>&);                                             \
    template void Class::Method(__VA_ARGS__ VectorXd& );                                            


#define INSTANTIATE_VEC_TEMPLATE(Class, Method)                                                      \
    INSTANTIATE_VEC_TEMPLATE_BASE(Class, Method,)

#define INSTANTIATE_VEC_TEMPLATE_ARGS(Class, Method, ...)                                            \
    INSTANTIATE_VEC_TEMPLATE_BASE(Class, Method, __VA_ARGS__,)
