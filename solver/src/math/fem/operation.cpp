#include "math/fem/operation.h"

using namespace simu;

/**
 * @brief Apply dof transformation to element matrix, which make sure the dof direction are aligned globally.
 * 
 * 
 * @param e_data element data struct, contain flag of transformation type.
 * @param element_matrix local element matrix.
 * 
 * @example if test and trial space are both in Hcurl, transformation P * element_matrix * P.transpose() will be applied,
 * where P is transformation matrix in Hcurl space.
 * 
 */
template<int phy_dim, int ref_dim, typename Mat_Type>
void Operation::dof_transformation(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;
    const auto rows = element_matrix.rows();
    const auto cols = element_matrix.cols();

    if (e_data.flag_H_curl_space_1 || e_data.flag_H_div__space_1) {
        Matrix<R, R> P_1;
        P_1.setZero(rows, rows);
        if (e_data.flag_H_curl_space_1) 
            e_data.e->compute_dof_transformation_H_curl(*e_data.mesh, P_1);
        else 
            e_data.e->compute_dof_transformation_H_div(*e_data.mesh, P_1);
        
        element_matrix = P_1 * element_matrix;
    }

    if (e_data.flag_H_curl_space_2 || e_data.flag_H_div__space_2) {
        Matrix<C, C> P_2;
        P_2.setZero(cols, cols);
        MatrixXd P_trial;
        if (e_data.flag_H_curl_space_2) 
            e_data.e->compute_dof_transformation_H_curl(*e_data.mesh, P_2);
        else 
            e_data.e->compute_dof_transformation_H_div(*e_data.mesh, P_2);
        
        element_matrix = element_matrix * P_2.transpose();
    }
}
INSTANTIATE_OPERATION_TEMPLATE(Operation, dof_transformation)
