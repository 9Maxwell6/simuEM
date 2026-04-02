#include "math/fem/operation.h"

#include "math/fem/assemble_data.h"


using namespace simu;

/**
 * @brief Apply dof transformation to element matrix, which make sure the dof direction are aligned globally.
 * 
 * @tparam phy_dim        fixed or dynamic matrix type
 * @tparam ref_dim        fixed or dynamic matrix type
 * @tparam Mat_Type        fixed or dynamic matrix type
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

    if (e_data.space_1 == Space::H_curl || e_data.space_1 == Space::H_div) {
        Matrix<R, R> P_1(e_data.rows, e_data.rows);
        if (e_data.space_1 == Space::H_curl) 
            e_data.e->compute_dof_transformation_H_curl(*e_data.mesh, P_1);
        else 
            e_data.e->compute_dof_transformation_H_div(*e_data.mesh, P_1);
        
        element_matrix = P_1 * element_matrix;
    }

    if (e_data.space_2 == Space::H_curl || e_data.space_2 == Space::H_div) {
        Matrix<C, C> P_2(e_data.cols, e_data.cols);
        if (e_data.space_2 == Space::H_curl) 
            e_data.e->compute_dof_transformation_H_curl(*e_data.mesh, P_2);
        else 
            e_data.e->compute_dof_transformation_H_div(*e_data.mesh, P_2);
        
        element_matrix = element_matrix * P_2.transpose();
    }
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Operation, dof_transformation)




/**
 * Insert entries of a local element matrix into a global matrix.
 *
 * @tparam Mat_Type        fixed or dynamic matrix type
 * @param data             data package containing global row/column indices and pointer to global matrix
 * @param element_matrix   Local matrix to insert
 */
template<typename Mat_Type>
void Operation::add_to_global(const Assemble_Data& data, Mat_Type& element_matrix)
{
    const auto rows = element_matrix.rows();
    const auto cols = element_matrix.cols();

#ifdef LOAD_PETSC
    // G_Matrix: using petsc Mat.
    PetscCallVoid(MatSetValues(*data.block_matrix, 
                            rows, &(*data.row_dof)[data.row_dof_offset], 
                            cols, &(*data.col_dof)[data.col_dof_offset], 
                            element_matrix.data(), ADD_VALUES));


#elif
    // G_Matrix: using eigen sparse matrix.
#endif
}
INSTANTIATE_MAT_TEMPLATE_ARGS(Operation, add_to_global, const Assemble_Data&)
