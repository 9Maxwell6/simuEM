#include "math/operator/operator.h"

#include "math/fem/assemble_data.h"
#include "math/fem/element_data.h"


using namespace simu;


const size_d* Operator::adjust_dof(int vdim, bool layout, size_d n_dof, size_d total_n_dof, const size_d dof_list[], size_d output[])
{
    if(vdim == 1) return dof_list;

    // local matrix always use layout 0.
    if(layout == 0){
        // [x1..xn, y1..yn, z1..zn]
        for (int c = 0; c < vdim; ++c)
            for (size_d i = 0; i < n_dof; ++i)
                output[c * n_dof + i] = dof_list[i] + c * total_n_dof;
    }else if(layout == 1){
        // [x1,y1,z1, x2,y2,z2, ...]
        for (size_d i = 0; i < n_dof; ++i)
            for (int c = 0; c < vdim; ++c)
                output[c * n_dof + i] = vdim * dof_list[i] + c;
    }
    return output;
}

/**
 * @brief Apply dof transformation to element matrix, which make sure the dof direction are aligned globally.
 * 
 * @tparam phy_dim        physical dimension of the space.
 * @tparam ref_dim        dimension of the element.
 * @tparam Mat_Type       fixed or dynamic matrix type.
 * @param e_data          element data struct, contain flag of transformation type.
 * @param element_matrix  local element matrix.
 * 
 * @example if test and trial space are both in Hcurl, transformation P * element_matrix * P.transpose() will be applied,
 * where P is transformation matrix in Hcurl space.
 * 
 */
template<int phy_dim, int ref_dim, typename Mat_Type>
void Operator::dof_transformation_mat(Element_Data<phy_dim, ref_dim>& e_data, Mat_Type& element_matrix)
{
    constexpr int R = Mat_Type::RowsAtCompileTime;
    constexpr int C = Mat_Type::ColsAtCompileTime;

    Matrix<R, R> P_1(e_data.rows, e_data.rows);
    Matrix<C, C> P_2(e_data.cols, e_data.cols);

    e_data.shape_space_1->dof_transformation(e_data.e->get_node_idx(), P_1);
    e_data.shape_space_2->dof_transformation(e_data.e->get_node_idx(), P_2);

    element_matrix = P_1 * element_matrix * P_2;
}
INSTANTIATE_ELEMENT_MAT_TEMPLATE(Operator, dof_transformation_mat)





/**
 * @brief Apply dof transformation to element vector, which make sure the dof direction are aligned globally.
 * 
 * @tparam phy_dim        physical dimension of the space.
 * @tparam ref_dim        dimension of the element.
 * @tparam Vec_Type       fixed or dynamic vector type.
 * @param e_data          element data struct, contain flag of transformation type.
 * @param element_vector  local element matrix.
 * 
 * @example if test space is not in H1, we have transformation P * element_vector 
 * 
 */
template<int phy_dim, int ref_dim, typename Vec_Type>
void Operator::dof_transformation_vec(Element_Data<phy_dim, ref_dim>& e_data, Vec_Type& element_vector)
{
    constexpr int R = Vec_Type::RowsAtCompileTime;

    Space s_1 = e_data.space_1->get_function_space();
    
    Matrix<R, R> P_1(e_data.rows, e_data.rows);

    e_data.shape_space_1->dof_transformation(e_data.e->get_node_idx(), P_1);

    element_vector = P_1 * element_vector;
}
INSTANTIATE_ELEMENT_VEC_TEMPLATE(Operator, dof_transformation_vec)





/**
 * Insert entries of a local element matrix into a global matrix.
 *
 * @tparam Mat_Type        fixed or dynamic matrix type
 * @param data             data package containing global row/column indices and pointer to global matrix
 * @param element_matrix   Local matrix to insert
 */
template<typename Mat_Type>
void Operator::add_to_global_mat(const Assemble_Data& data, Mat_Type& element_matrix)
{
    int  vdim_1   = data.space_1->get_vdim();
    bool layout_1 = data.space_1->get_layout();

    int  vdim_2   = data.space_2->get_vdim();
    bool layout_2 = data.space_2->get_layout();

    std::vector<dof_idx> row_buf(element_matrix.rows());
    std::vector<dof_idx> col_buf(element_matrix.cols());
    
    // TODO: optimize this algorithm!
    // local element matrix is always use layout 1 [x1,y1,z1, x2,y2,z2, ...]
    const size_d* row_ptr = adjust_dof(vdim_1, layout_1, element_matrix.rows()/vdim_1 , data.row_size/vdim_1, &(*data.row_dof)[data.row_dof_offset], row_buf.data());
    const size_d* col_ptr = adjust_dof(vdim_2, layout_2, element_matrix.cols()/vdim_2 , data.col_size/vdim_2, &(*data.col_dof)[data.col_dof_offset], col_buf.data());

    

    la_kernel::add_to_mat(element_matrix.rows(), row_ptr, 
                          element_matrix.cols(), col_ptr, 
                          element_matrix.data(), data.block_matrix);


    //la_kernel::add_to_mat(element_matrix.rows(), &(*data.row_dof)[data.row_dof_offset], 
    //                      element_matrix.cols(), &(*data.col_dof)[data.col_dof_offset], 
    //                      element_matrix.data(), data.block_matrix);



#ifdef LOAD_PETSC
    // G_Matrix: using petsc Mat.
    //PetscCallVoid(MatSetValues(*data.block_matrix, 
    //                            rows, &(*data.row_dof)[data.row_dof_offset], 
    //                            cols, &(*data.col_dof)[data.col_dof_offset], 
    //                            element_matrix.data(), ADD_VALUES));


#elif
    // G_Matrix: using eigen sparse matrix.
    Logger::error("Operator::add_to_global_mat - Eigen::SparseMatrix assemble not available yet.");
#endif
}
INSTANTIATE_MAT_TEMPLATE_ARGS(Operator, add_to_global_mat, const Assemble_Data&)



/**
 * Insert entries of a local element vector into a global vector.
 *
 * @tparam Vec_Type        fixed or dynamic vector type
 * @param data             data package containing global row/column indices and pointer to global matrix
 * @param element_vector   Local vertor to insert
 */
template<typename Vec_Type>
void Operator::add_to_global_vec(const Assemble_Data& data, Vec_Type& element_vector)
{
    const auto rows = element_vector.size();

    la_kernel::add_to_vec(element_vector.size(), &(*data.row_dof)[data.row_dof_offset], element_vector.data(), data.block_vector);

#ifdef LOAD_PETSC
    // G_Vector using petsc Vec.
    //PetscCallVoid(VecSetValues(*data.block_vector, 
    //                           rows, &(*data.row_dof)[data.row_dof_offset], 
    //                            element_vector.data(), ADD_VALUES));
    





#elif
    // G_Vector: using eigen vectorXd.
    Logger::error("Operator::add_to_global_vec - Eigen global vector assemble not available yet.");
#endif
}
INSTANTIATE_VEC_TEMPLATE_ARGS(Operator, add_to_global_vec, const Assemble_Data&)