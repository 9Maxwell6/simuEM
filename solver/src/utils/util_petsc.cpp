#include "utils/util_petsc.h"

#ifdef LOAD_PETSC


namespace petsc_util
{

/**
 * @brief Initialize PETSc runtime.
 *
 * @param argc  pointer to command line argument count.
 * @param argv  pointer to command line argument array.
 * 
 * @return PetscErrorCode  PETSC_SUCCESS on success.
 */
PetscErrorCode petsc_initialize(int* argc, char*** argv) 
{
    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(argc, argv, nullptr, nullptr));
    PetscFunctionReturn(PETSC_SUCCESS);
}




/**
 * @brief Finalize PETSc runtime.
 *
 * @return PetscErrorCode  PETSC_SUCCESS on success.
 */
PetscErrorCode petsc_finalize() 
{
    PetscFunctionBeginUser;
    PetscCall(PetscFinalize());
    PetscFunctionReturn(PETSC_SUCCESS);
}




/**
 * Creates and preallocates a PETSc sparse matrix.
 *
 * @note Sets MAT_NEW_NONZERO_ALLOCATION_ERR so that any insertion outside
 * the preallocated sparsity pattern triggers an error.
 * 
 * TODO: currently only support sequential sparse matrix, extend it to mpi-version.
 *
 * @param rows  total number of rows.
 * @param cols  total number of columns.
 * @param nnz   number of nonzeros per row (length must equal rows).
 * @param mat   PETSc matrix to be created and preallocated.
 *
 * @return PETSc error code.
 */
PetscErrorCode petsc_init_mat(PetscInt rows, PetscInt cols, const std::vector<PetscInt>& nnz, Mat& mat)
{
    PetscFunctionBeginUser;
    PetscCall(MatCreate(PETSC_COMM_WORLD, &mat));
    PetscCall(MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols));
    PetscCall(MatSetType(mat, MATAIJ));
    PetscCall(MatSeqAIJSetPreallocation(mat, 0, nnz.data()));
    PetscCall(MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode petsc_init_vec(PetscInt size, Vec& vec)
{
    PetscFunctionBeginUser;
    PetscCall(VecCreate(PETSC_COMM_WORLD, &vec));
    PetscCall(VecSetSizes(vec, PETSC_DECIDE, size));
    PetscCall(VecSetFromOptions(vec));
    PetscCall(VecZeroEntries(vec));
    PetscFunctionReturn(PETSC_SUCCESS);
}




/**
 * @brief Set values into a PETSc matrix.
 *
 * @param row_size number of rows.
 * @param rows     row indices.
 * @param col_size number of columns.
 * @param cols     column indices.
 * @param values   values to insert.
 * @param mat      PETSc matrix to insert values into.
 *
 * @return PetscErrorCode  PETSC_SUCCESS on success.
 */
PetscErrorCode petsc_add_to_mat(PetscInt row_size, const PetscInt rows[],
                                PetscInt col_size, const PetscInt cols[],
                                const PetscScalar values[], Mat mat) 
{
    PetscFunctionBeginUser;
    PetscCall(MatSetValues(mat, row_size, rows, col_size, cols, values, ADD_VALUES));
    PetscFunctionReturn(PETSC_SUCCESS);
}




/**
 * @brief Set values into a PETSc vector.
 *
 * @param row_size number of entries.
 * @param rows     row indices.
 * @param values   values to insert.
 * @param vec      PETSc vector to insert values into.
 * 
 * @return PetscErrorCode  PETSC_SUCCESS on success.
 */
PetscErrorCode petsc_add_to_vec(PetscInt row_size, const PetscInt rows[], const PetscScalar values[], Vec vec) 
{
    PetscFunctionBeginUser;
    PetscCall(VecSetValues(vec, row_size, rows, values, ADD_VALUES));
    PetscFunctionReturn(PETSC_SUCCESS);
}




PetscErrorCode petsc_finalize_mat(Mat mat) 
{
    PetscFunctionBeginUser;
    PetscCall(MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode petsc_finalize_vec(Vec vec) 
{
    PetscFunctionBeginUser;
    PetscCall(VecAssemblyBegin(vec));
    PetscCall(VecAssemblyEnd(vec));
    PetscFunctionReturn(PETSC_SUCCESS);
}




PetscErrorCode petsc_save_mat(Mat mat, const std::string& file_name)
{
    PetscFunctionBeginUser;
    PetscViewer viewer_out;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filepath.c_str(), FILE_MODE_WRITE, &viewer_out));
    PetscCall(MatView(mat, viewer_out));
    PetscCall(PetscViewerDestroy(&viewer_out));
    PetscFunctionReturn(PETSC_SUCCESS);
}



PetscErrorCode petsc_load_mat(Mat mat, const std::string& file_name)
{
    PetscFunctionBeginUser;
    PetscViewer viewer_in;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filepath.c_str(), FILE_MODE_READ, &viewer_in));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &mat));
    PetscCall(MatSetFromOptions(mat));
    PetscCall(MatLoad(mat, viewer_in));
    PetscCall(PetscViewerDestroy(&viewer_in));
    PetscFunctionReturn(PETSC_SUCCESS);
}




PetscErrorCode petsc_save_vec(Vec vec, const std::string& file_name) 
{
    PetscFunctionBeginUser;
    PetscViewer viewer_out;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filepath.c_str(), FILE_MODE_WRITE, &viewer_out));
    PetscCall(VecView(vec, viewer_out));
    PetscCall(PetscViewerDestroy(&viewer_out));
    PetscFunctionReturn(PETSC_SUCCESS);
}




PetscErrorCode petsc_load_vec(Vec vec, const std::string& file_name) 
{
    PetscFunctionBeginUser;
    PetscViewer viewer_in;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filepath.c_str(), FILE_MODE_READ, &viewer_in));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &vec));
    PetscCall(VecLoad(vec, viewer_in));
    PetscCall(PetscViewerDestroy(&viewer_in));
    PetscFunctionReturn(PETSC_SUCCESS);
}



PetscErrorCode petsc_destroy_mat(Mat& mat) 
{
    PetscFunctionBeginUser;
    PetscCall(MatDestroy(&mat));
    PetscFunctionReturn(PETSC_SUCCESS);
}



PetscErrorCode petsc_destroy_vec(Vec& vec) 
{
    PetscFunctionBeginUser;
    PetscCall(VecDestroy(&vec));
    PetscFunctionReturn(PETSC_SUCCESS);
}



/**
 * @brief Resize a vector of PETSc Mat objects, destroying any truncated matrices.
 *
 * When shrinking, all Mat objects beyond the new size are destroyed via MatDestroy().
 * When growing, new entries are initialized to nullptr.
 * When the size is unchanged, this is a no-op.
 *
 * @param[in,out] mat_list  Vector of PETSc Mat handles to resize.
 * @param[in]     size      Desired new size of the vector.
 * @return PetscErrorCode   PETSC_SUCCESS on success, or an error code if any MatDestroy() fails.
 */
PetscErrorCode petsc_resize_mat_list(std::vector<Mat>& mat_list, size_t size)
{
    PetscFunctionBeginUser;
    for (size_t i = size; i < mat_list.size(); ++i)
    {
        if (mat_list[i]) PetscCall(MatDestroy(&mat_list[i]));
    }
    mat_list.resize(size, nullptr);
    PetscFunctionReturn(PETSC_SUCCESS);
}



/**
 * @brief Resize a vector of PETSc Vec objects, destroying any truncated vectors.
 *
 * When shrinking, all Vec objects beyond the new size are destroyed via VecDestroy().
 * When growing, new entries are initialized to nullptr.
 * When the size is unchanged, this is a no-op.
 *
 * @param[in,out] vec_list  Vector of PETSc Vec handles to resize.
 * @param[in]     size      Desired new size of the vector.
 * @return PetscErrorCode   PETSC_SUCCESS on success, or an error code if any VecDestroy() fails.
 */
PetscErrorCode petsc_resize_vec_list(std::vector<Vec>& vec_list, size_t size)
{
    PetscFunctionBeginUser;
    for (size_t i = size; i < vec_list.size(); ++i) 
    {
        if (vec_list[i]) PetscCall(VecDestroy(&vec_list[i]));
    }
    vec_list.resize(size, nullptr);
    PetscFunctionReturn(PETSC_SUCCESS);
}






// ==================================================== for debug ======================================================

PetscErrorCode petsc_print_mat(Mat mat, const std::string& name) 
{
    PetscFunctionBeginUser;
    if (!name.empty()) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "--- %s ---\n", name.c_str()));
    PetscCall(MatView(mat, PETSC_VIEWER_STDOUT_WORLD));
    return PETSC_SUCCESS;
}

PetscErrorCode petsc_print_vec(Vec vec, const std::string& name) 
{
    PetscFunctionBeginUser;
    if (!name.empty()) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "--- %s ---\n", name.c_str()));
    PetscCall(VecView(vec, PETSC_VIEWER_STDOUT_WORLD));
    return PETSC_SUCCESS;
}


PetscErrorCode petsc_save_ascii_mat(Mat mat, const std::string& file_name) 
{
    PetscFunctionBeginUser;
    PetscViewer viewer_out;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filepath.c_str(), &viewer_out));
    PetscCall(MatView(mat, viewer_out));
    PetscCall(PetscViewerDestroy(&viewer_out));
    return PETSC_SUCCESS;
}


PetscErrorCode petsc_save_ascii_vec(Vec vec, const std::string& file_name) 
{
    PetscFunctionBeginUser;
    PetscViewer viewer_out;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filepath.c_str(), &viewer_out));
    PetscCall(VecView(vec, viewer_out));
    PetscCall(PetscViewerDestroy(&viewer_out));
    return PETSC_SUCCESS;
}


}

#endif