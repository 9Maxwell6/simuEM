#include "utils/util_petsc.h"

#ifdef LOAD_PETSC


namespace petsc_util{

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
PetscErrorCode init_petsc_matrix(PetscInt rows, PetscInt cols, const std::vector<PetscInt>& nnz, Mat& mat)
{
    PetscFunctionBeginUser;
    PetscCall(MatCreate(PETSC_COMM_WORLD, &mat));
    PetscCall(MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols));
    PetscCall(MatSetType(mat, MATAIJ));
    PetscCall(MatSeqAIJSetPreallocation(mat, 0, nnz.data()));
    PetscCall(MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode init_petsc_vector(PetscInt size, Vec& vec)
{
    PetscFunctionBeginUser;
    PetscCall(VecCreate(PETSC_COMM_WORLD, &vec));
    PetscCall(VecSetSizes(vec, PETSC_DECIDE, size));
    PetscCall(VecSetFromOptions(vec));
    PetscCall(VecZeroEntries(vec));
    PetscFunctionReturn(PETSC_SUCCESS);
}





PetscErrorCode finalize_matrix(Mat& mat) 
{
    PetscCall(MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY));
    return PETSC_SUCCESS;
}

PetscErrorCode finalize_vector(Vec& vec) 
{
    PetscCall(VecAssemblyBegin(vec));
    PetscCall(VecAssemblyEnd(vec));
    return PETSC_SUCCESS;
}




PetscErrorCode save_matrix(Mat& mat, const std::string& file_name)
{
    PetscViewer viewer_out;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filepath.c_str(), FILE_MODE_WRITE, &viewer_out));
    PetscCall(MatView(mat, viewer_out));
    PetscCall(PetscViewerDestroy(&viewer_out));
    return PETSC_SUCCESS;
}



PetscErrorCode load_matrix(Mat& mat, const std::string& file_name) {
    PetscViewer viewer_in;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filepath.c_str(), FILE_MODE_READ, &viewer_in));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &mat));
    PetscCall(MatSetFromOptions(mat));
    PetscCall(MatLoad(mat, viewer_in));
    PetscCall(PetscViewerDestroy(&viewer_in));
    return PETSC_SUCCESS;
}




PetscErrorCode save_vector(Vec& vec, const std::string& file_name) 
{
    PetscViewer viewer_out;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filepath.c_str(), FILE_MODE_WRITE, &viewer_out));
    PetscCall(VecView(vec, viewer_out));
    PetscCall(PetscViewerDestroy(&viewer_out));
    return PETSC_SUCCESS;
}




PetscErrorCode load_vector(Vec& vec, const std::string& file_name) 
{
    PetscViewer viewer_in;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filepath.c_str(), FILE_MODE_READ, &viewer_in));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &vec));
    PetscCall(VecLoad(vec, viewer_in));
    PetscCall(PetscViewerDestroy(&viewer_in));
    return PETSC_SUCCESS;
}




PetscErrorCode destroy_petsc_matrix(Mat& mat)
{
    PetscCall(MatDestroy(&mat));
}






// ==================================================== for debug ======================================================

PetscErrorCode print_matrix(Mat& mat, const std::string& name) 
{
    if (!name.empty()) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "--- %s ---\n", name.c_str()));
    PetscCall(MatView(mat, PETSC_VIEWER_STDOUT_WORLD));
    return PETSC_SUCCESS;
}

PetscErrorCode print_vec(Vec& vec, const std::string& name) 
{
    if (!name.empty()) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "--- %s ---\n", name.c_str()));
    PetscCall(VecView(vec, PETSC_VIEWER_STDOUT_WORLD));
    return PETSC_SUCCESS;
}


PetscErrorCode save_ascii_mat(Mat& mat, const std::string& file_name) 
{
    PetscViewer viewer_out;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filepath.c_str(), &viewer_out));
    PetscCall(MatView(mat, viewer_out));
    PetscCall(PetscViewerDestroy(&viewer_out));
    return PETSC_SUCCESS;
}


PetscErrorCode save_ascii_vec(Vec& vec, const std::string& file_name) 
{
    PetscViewer viewer_out;
    std::string filepath = std::string(PETSC_DATA_OUTPUT_DIR) + "/" + file_name;

    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filepath.c_str(), &viewer_out));
    PetscCall(VecView(vec, viewer_out));
    PetscCall(PetscViewerDestroy(&viewer_out));
    return PETSC_SUCCESS;
}


}

#endif