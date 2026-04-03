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

    

    PetscErrorCode destroy_petsc_matrix(Mat& mat)
    {
        PetscCall(MatDestroy(&mat));
    }


}

#endif