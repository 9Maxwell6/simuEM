#include "utils/util_petsc.h"

#ifdef LOAD_PETSC


namespace petsc_util{

    PetscErrorCode init_block_matrix(PetscInt rows, PetscInt cols, PetscInt nz, Mat& mat)
    {
        PetscFunctionBeginUser;
        PetscCall(MatCreate(PETSC_COMM_WORLD, &mat));
        PetscCall(MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols));
        PetscCall(MatSetType(mat, MATAIJ));
        PetscCall(MatSeqAIJSetPreallocation(mat, nz, NULL));
        PetscFunctionReturn(PETSC_SUCCESS);
    }

    

    PetscErrorCode destroy_block_matrix(Mat& mat)
    {
        PetscCall(MatDestroy(&mat));
    }


}

#endif