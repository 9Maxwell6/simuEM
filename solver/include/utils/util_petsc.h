#pragma once

#ifdef LOAD_PETSC
    #include <petsc.h>



namespace petsc_util 
{

    PetscErrorCode init_block_matrix(PetscInt rows, PetscInt cols, PetscInt nnz_est, Mat& mat);

    PetscErrorCode destroy_block_matrix(Mat& mat);



}




#endif