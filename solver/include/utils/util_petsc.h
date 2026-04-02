#pragma once

#ifdef LOAD_PETSC
    #include <petsc.h>

#include <vector>

namespace petsc_util 
{

    PetscErrorCode init_petsc_matrix(PetscInt rows, PetscInt cols, const std::vector<PetscInt>& nnz, Mat& mat);

    PetscErrorCode destroy_petsc_matrix(Mat& mat);



}




#endif