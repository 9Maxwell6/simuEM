#pragma once

#ifdef LOAD_PETSC
    #include <petsc.h>

#include <vector>

namespace petsc_util 
{

PetscErrorCode init_petsc_matrix(PetscInt rows, PetscInt cols, const std::vector<PetscInt>& nnz, Mat& mat);

PetscErrorCode init_petsc_vector(PetscInt size, Vec& vec);


PetscErrorCode finalize_matrix(Mat& mat);
PetscErrorCode finalize_vector(Vec& vec);


PetscErrorCode destroy_petsc_matrix(Mat& mat);

PetscErrorCode save_matrix(Mat& mat, const std::string& file_name);
PetscErrorCode load_matrix(Mat& mat, const std::string& file_name);
PetscErrorCode save_vector(Vec& vec, const std::string& file_name);
PetscErrorCode load_vector(Vec& vec, const std::string& file_name);



    // for debug

PetscErrorCode print_matrix(Mat& mat, const std::string& name = "");
PetscErrorCode print_vector(Mat& mat, const std::string& name = "");
PetscErrorCode save_ascii_mat(Mat& mat, const std::string& file_name);
PetscErrorCode save_ascii_vec(Mat& mat, const std::string& file_name);


}




#endif