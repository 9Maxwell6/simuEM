#pragma once

#ifdef LOAD_PETSC
    #include <petsc.h>
    #include <petscsys.h>


#include "math/data_format.h"

#include <vector>
#include <iostream>
namespace petsc_util {

PetscErrorCode petsc_initialize(int* argc, char*** argv);
PetscErrorCode petsc_finalize();


PetscErrorCode petsc_init_mat(PetscInt row_size, PetscInt col_size, const std::vector<PetscInt>& nnz, Mat& mat);

PetscErrorCode petsc_init_vec(PetscInt size, Vec& vec);

PetscErrorCode petsc_create_transpose(const Mat mat_A, Mat &mat_B);

PetscErrorCode petsc_create_nest_mat(PetscInt b_row_size, PetscInt b_col_size, const std::vector<Mat> &block_mat, Mat& mat);
PetscErrorCode petsc_create_nest_vec(const std::vector<Vec> &block_vec, Vec &vec);

PetscErrorCode petsc_get_local_size_mat(Mat mat, PetscInt* n_row, PetscInt* n_col);

PetscErrorCode petsc_extract_block_mat(const std::vector<PetscInt>& row_local_sizes, const std::vector<PetscInt>& col_local_sizes, Mat mat, std::vector<Mat>& blocks);
PetscErrorCode petsc_extract_block_vec(const std::vector<PetscInt>& row_local_sizes, Vec vec, std::vector<Vec>& blocks);


PetscErrorCode petsc_add_to_mat(PetscInt row_size, const PetscInt rows[],
                                PetscInt col_size, const PetscInt cols[],
                                const PetscScalar values[], Mat mat);

PetscErrorCode petsc_add_to_vec(PetscInt row_size, const PetscInt rows[], const PetscScalar values[], Vec vec);

PetscErrorCode petsc_zero_row_col_mat(const std::vector<PetscInt>& dofs, PetscScalar diag_val, Mat mat, Vec x, Vec b);


PetscErrorCode petsc_set_value_vec(const std::vector<PetscInt>& dofs, const std::vector<PetscScalar>& values, Vec vec);

PetscErrorCode petsc_extract_vec(Vec vec, std::vector<simu::scalar_t> &out);

PetscErrorCode petsc_finalize_mat(Mat mat);
PetscErrorCode petsc_finalize_vec(Vec vec);

PetscErrorCode petsc_is_ready_mat(const Mat mat, PetscBool *ready);
PetscErrorCode petsc_is_ready_vec(const Vec vec, PetscBool *ready);




PetscErrorCode petsc_destroy_mat(Mat& mat);
PetscErrorCode petsc_destroy_vec(Vec& vec);


PetscErrorCode petsc_save_mat(Mat mat, const std::string& file_name);
PetscErrorCode petsc_load_mat(Mat mat, const std::string& file_name);
PetscErrorCode petsc_save_vec(Vec vec, const std::string& file_name);
PetscErrorCode petsc_load_vec(Vec vec, const std::string& file_name);


// resize petsc vec: need destory un-used petsc mat or vec
PetscErrorCode petsc_resize_mat_list(std::vector<Mat>& mat_list, size_t size);
PetscErrorCode petsc_resize_vec_list(std::vector<Vec>& vec_list, size_t size);


PetscErrorCode petsc_ksp_convergence(KSP ksp, bool *successful, int *iterations, const std::string &label = "KSP");

// for debug

PetscErrorCode petsc_print_mat(Mat mat, const std::string& name = "");
PetscErrorCode petsc_print_vec(Mat mat, const std::string& name = "");
PetscErrorCode petsc_save_ascii_mat(Mat mat, const std::string& file_name);
PetscErrorCode petsc_save_ascii_vec(Vec vec, const std::string& file_name);

PetscErrorCode petsc_print_memory_usage(const char *label);

}




#endif