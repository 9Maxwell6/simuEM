#pragma once

#ifdef LOAD_PETSC
    #include "utils/util_petsc.h"
#endif

#include <vector>

#include "math/data_format.h"


namespace la_kernel {

using namespace simu;


void initialize(int* argc, char*** argv);
void finalize();

void init_mat(size_t row_size, size_t col_size, std::vector<size_d> nnz, G_Matrix& mat);
void init_vec(size_t row_size, G_Vector& vec);


void destroy_mat(G_Matrix& mat);
void destroy_vec(G_Vector& vec);

void add_to_mat(size_d row_size, const size_d rows[], size_d col_size, const size_d cols[], const scalar_t values[], G_Matrix mat);
void add_to_vec(size_d row_size, const size_d rows[], const scalar_t values[], G_Vector vec);

void resize_mat_list(std::vector<G_Matrix>& mat_list, size_t size);
void resize_vec_list(std::vector<G_Vector>& vec_list, size_t size);



}