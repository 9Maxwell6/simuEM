#pragma once

#ifdef LOAD_PETSC
    #include "utils/util_petsc.h"
#endif


#include "utils/logger.h"
#include "math/data_format.h"

#include <vector>


/**
 * @brief linear algebra kernel interface. 
 * 
 * @note functions will pick operations from the specified library based on the solver configureation.
 * 
 * Currently support PETSc, Eigen(not complete).
 */
namespace la_kernel {

using namespace simu;


void initialize(int* argc, char*** argv);
void finalize();

void init_mat(size_t row_size, size_t col_size, std::vector<size_d> nnz, G_Matrix& mat);
void init_vec(size_t row_size, G_Vector& vec);

void create_transpose(const G_Matrix mat_A, G_Matrix &mat_B);

void create_nest_mat(size_d b_row_size, size_d b_col_size, const std::vector<G_Matrix> &block_mat, G_Matrix& mat);
void create_nest_vec(const std::vector<G_Vector>& block_vec, G_Vector &vec);


void destroy_mat(G_Matrix& mat);
void destroy_vec(G_Vector& vec);

void add_to_mat(size_d row_size, const size_d rows[], size_d col_size, const size_d cols[], const scalar_t values[], G_Matrix mat);
void add_to_vec(size_d row_size, const size_d rows[], const scalar_t values[], G_Vector vec);

void finalize_mat(G_Matrix mat);
void finalize_vec(G_Vector vec);

bool is_ready_mat(const G_Matrix mat);
bool is_ready_vec(const G_Vector vec);


void resize_mat_list(std::vector<G_Matrix>& mat_list, size_t size);
void resize_vec_list(std::vector<G_Vector>& vec_list, size_t size);



}