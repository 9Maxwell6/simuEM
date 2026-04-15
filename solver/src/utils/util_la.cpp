#include "utils/util_la.h"


namespace la_kernel
{

using namespace simu;


void initialize(int* argc, char*** argv)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_initialize(argc, argv);
    #endif
}


void finalize()
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_finalize();
    #endif
}


void init_mat(size_t row_size, size_t col_size, std::vector<size_d> nnz, G_Matrix& mat)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_init_mat(row_size, col_size, nnz, mat);
    #else
        mat = std::make_shared<Eigen::SparseMatrix<scalar_t>>(row_size, col_size);
        Eigen::VectorXi nnz_eigen = Eigen::Map<Eigen::VectorXi>(nnz.data(), nnz.size());
        mat->reserve(nnz_eigen);
    #endif
}



void init_vec(size_t row_size, G_Vector& vec)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_init_vec(row_size, vec);
    #else
        vec = std::make_shared<Eigen::VectorXd>(Eigen::VectorXd::Zero(row_size));
    #endif
}


void create_transpose(const G_Matrix mat_A, G_Matrix &mat_B)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_create_transpose(mat_A, mat_B);  // no copy
    #else
        mat_B = std::make_shared<Eigen::SparseMatrix<scalar_t>>(mat_A->transpose());  // copy
    #endif
}



void create_nest_mat(size_d b_row_size, size_d b_col_size, const std::vector<G_Matrix> &block_mat, G_Matrix& mat)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_create_nest_mat(b_row_size, b_col_size, block_mat, mat);
    #else
        // TODO: implement with eigen library.
        Logger::error("la_kernel::create_nest_mat: default implementation not ready, only petsc version available.");
    #endif
}



void create_nest_vec(const std::vector<G_Vector>& block_vec, G_Vector &vec)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_create_nest_vec(block_vec, vec);
    #else
        // TODO: implement with eigen library.
        Logger::error("la_kernel::create_nest_vec: default implementation not ready, only petsc version available.");
    #endif
    
}



void destroy_mat(G_Matrix& mat)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_destroy_mat(mat);
    #else
        mat.reset();
    #endif
}



void destroy_vec(G_Vector& vec)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_destroy_vec(vec);
    #else
        vec.reset();
    #endif
}




void add_to_mat(size_d row_size, const size_d rows[], size_d col_size, const size_d cols[], const scalar_t values[], G_Matrix mat)
{
    #ifdef LOAD_PETSC
        // G_Matrix: using petsc Mat.
        petsc_util::petsc_add_to_mat(row_size, rows, col_size, cols, values, mat);
    #else
        // G_Matrix: using eigen sparse matrix.
        for (int i = 0; i < row_size; ++i)
            for (int j = 0; j < col_size; ++j)
                mat->coeffRef(rows[i], cols[j]) += values[i*col_size + j];
    #endif
}




void add_to_vec(size_d row_size, const size_d rows[], const scalar_t values[], G_Vector vec)
{
    #ifdef LOAD_PETSC
        // G_Vector: using petsc Vec.
        petsc_util::petsc_add_to_vec(row_size, rows, values, vec);
    #else
        // G_Vector: using eigen VectorXd.
        for (int i = 0; i < row_size; ++i) (*vec)(rows[i]) += values[i];
    #endif
}




void finalize_mat(G_Matrix mat)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_finalize_mat(mat);
    #else
        mat->makeCompressed();      
    #endif
}



void finalize_vec(G_Vector vec)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_finalize_vec(vec);
    #endif
}



bool is_ready_mat(const G_Matrix mat)
{
    #ifdef LOAD_PETSC
        // strong check
        PetscBool ready;
        petsc_util::petsc_is_ready_mat(mat, &ready);
        return (ready == PETSC_TRUE);
    #else
        // strong check
        return mat && mat->rows()>0 && mat->cols()>0 && mat->isCompressed();
    #endif

}


bool is_ready_vec(const G_Vector vec)
{
    #ifdef LOAD_PETSC
        // weak check
        PetscBool ready;
        petsc_util::petsc_is_ready_vec(vec, &ready);
        return (ready == PETSC_TRUE);
    #else
        // weak check
        return vec && vec->size() > 0;
    #endif
}



void resize_mat_list(std::vector<G_Matrix>& mat_list, size_t size)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_resize_mat_list(mat_list, size);
    #else
        mat_list.resize(size);
    #endif
}




void resize_vec_list(std::vector<G_Vector>& vec_list, size_t size)
{
    #ifdef LOAD_PETSC
        petsc_util::petsc_resize_vec_list(vec_list, size);
    #else
        vec_list.resize(size);
    #endif
}



}
