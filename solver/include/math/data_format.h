#pragma once

#include "config.h"

#include <complex>

#include <Eigen/Dense>
#include <Eigen/Sparse>


#ifdef LOAD_PETSC
    #include <petsc.h>
#endif

namespace simu {

template <int R>
using Vector = Eigen::Matrix<double, R, 1>;

using VectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;



#ifdef LOAD_PETSC
    // local element matrix: row major 
    template<int R, int C>
    using Matrix = Eigen::Matrix<double, R, C, (C==1) ? Eigen::ColMajor : Eigen::RowMajor>;
    using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    // global matrix
    using G_Matrix = Mat;   
    using G_Vector = Vec;
    using dof_idx = PetscInt; 
    using size_d = PetscInt;

    using scalar_t  = PetscScalar;

    #if defined(PETSC_USE_COMPLEX)
        
        using real_t  = PetscScalar;  // std::complex<double>
        using cplx_t  = PetscScalar;  // std::complex<double>
    #else
        using real_t  = PetscScalar;  // double
        using cplx_t  = std::complex<PetscReal>;
    #endif
#else 
    // local element matrix: column major 
    template<int R, int C>
    using Matrix = Eigen::Matrix<double, R, C>;
    using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

    using G_Matrix = std::shared_ptr<Eigen::SparseMatrix<double>>;
    using G_Vector = std::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, 1>>;
    using dof_idx = size_t;   
    using size_d = int;

    using scalar_t  = double;

    using real_t  = double;
    using cplx_t  = std::complex<double>;
#endif

}