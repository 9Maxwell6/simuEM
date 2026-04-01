#pragma once

#include "config.h"


#include <Eigen/Dense>

#ifdef LOAD_PETSC
    #include <petsc.h>
#endif

namespace simu {

template <int R>
using Vector = Eigen::Matrix<double, R, 1>;

using VectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;



template <int R, int C>
//using Matrix = Eigen::Matrix<double, R, C, Eigen::RowMajor>;
using Matrix = Eigen::Matrix<double, R, C>;


//using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;


// global sparse matrix
#ifdef LOAD_PETSC
    using G_Matrix = Mat;   
    using dof_idx = PetscInt;    
#else 
    using G_Matrix = std::shared_ptr<Eigen::SparseMatrix<double>>;
    using dof_idx = size_t;   
#endif

}