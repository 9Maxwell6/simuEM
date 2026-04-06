#pragma once

#include "config.h"


#include <Eigen/Dense>
#include <Eigen/Sparse>


#ifdef LOAD_PETSC
    #include <petsc.h>
#endif

namespace simu {

template <int R>
using Vector = Eigen::Matrix<double, R, 1>;

using VectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;

using S_Field_3d = std::function<double(Vector<3>)>;
using S_Field_2d = std::function<double(Vector<2>)>;
using S_Field_1d = std::function<double(Vector<1>)>;
using V_Field_3d = std::function<Vector<3>(Vector<3>)>;
using V_Field_2d = std::function<Vector<2>(Vector<2>)>;
using V_Field_1d = std::function<Vector<1>(Vector<1>)>;


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
#else 
    // local element matrix: column major 
    template<int R, int C>
    using Matrix = Eigen::Matrix<double, R, C>;
    using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

    using G_Matrix = std::shared_ptr<Eigen::SparseMatrix<double>>;
    using G_Vector = std::shared_ptr<Eigen::Matrix<double, Eigen::Dynamic, 1>>;
    using dof_idx = size_t;   
    using size_d = int;
#endif

}