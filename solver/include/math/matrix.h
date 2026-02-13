#pragma once

#include <Eigen/Dense>

namespace simu {

template <int R>
using Vector = Eigen::Matrix<double, R, 1>;

using VectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;



template <int R, int C>
using Matrix = Eigen::Matrix<double, R, C, Eigen::RowMajor>;


using MatrixXd = Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;


}