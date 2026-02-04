#pragma once

#include "fem_space.h"

class H1_tetrahedron : public FEM_Space {
private:
    std::vector<std::vector<double>> basis_cache;

public:
    H1_tetrahedron(int p = 1);

    ~H1_tetrahedron(){};

    void get_basis(Integration_Point p, const Eigen::Ref<Eigen::MatrixXd> basis) const override {
        const double * ss = basis.data();
    };
};