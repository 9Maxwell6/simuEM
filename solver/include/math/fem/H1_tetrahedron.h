#pragma once

#include "fem_space.h"

class H1_tetrahedron : public FEM_Space {
private:
    std::vector<std::vector<double>> basis_cache;

public:
    H1_tetrahedron(int p = 1, const QuadratureRule& qrule=null);


}