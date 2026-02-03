#pragma once

#include "fem_space.h"


class Hcurl_tetrahedron : public FEM_Space {
private:
    std::vector<std::vector<double>> basis_cache;

public:
    Hcurl_tetrahedron(int p = 1, const QuadratureRule& qrule=null);


}