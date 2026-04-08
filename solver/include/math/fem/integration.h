#pragma once

#include "utils/logger.h"

#include <vector>

namespace simu {

// integration point on reference element
struct Ref_Coord
{
    double x, y, z;
};

struct Integration_Point {
    Ref_Coord coord;
    double weight;
};

enum class Basis_Shape;

namespace Integration{

    const std::vector<Integration_Point>& integration_rule_update(std::vector<const std::vector<Integration_Point>*>& i_r_list, Basis_Shape& b_shape, int order);

    const std::vector<Integration_Point>& get_integration_points(Basis_Shape& b_shape, int order);

    const std::vector<Integration_Point>& get_integration_points_triangle(int order);

    const std::vector<Integration_Point>& get_integration_points_tetrahedron(int order);

};

}


