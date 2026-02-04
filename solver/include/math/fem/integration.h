#pragma once

#include <vector>

// integration point on reference element
struct Integration_Point {
    double x, y, z;
    double weight;
};


namespace Integration{

    const std::vector<Integration_Point>& get_integrationPoints_triangle(int order);


    const std::vector<Integration_Point>& get_integrationPoints_tetrahedron(int order);

};


