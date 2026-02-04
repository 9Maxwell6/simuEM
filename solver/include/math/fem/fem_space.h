#pragma once

#include "integration.h"

#include <vector>
#include <Eigen/Dense>

enum class Space { 
    H1, 
    HCurl 
};

enum class Geometry {
    Tetrahedron, 
    // currently not support other element geometry
};


class FEM_Space
{
private: 
    int dim_;
    int p_;  // polynomial order


public:
    virtual ~FEM_Space() {}

    // Returns the number of DOFs per element
    virtual int get_element_dof() const = 0;
    virtual int get_order() const = 0;

    // Returns basis values at a point in the unit tetrahedron
    // For H1: Scalars. For HCurl: Vectors.
    virtual void get_basis(Integration_Point p, const Eigen::Ref<Eigen::MatrixXd> basis) const = 0;



    virtual Space get_function_space() const = 0;
};