#pragma once

// integration point of reference element
struct Integration_Point {
    double x_, y_, z_;
    double weight;
}

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
    int order_;


public:
    virtual ~FEM_Space() {}

    // Returns the number of DOFs per element
    virtual int get_element_dof() const = 0;
    virtual int get_order() const = 0;

    // Returns basis values at a point in the unit tetrahedron
    // For H1: Scalars. For HCurl: Vectors.
    virtual void get_basis(int id, double* basis) const = 0;



    virtual Space get_function_space() const = 0;
}