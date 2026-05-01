#pragma once

    
#include "world/mesh/e_collection.h"

namespace simu {

enum class Basis_Shape {
    EDGE,
    TRIANGLE,
    TETRAHEDRON, 
    // currently not support other element geometry
};

constexpr size_t N_SHAPE = 3;
static_assert(static_cast<size_t>(Basis_Shape::TETRAHEDRON) + 1 == N_SHAPE, "N_SHAPES is out of sync with Basis_Shape enum");


struct Shape_Hash {
    size_t operator()(Basis_Shape s) const {
        return std::hash<int>{}(static_cast<int>(s));
    }
};


inline Basis_Shape to_basis_shape(Geometry t)
{
    switch (t) {
        case Geometry::TRIANGLE: return Basis_Shape::TRIANGLE;
        case Geometry::TETRAHEDRON: return Basis_Shape::TETRAHEDRON;
        default:
        {
            Logger::error("FEM_System::to_basis_shape - type not supported yet.");
            throw std::invalid_argument("geometry not supported yet.");
        }
    }
}

inline Geometry to_element_geometry(Basis_Shape g)
{
    switch (g) {
        case Basis_Shape::TRIANGLE: return Geometry::TRIANGLE;
        case Basis_Shape::TETRAHEDRON: return Geometry::TETRAHEDRON;
        default:
        {
            Logger::error("FEM_System::to_basis_geometry - geometry not supported yet.");
            throw std::invalid_argument("geometry not supported yet.");
        }
    }
}

}