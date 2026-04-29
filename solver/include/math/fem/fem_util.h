#pragma once

    
#include "world/mesh/e_collection.h"
#include "math/fem/space_collection.h"

namespace simu {

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