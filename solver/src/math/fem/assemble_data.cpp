#include "math/fem/assemble_data.h"

using namespace simu;


Basis_Shape to_basis_shape(Geometry t)
{
    switch (t) {
        case Geometry::TETRAHEDRON: return Basis_Shape::TETRAHEDRON;
        default:
        {
            Logger::error("FEM_System::to_basis_geometry - type not supported yet.");
            throw std::invalid_argument("geometry not supported yet.");
        }
    }
}

Geometry to_element_geometry(Basis_Shape g)
{
    switch (g) {
        case Basis_Shape::TETRAHEDRON: return Geometry::TETRAHEDRON;
        default:
        {
            Logger::error("FEM_System::to_basis_geometry - geometry not supported yet.");
            throw std::invalid_argument("geometry not supported yet.");
        }
    }
}