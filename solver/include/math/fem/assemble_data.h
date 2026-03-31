#pragma once
#include "entity/mesh/e_collection.h"
#include "math/fem/space_collection.h"
#include "entity/mesh/mesh.h"
#include "math/fem/integrator.h"

#include <mutex>

namespace simu {


struct Assemble_Data
{
    int mesh_dim;
    int element_dim;

    size_t row_size;
    size_t col_size;

    const Mesh* mesh;

    const FEM_Space* space_1;  // trial space
    const FEM_Space* space_2;  // test  space

    const std::vector<Element*>* elements;

    const std::vector<size_t>* row_dof;
    const std::vector<size_t>* col_dof;

    std::unordered_map<Basis_Shape , std::vector<const std::vector<Integration_Point>*>, Shape_Hash>& integration_rule;

    mutable std::array<std::once_flag, Integrator::SIZE> integrator_check_flags;
    
};





/**
 * Element_Data - per-element data shared across integrators.
 *
 * Available members:
 *   e             - const Element*, current element
 *   J             - Matrix<phy_dim, ref_dim>, Jacobian at current quad point
 *   inv_J         - Matrix<ref_dim, phy_dim>, inverse/pseudo-inverse of Jacobian
 *   det_J         - double, determinant of Jacobian
 *   b_shape       - Basis_Shape, element geometry
 *   shape_space_1 - const FEM_Space*, trial function space
 *   shape_space_2 - const FEM_Space*, test function space
 *   i_r_list      - integration rules per order
 *
 *
 * Template parameters:
 *   phy_dim - physical space dimension (1, 2, or 3)
 *   ref_dim - reference element dimension (1, 2, or 3), ref_dim <= phy_dim
 *
 * Note: accessed via auto& in user callbacks. Use e_data.e, e_data.J, etc.
 */
template<int phy_dim, int ref_dim>
struct Element_Data
{
    const Mesh* mesh;

    const Element* e;

    std::vector<const std::vector<Integration_Point>*>* i_r_list;
    std::array<std::once_flag, Integrator::SIZE>* integrator_check;

    const FEM_Space* shape_space_1;  // trial space of specific geometry shape
    const FEM_Space* shape_space_2;  // test  space of specific geometry shape

    // flag for transformation from conventional dof direction on reference element to dof direction on actual element.
    Space space_1;
    Space space_2;

    Basis_Shape b_shape;             // geometry shape

    Matrix<phy_dim, ref_dim>      J;
    Matrix<ref_dim, phy_dim>  inv_J;
    double                    det_J;

    size_t rows;
    size_t cols;

    // flag for whether the Jacobian data can be reused.
    bool flag_J     = false;
    bool flag_inv_J = false;
    bool flag_det_J = false;

    const Matrix<phy_dim, ref_dim>& get_J(const Integration_Point& i_p);
    const Matrix<ref_dim, phy_dim>& get_inv_J(const Integration_Point& i_p);
    double                          get_det_J(const Integration_Point& i_p);

    void reset_flag() 
    { 
        flag_J = false;
        flag_inv_J = false;
        flag_det_J = false;
    }
};

}
