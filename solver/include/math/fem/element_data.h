#pragma once
#include "world/mesh/e_collection.h"
#include "math/fem/space_collection.h"
#include "math/operator/operator_collection.h"
#include "world/mesh/mesh.h"
#include "math/fem/block.h"
#include "math/dof/dof_manager.h"



namespace simu {



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

    // [n_node, n_edge, n_face, n_cell]  used for initialize dof_manager
    const std::array<size_t,4>* entity_size;

    const Element* e;

    // to check if the constraints of each operator applied is satisfied.
    // e.g. FEM_Space space_1 must be Hcurl space.
    Check_Flag* check;

    const FEM_Space* shape_space_1;  // trial space of specific geometry shape
    const FEM_Space* shape_space_2;  // test  space of specific geometry shape

    const FEM_Space* space_1;
    const FEM_Space* space_2;

    mutable DoF_Manager* dof_manager;

    Basis_Shape b_shape;             // geometry shape


    // ---------- used for post process only ------------
    const std::vector<const FEM_Space*>* space_list; // function space each block belong to.
    const std::vector<std::vector<scalar_t>>* dof_value_list; // value of each dof on this element for each block.

    // --------------------------------------------------


    Matrix<phy_dim, ref_dim>      J;
    Matrix<ref_dim, phy_dim>  inv_J;
    double                    det_J;

    size_t rows;
    size_t cols;

    // flag for whether the Jacobian data can be reused.
    bool flag_J     = false;
    bool flag_inv_J = false;
    bool flag_det_J = false;

    const Matrix<phy_dim, ref_dim>& get_J(const Ref_Coord& ref_coord);
    const Matrix<ref_dim, phy_dim>& get_inv_J(const Ref_Coord& ref_coord);
    double                          get_det_J(const Ref_Coord& ref_coord);

    Vector<phy_dim> physical_point(const Ref_Coord& ref_coord) const;

    void reset_flag() 
    { 
        flag_J = false;
        flag_inv_J = false;
        flag_det_J = false;
    }
};







}
