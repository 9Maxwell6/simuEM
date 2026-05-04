#pragma once

#include "utils/logger.h"

#include "math/ref_coord.h"
#include "math/data_format.h"
#include "math/fem/shape.h"


#include <vector>
#include <Eigen/Dense>

namespace simu {


enum class Space { 
    H_1, 
    H_curl,
    H_div 
};



class FEM_Space
{
protected: 
    int dim_; // space dimension
    int p_;   // polynomial order

    size_t n_node_;
    size_t n_edge_;
    size_t n_face_;
    size_t n_volume_;

    int n_dof_;
    int n_dof_per_node_;
    int n_dof_per_edge_;
    int n_dof_per_face_;
    int n_dof_per_volume_;

    FEM_Space() = default;
    FEM_Space(int dim, int p);
    virtual ~FEM_Space() {}

public:
    
    //FEM_Space(int dim, int p, 
    //            int n_node, int n_edge, int n_face, int n_volume, 
    //            int n_dof, int n_dof_per_node, int n_dof_per_edge, int n_dof_per_face, int n_dof_per_volume);
    


    virtual bool add_basis_shape(Basis_Shape g) = 0;

    


    virtual Space get_function_space() const = 0;

    virtual FEM_Space * get_basis_space(Basis_Shape s) const = 0;
    virtual const std::vector<Basis_Shape>& get_basis_shapes() const = 0;

    // Returns basis values at a point in the unit tetrahedron (e.g. implementation from Hcurl_Space)
    // For H1: Scalars. For HCurl: Vectors.
    //virtual void get_basis_v(Basis_Shape s, const Integration_Point& p, Eigen::Ref<MatrixXd> basis) const = 0;
    //virtual void get_basis_s(Basis_Shape s, const Integration_Point& p, Eigen::Ref<VectorXd> basis) const = 0;

    // Return vector proxy of Exterior Derivative of basis of the corresponding form.
    //virtual void get_ED_basis_v(Basis_Shape s, const Integration_Point& p, Eigen::Ref<MatrixXd> basis) const = 0;
    //virtual void get_ED_basis_s(Basis_Shape s, const Integration_Point& p, Eigen::Ref<VectorXd> basis) const = 0;


    int get_basis_order() const { return p_;}
    int get_dim() const { return dim_;}
    int get_n_node() const { return n_node_;}
    int get_n_edge() const { return n_edge_;}
    int get_n_face() const { return n_face_;}
    int get_n_volume() const { return n_volume_;}
    int get_n_dof() const { return n_dof_;}
    int get_n_dof_per_node() const { return n_dof_per_node_;}
    int get_n_dof_per_edge() const { return n_dof_per_edge_;}
    int get_n_dof_per_face() const { return n_dof_per_face_;}
    int get_n_dof_per_volume() const { return n_dof_per_volume_;}



    // in case Basis_Shape is known (e.g. implementation from Hcurl_tetrahedron)
    virtual void get_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {};
    virtual void get_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const {};

    virtual void get_ED_basis_s(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {};
    virtual void get_ED_basis_v(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const {};

    /**
     * @brief compute transformation matrix to correct dof direction on reference element to actual element.
     * 
     * @param node_idx list of global node index in mesh.
     * @param P transformation matrix to be filled.
     */
    virtual void dof_transformation(const size_t* node_idx, Eigen::Ref<MatrixXd> P) const = 0;

    
    virtual void get_edge_dual_basis(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {};
    virtual void get_face_dual_basis(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {};
    virtual void get_volume_dual_basis(const Ref_Coord& coord, Eigen::Ref<VectorXd> basis) const {};

    virtual void get_dof_signatures(const Ref_Coord& coord, Eigen::Ref<MatrixXd> basis) const {}
    



    /**
     * @brief Evaluates the dual basis "signatures" (kernels) for all DOFs on a specific topological entity.
     * 
     * This function provides the mathematical weights (scalars or vectors) used to define 
     * a Degree of Freedom (DOF) via an integral moment or point evaluation. For H1 spaces, 
     * this typically returns a scalar 1.0; for H(curl) or H(div), it returns the 
     * constant reference tangent or normal vector scaled by the appropriate polynomial.
     * 
     * @param[in] entity_dim dimension of the entity (0: Node, 1: Edge, 2: Face, 3: Cell).
     * @param[in] entity_idx local index of the entity within the element (e.g., edge index 0-5 for a Tet).
     * @param[in] coord The local coordinate relative to the entity's own reference space.
     *                  - Dim 1 (Edge): 1D coordinate s in [0, 1].
     *                  - Dim 2 (Face): 2D coordinates (s0, s1) on the reference triangle/square.
     *                  - Dim 3 (Cell): 3D coordinates (s0, s1, s2) on the reference tetrahedron/hexahedron.
     * @param[out] kernels A matrix where each column corresponds to a single DOF on the entity.
     *                     - Rows: Match the geometric dimension (1 for scalar H1, 2/3 for vector Hcurl/Hdiv).
     *                     - Cols: Number of DOFs associated with this specific entity.
     * 
     * @note For curved elements, these signatures represent the reference space directions; 
     *       the physical transformation (Jacobian) should be applied during integration.
     */
    virtual void get_dof_signature(int entity_dim, int entity_idx, const Ref_Coord& coord, Eigen::Ref<MatrixXd> kernels) const = 0;
    
};

}