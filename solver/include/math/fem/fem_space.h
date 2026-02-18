#pragma once

#include "integration.h"
#include "../matrix.h"

#include <vector>
#include <Eigen/Dense>

namespace simu {

enum class Space { 
    H_1, 
    H_curl 
};

enum class Geometry {
    TETRAHEDRON, 
    // currently not support other element geometry
};


class FEM_Space
{
protected: 
    int dim_; // space dimension
    int p_;   // polynomial order

    int n_dof_;
    int n_dof_node_offset_;
    int n_dof_edge_offset_;
    int n_dof_face_offset_;
    int n_dof_volume_offset_;

public:
    FEM_Space(int p);
    virtual ~FEM_Space() {}

    // Returns the number of DOFs per element
    virtual int get_element_dof() const = 0;
    virtual Space get_function_space() const = 0;

    // Returns basis values at a point in the unit tetrahedron
    // For H1: Scalars. For HCurl: Vectors.
    virtual void get_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const = 0;
    virtual void get_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const = 0;

    // Return vector proxy of Exterior Derivative of basis of the corresponding form.
    virtual void get_ED_basis_v(Integration_Point p, Eigen::Ref<MatrixXd> basis) const = 0;
    virtual void get_ED_basis_s(Integration_Point p, Eigen::Ref<VectorXd> basis) const = 0;


    

    int get_order() const { return p_;}
    int get_dim() const { return dim_;}
};

}