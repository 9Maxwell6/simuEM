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

    int n_node_;
    int n_edge_;
    int n_face_;
    int n_volume_;

    int n_dof_;
    int n_dof_per_node_;
    int n_dof_per_edge_;
    int n_dof_per_face_;
    int n_dof_per_volume_;
    

    

public:
    FEM_Space(int dim, int p);
    //FEM_Space(int dim, int p, 
    //            int n_node, int n_edge, int n_face, int n_volume, 
    //            int n_dof, int n_dof_per_node, int n_dof_per_edge, int n_dof_per_face, int n_dof_per_volume);
    
    virtual ~FEM_Space() {}

    const int dim_; // space dimension
    const int p_;   // polynomial order

    

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
    int get_n_node() const { return n_node_;}
    int get_n_edge() const { return n_edge_;}
    int get_n_face() const { return n_face_;}
    int get_n_volume() const { return n_volume_;}
    int get_n_dof() const { return n_dof_;}
    int get_n_dof_per_node() const { return n_dof_per_node_;}
    int get_n_dof_per_edge() const { return n_dof_per_edge_;}
    int get_n_dof_per_face() const { return n_dof_per_face_;}
    int get_n_dof_per_volume() const { return n_dof_per_volume_;}
};

}