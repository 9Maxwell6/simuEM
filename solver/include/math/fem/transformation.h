#pragma once

#include "world/mesh/e_collection.h"
#include "math/fem/space_collection.h"

#include "math/fem/block_rack.h"
#include "math/fem/integration.h"
#include "math/matrix.h"

#include "utils/logger.h"


#include <Eigen/Dense>

namespace simu {



template <int phy_dim, int ref_dim, int node_num>
struct Element_Transformation
{

protected: 
    Matrix<phy_dim, ref_dim> J_;
    Matrix<ref_dim, phy_dim> inv_J_;
    double det_J_;

    const Integration_Point * i_p_;

    bool reset_ = true;
    

    void get_D_shape(Geometry geometry, int order, Matrix<node_num, ref_dim>& d_shape);


    void compute_Jacobian(const Integration_Point& i_p, const Matrix<phy_dim, node_num>& node_matrix, Matrix<phy_dim, ref_dim>& J);



    // TODO: implement general case of compute_Jacobian
    /*
    template <int phy_dim, int ref_dim, int node_dim>
    void compute_Jacobian(const Integration_Point& i_p, const Matrix<phy_dim, node_dim>& node_matrix, Matrix<phy_dim, ref_dim>& J)
    {
        // reference shape derivatives: (same as exterior derivative in H1/Lagrange-space)
        Matrix<node_dim, ref_dim> d_shape;

        calc_geometry_derivative<ref_dim, node_dim>(i_p, d_shape);

        // 2. Compute J = node_matrix * dshape
        // This is now a multiplication of two fixed-size matrices.
        // For p=1 Tet, this is (3x4) * (4x3) = (3x3).
        J = node_matrix * dshape;
    }

    void calc_geometry_derivative(Geometry g, int p, const Integration_Point& i_p, Matrix<node_dim, phy_dim>& d_shape);

    */

    void compute_inverse_Jacobian(Matrix<phy_dim, ref_dim>& J, Matrix<ref_dim, phy_dim>& inv_J);
    
};




using Transformation_Triangle_2D_o1        = Element_Transformation<2,2,3>;  
using Transformation_Triangle_3D_o1        = Element_Transformation<3,2,3>; 
using Transformation_Tetreahedron_3D_o1    = Element_Transformation<3,3,4>; 
using Transformation_Tetreahedron_3D_o2    = Element_Transformation<3,3,10>; 

using Transformation_general  = Element_Transformation<Eigen::Dynamic, Eigen::Dynamic, Eigen::Dynamic>; 








namespace Transformation {

// explicitly for polynomial order-1
void compute_Jacobian_Triangle_2D_p1(const Node& n1, const Node& n2, const Node& n3, Matrix<2,2>& J);
void compute_Jacobian_Triangle_3D_p1(const Node& n1, const Node& n2, const Node& n3, Matrix<3,2>& J);
void compute_Jacobian_Tetrahedron_3D_p1(const Node& n1, const Node& n2, const Node& n3, const Node& n4, Matrix<3,3>& J);

}



}