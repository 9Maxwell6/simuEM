#pragma once

#include "entity/mesh/e_collection.h"
#include "math/fem/block_rack.h"
#include "math/fem/integration.h"
#include "math/matrix.h"

#include <Eigen/Dense>

namespace simu {



template <int phy_dim, int ref_dim, bool linear_mapping>
struct Element_Transformation
{
    Matrix<phy_dim, ref_dim> J_;
    Matrix<ref_dim, phy_dim> inv_J_;
    double det_J_;

    const Integration_Point * i_p_;

    // Jacobian is constant for linear mapping (e.g. on triangle, tetrahedron)
    


    
};


using Transformation_Triangle_2D       = Element_Transformation<2,2, true>;  

using Transformation_Triangle_3D       = Element_Transformation<3,2, true>; 
using Transformation_Tetreahedron_3D   = Element_Transformation<3,3, true>; 





namespace Transformation{

    void compute_Jacobian_Triangle_2D(const Node& n1, const Node& n2, const Node& n3, Matrix<2,2>& J);
    void compute_Jacobian_Triangle_3D(const Node& n1, const Node& n2, const Node& n3, Matrix<3,2>& J);

    void compute_Jacobian_Tetrahedron_3D(const Node& n1, const Node& n2, const Node& n3, const Node& n4, Matrix<3,3>& J);



    template <int phy_dim, int ref_dim>
    void compute_inverse_Jacobian(Matrix<phy_dim, ref_dim>& J, Matrix<ref_dim, phy_dim>& inv_J);

}


}