#include "math/fem/transformation.h"


using namespace simu;


namespace Transformation{

    void compute_Jacobian_Triangle_2D(const Node& n1, const Node& n2, const Node& n3, Matrix<2,2>& J)
    {
        J << n2.x - n1.x, n3.x - n1.x,
             n2.y - n1.y, n3.y - n1.y;
    }

    void compute_Jacobian_Triangle_3D(const Node& n1, const Node& n2, const Node& n3, Matrix<3,2>& J)
    {
        J << n2.x - n1.x, n3.x - n1.x,
             n2.y - n1.y, n3.y - n1.y,
             n2.z - n1.z, n3.z - n1.z;
    }


    void compute_Jacobian_Tetrahedron_3D(const Node& n1, const Node& n2, const Node& n3, const Node& n4, Matrix<3,3>& J)
    {
        J << n2.x - n1.x,  n3.x-n1.x, n4.x-n1.x,
             n2.y - n1.y,  n3.y-n1.y, n4.y-n1.y,
             n2.z - n1.z,  n3.z-n1.z, n4.z-n1.z;
    }

    

    template <int phy_dim, int ref_dim>
    void compute_inverse_Jacobian(Matrix<phy_dim, ref_dim>& J, Matrix<ref_dim, phy_dim>& inv_J)
    {
        if constexpr (phy_dim == ref_dim)
            inv_J = J.inverse();
        else
            inv_J = (J.transpose() * J).inverse() * J.transpose();
    }


}
