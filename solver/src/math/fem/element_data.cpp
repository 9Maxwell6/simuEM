#include "math/fem/element_data.h"

using namespace simu;

template<int phy_dim, int ref_dim>
const Matrix<phy_dim, ref_dim>& Element_Data<phy_dim, ref_dim>::get_J(const Ref_Coord& ref_coord)
{
    if(flag_J) return J;

    flag_J = e->compute_Jacobian(*mesh, ref_coord, J);

    return J;
}

template<int phy_dim, int ref_dim>
const Matrix<ref_dim, phy_dim>& Element_Data<phy_dim, ref_dim>::get_inv_J(const Ref_Coord& ref_coord)
{
    if(flag_inv_J) return inv_J;

    if(!flag_J) flag_J = e->compute_Jacobian(*mesh, ref_coord, J);
    flag_inv_J = flag_J;
    
    if constexpr (phy_dim == ref_dim)
        inv_J = J.inverse();
    else
        inv_J = (J.transpose() * J).inverse() * J.transpose();
    return inv_J;
}






template<int phy_dim, int ref_dim>
double Element_Data<phy_dim, ref_dim>::get_det_J(const Ref_Coord& ref_coord)
{
    if(flag_det_J) return det_J;

    if(!flag_J) flag_J = e->compute_Jacobian(*mesh, ref_coord, J);
    flag_det_J = flag_J;


    if constexpr (phy_dim == ref_dim) {
        det_J = J.determinant();
    } 
    else if constexpr (ref_dim == 1) {
        // reference segment 
        det_J = J.col(0).norm();
    } 
    else if constexpr (phy_dim == 3 && ref_dim == 2) {
        // Surface in 3D
        double E = J.col(0).squaredNorm();
        double G = J.col(1).squaredNorm();
        double F = J.col(0).dot(J.col(1));
        det_J = std::sqrt(E * G - F * F);
    }
    else {
        det_J = std::sqrt((J.transpose() * J).determinant());
    }
    return det_J;
}


template<int phy_dim, int ref_dim>
Vector<phy_dim> Element_Data<phy_dim, ref_dim>::physical_point(const Ref_Coord& ref_coord) const
{

    Vector<phy_dim> phy_coord = Vector<phy_dim>::Zero();

    int node_num = e->get_node_num();
    const size_t * node_idx =  e->get_node_idx();

    int geo_node_num = e->get_geometry_node_num();
    // TODO: higher geometry order not support yet, mesh object should provide extra nodes for higher order elements.
    // possible idea: using Bounding volume hierarchy (BVH) algorithm to quickly locate current element in mesh. then get all the higher order nodes.

    VectorXd shape(geo_node_num);

    e->compute_shape(ref_coord, shape);

    for(int i=0; i<node_num; ++i)
    {
        const Node& n = mesh->get_node(node_idx[i]);
        if      constexpr (phy_dim == 3){ phy_coord += shape[i] * Vector<3>{n.x, n.y, n.z}; }
        else if constexpr (phy_dim == 2){ phy_coord += shape[i] * Vector<2>{n.x, n.y};      }
        else                            { phy_coord += shape[i] * Vector<1>{n.x};           }
    }

    for(int j=node_num; j<geo_node_num; ++j ){
        Logger::error("Element_Data::physical_point - higher geometry order not support yet");
    }

    return phy_coord;

    
}



template struct Element_Data<3, 3>;
template struct Element_Data<3, 2>;
template struct Element_Data<3, 1>;
template struct Element_Data<2, 2>;
template struct Element_Data<2, 1>;
template struct Element_Data<1, 1>;
