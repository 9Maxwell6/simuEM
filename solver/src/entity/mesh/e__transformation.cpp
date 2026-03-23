#include "entity/mesh/e__transformation.h"

using namespace simu;


void Triangle_o1::D_shape(const Element * e, const Integration_Point& i_p, Matrix<3, 2>& d_shape) 
{
    d_shape << -1, -1,
                1,  0,
                0,  1;
}

void Tetrahedron_o1::D_shape(const Element * e, const Integration_Point& i_p, Matrix<4, 3>& d_shape) {
        d_shape << -1, -1, -1,
                    1,  0,  0,
                    0,  1,  0,
                    0,  0,  1;
    }



template <int phy_dim>
void Triangle_o1::Jacobian(const Mesh* mesh, const Element * e, Matrix<phy_dim, 2>& J)
{
   static_assert(phy_dim == 2 || phy_dim == 3, "Triangle_o1::Jacobian only supports phy_dim = 2 or 3.");
   const Triangle* tri = static_cast<const Triangle*>(e);

   const size_t * node_idx = tri->get_node_idx();
   const Node& n0 = mesh->get_node(node_idx[0]);
   const Node& n1 = mesh->get_node(node_idx[1]);
   const Node& n2 = mesh->get_node(node_idx[2]);

   if constexpr (phy_dim == 2) {
      J << n1.x-n0.x,  n2.x-n0.x,
           n1.y-n0.y,  n2.y-n0.y;
   } else {
      J << n1.x-n0.x,  n2.x-n0.x,
           n1.y-n0.y,  n2.y-n0.y,
           n1.z-n0.z,  n2.z-n0.z;
   } 
}


template <int phy_dim>
void Tetrahedron_o1::Jacobian(const Mesh* mesh, const Element * e, Matrix<phy_dim, 3>& J)
{
   static_assert(phy_dim == 3, "Tetrahedron_o1::Jacobian only supports phy_dim = 3.");
   const Tetrahedron* tri = static_cast<const Tetrahedron*>(e);

   const size_t * node_idx = tri->get_node_idx();
   const Node& n0 = mesh->get_node(node_idx[0]);
   const Node& n1 = mesh->get_node(node_idx[1]);
   const Node& n2 = mesh->get_node(node_idx[2]);
   const Node& n3 = mesh->get_node(node_idx[3]);

   J << n1.x-n0.x,  n2.x-n0.x,  n3.x-n0.x,
        n1.y-n0.y,  n2.y-n0.y,  n3.y-n0.y,
        n1.z-n0.z,  n2.z-n0.z,  n3.z-n0.z;
}


template <int phy_dim>
void General_Element::Jacobian(const Mesh* mesh, const Element * e, Matrix<phy_dim, ref_dim>& J)
{
   Logger::error("General_Element::Jacobian: not available.");
}


template <typename Ref_Element, int phy_dim>
const Matrix<phy_dim, Ref_Element::ref_dim>&  Element_Transformation<Ref_Element, phy_dim>::compute_Jacobian(const Integration_Point& i_p, const Element * e, const Matrix<phy_dim, node_num>& node_matrix)
{
   Matrix<node_num, ref_dim> d_shape;
   Ref_Element::D_shape(e, i_p, d_shape); 
   J_ = node_matrix * d_shape;          
   
   return J_;
}


template <typename Ref_Element, int phy_dim>
const Matrix<phy_dim, Ref_Element::ref_dim>&  Element_Transformation<Ref_Element, phy_dim>::compute_Jacobian(const Element * e)
{   
   Ref_Element::Jacobian(mesh_, e, J_);
   return J_;
}


template <typename Ref_Element, int phy_dim>
const Matrix<Ref_Element::ref_dim, phy_dim>&  Element_Transformation<Ref_Element, phy_dim>::compute_inv_Jacobian()
{
   if constexpr (phy_dim == ref_dim)
      inv_J_ = J_.inverse();
   else
      inv_J_ = (J_.transpose() * J_).inverse() * J_.transpose();
   return inv_J_;
}


template <typename Ref_Element, int phy_dim>
double Element_Transformation<Ref_Element, phy_dim>::compute_det_Jacobian()
{
   if constexpr (phy_dim == ref_dim) {
      det_J_ = J_.determinant();
   } 
   else if constexpr (ref_dim == 1) {
      // reference segment 
      det_J_ = J_.col(0).norm();
   } 
   else if constexpr (phy_dim == 3 && ref_dim == 2) {
      // Surface in 3D
      double E = J_.col(0).squaredNorm();
      double G = J_.col(1).squaredNorm();
      double F = J_.col(0).dot(J_.col(1));
      det_J_ = std::sqrt(E * G - F * F);
   }
   else {
      det_J_ = std::sqrt((J_.transpose() * J_).determinant());
   }
   return det_J_;
}


template struct Element_Transformation<Triangle_o1, 2>;
template struct Element_Transformation<Triangle_o1, 3>;
template struct Element_Transformation<Tetrahedron_o1, 3>;

template struct Element_Transformation<General_Element, 2>;
template struct Element_Transformation<General_Element, 3>;

