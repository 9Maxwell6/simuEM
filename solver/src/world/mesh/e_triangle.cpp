#include "world/mesh/mesh.h"
#include "world/mesh/e_triangle.h"

using namespace simu;

Geometry Triangle::get_geometry() const { return Geometry::TRIANGLE; }

Triangle::Triangle(const size_t *node_idx, size_t id, size_t property_id, int o) : Element(id, property_id, o)
{
   property_id_ = property_id;
   for (int i = 0; i < 3; i++)
   {
      node_idx_[i] = node_idx[i];
   }
}


Triangle::Triangle(std::vector<std::size_t> node_idx, size_t id, size_t property_id, int o) : Element(id, property_id, o)
{
   property_id_ = property_id;
   for (int i = 0; i < 3; i++)
   {
      node_idx_[i] = node_idx[i];
   }
}


Triangle::Triangle(size_t idx_1, size_t idx_2, size_t idx_3, size_t id, size_t property_id, int o) : Element(id, property_id, o)
{
    property_id_  = property_id;
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
    node_idx_[2] = idx_3;
}

void Triangle::set_node_idx(const size_t *idx)
{
   for (int i = 0; i < 3; i++)
   {
      node_idx_[i] = idx[i];
   }
}

void Triangle::set_node_idx(size_t idx_1, size_t idx_2, size_t idx_3)
{
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
    node_idx_[2] = idx_3;
}

const size_t * Triangle::get_node_idx() const
{
   return node_idx_;
}


bool Triangle::compute_Jacobian(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> J) const 
{
   switch (o_)
   {
   case 1:
   {
      const Node& n0 = mesh.get_node(node_idx_[0]);
      const Node& n1 = mesh.get_node(node_idx_[1]);
      const Node& n2 = mesh.get_node(node_idx_[2]);
      if (mesh.get_mesh_dimension() == 2) {
         J << n1.x-n0.x,  n2.x-n0.x,
              n1.y-n0.y,  n2.y-n0.y;
      } else {
         J << n1.x-n0.x,  n2.x-n0.x,
              n1.y-n0.y,  n2.y-n0.y,
              n1.z-n0.z,  n2.z-n0.z;
      } 
      return true;
   }
   default:
      Logger::warning("Triangle::compute_Jacobian: higher order case not available.");
      return false;
      break;
   }
}


void Triangle::compute_shape(const Ref_Coord& coord, Eigen::Ref<VectorXd> shape) const
{
   double x = coord.x;
   double y = coord.y;
   switch(o_)
   {
      case 1:
         //  Barycentric coordinates λ
         //          λ0 = 1.0 - x - y;
         //          λ1 = x;
         //          λ2 = y;
         //
         if (shape.size() != 3)
         {
            throw std::invalid_argument("vector must be 3x1 for p-1 H(grad) Triangle.");
         }

         shape(0) = 1.0 - x - y;  // vertex 0
         shape(1) = x;            // vertex 1
         shape(2) = y;            // vertex 2
         break;
      default:
         throw std::invalid_argument("Nodal element not available for order: " + std::to_string(o_));
   }
}



void Triangle::compute_D_shape(const Ref_Coord& coord, Eigen::Ref<MatrixXd> d_shape) const
{
   switch (o_)
   {
   case 1:
      d_shape << -1, -1,
                  1,  0,
                  0,  1;
      break;
   
   default:
      Logger::warning("Triangle::compute_D_shape: higher order case not available.");
      break;
   }
}



std::vector<Ref_Coord> Triangle::edge_map(const Ref_Coord& edge_coord) const
{
   double x = edge_coord.x;
   // mapping: (1-x)*pi + x*pj
   // Edge 0: 0 -> 1
   // Edge 1: 0 -> 2
   // Edge 2: 1 -> 2
   return {{x  , 0, 0},
           {0  , x, 0},
           {1-x, x, 0},};
}

std::vector<Ref_Coord> Triangle::face_map(const Ref_Coord& face_coord) const
{
   double x = face_coord.x;
   double y = face_coord.y;
   // For a triangle element, the face is the element itself
   // so the mapping is the identity.
   return {{x, y, 0}};
}




void Triangle::tangent(Eigen::Ref<MatrixXd> t) const
{
   // 3 edges, 2D
   t << 1,  0,
        0,  1,
       -1,  1;
}


void Triangle::normal(Eigen::Ref<MatrixXd> n) const
{
   // 3 edges, 2D
   n << 1,  1,  // vertex 0 → normal of edge (1,2)
        1,  0,  // vertex 1 → normal of edge (0,2)
        0, -1;  // vertex 2 → normal of edge (0,1)
}