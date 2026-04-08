#include "entity/mesh/mesh.h"
#include "entity/mesh/e_tetrahedron.h"

using namespace simu;

Geometry Tetrahedron::get_geometry() const { return Geometry::TETRAHEDRON; }

Tetrahedron::Tetrahedron(const size_t *node_idx, size_t id, size_t property_id, int o) : Element(id, property_id, o)
{
   property_id_ = property_id;
   for (int i = 0; i < 4; i++)
   {
      node_idx_[i] = node_idx[i];
   }
}

Tetrahedron::Tetrahedron(std::vector<std::size_t> node_idx, size_t id, size_t property_id, int o) : Element(id, property_id, o)
{
   property_id_ = property_id;
   for (int i = 0; i < 4; i++)
   {
      node_idx_[i] = node_idx[i];
   }
}

Tetrahedron::Tetrahedron(size_t idx_1, size_t idx_2, size_t idx_3, size_t idx_4, size_t id, size_t property_id, int o) : Element(id, property_id, o)
{
    property_id_  = property_id;
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
    node_idx_[2] = idx_3;
    node_idx_[3] = idx_4;
}

void Tetrahedron::set_node_idx(const size_t *idx)
{
   for (int i = 0; i < 4; i++)
   {
      node_idx_[i] = idx[i];
   }
}

void Tetrahedron::set_node_idx(size_t idx_1, size_t idx_2, size_t idx_3, size_t idx_4)
{
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
    node_idx_[2] = idx_3;
    node_idx_[3] = idx_4;
}

const size_t * Tetrahedron::get_node_idx() const
{
   return node_idx_;
}


bool Tetrahedron::compute_Jacobian(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> J) const 
{
   switch (o_)
   {
   case 1:
   {
      const Node& n0 = mesh.get_node(node_idx_[0]);
      const Node& n1 = mesh.get_node(node_idx_[1]);
      const Node& n2 = mesh.get_node(node_idx_[2]);
      const Node& n3 = mesh.get_node(node_idx_[3]);

      J << n1.x-n0.x,  n2.x-n0.x,  n3.x-n0.x,
           n1.y-n0.y,  n2.y-n0.y,  n3.y-n0.y,
           n1.z-n0.z,  n2.z-n0.z,  n3.z-n0.z;
      
      return true;
   }
   default:
      Logger::warning("Tetrahedron::compute_Jacobian: higher order case not available.");
      return false;
      break;
   }
}



void Tetrahedron::compute_shape(const Ref_Coord& coord, Eigen::Ref<VectorXd> shape) const
{
   double x = coord.x;
   double y = coord.y;
   double z = coord.z;
   switch(o_)
   {
      case 1:
         //  Barycentric coordinates λ
         //          λ0 = 1.0 - x - y - z;
         //          λ1 = x;
         //          λ2 = y;
         //          λ3 = z;
         //
         if (shape.size() != 4) 
         {
            throw std::invalid_argument("vector must be 4x1 for p-1 H(grad) Tetrahedron.");
         }
         

         shape(0) = 1.0 - x - y - z;  // Vertex 0
         shape(1) = x;                // Vertex 1
         shape(2) = y;                // Vertex 2
         shape(3) = z;                // Vertex 3
         break;
      default:
         throw std::invalid_argument("Nodal element not available for order:  "+std::to_string(o_));
   }
}



void Tetrahedron::compute_D_shape(const Ref_Coord& coord, Eigen::Ref<MatrixXd> d_shape) const
{
   switch (o_)
   {
   case 1:
      d_shape << -1, -1, -1,
                  1,  0,  0,
                  0,  1,  0,
                  0,  0,  1;
      break;
   
   default:
      Logger::warning("Tetrahedron::compute_D_shape: higher order case not available.");
      break;
   }
}

void Tetrahedron::compute_dof_transformation_H_curl(const Mesh& mesh, Eigen::Ref<MatrixXd> P) const
{
   switch (o_)
   {
   case 1:
   {
      // P is 6x6 diagonal matrix of +1/-1
      // Edge 0: 0 -> 1
      // Edge 1: 0 -> 2
      // Edge 2: 0 -> 3
      // Edge 3: 1 -> 2
      // Edge 4: 1 -> 3
      // Edge 5: 2 -> 3
      size_t idx_0 = node_idx_[0];
      size_t idx_1 = node_idx_[1];
      size_t idx_2 = node_idx_[2];
      size_t idx_3 = node_idx_[3];

      double s0 = (idx_0 < idx_1) ? 1. : -1.;
      double s1 = (idx_0 < idx_2) ? 1. : -1.;
      double s2 = (idx_0 < idx_3) ? 1. : -1.;
      double s3 = (idx_1 < idx_2) ? 1. : -1.;
      double s4 = (idx_1 < idx_3) ? 1. : -1.;
      double s5 = (idx_2 < idx_3) ? 1. : -1.;

      P << s0,  0,  0,  0,  0,  0,
            0, s1,  0,  0,  0,  0,
            0,  0, s2,  0,  0,  0,
            0,  0,  0, s3,  0,  0,
            0,  0,  0,  0, s4,  0,
            0,  0,  0,  0,  0, s5;
      break;
   }
   default:
      Logger::warning("Tetrahedron::compute_dof_transformation_H_curl: higher order case not available.");
      break;
   }
}


void Tetrahedron::compute_dof_transformation_H_div(const Mesh& mesh, Eigen::Ref<MatrixXd> P) const
{
   // TODO...
   Logger::warning("Tetrahedron::compute_dof_transformation_H_div: not implemented.");
}

