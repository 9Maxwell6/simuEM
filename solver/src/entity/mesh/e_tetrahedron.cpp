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


bool Tetrahedron::compute_Jacobian(const Mesh& mesh, const Integration_Point& i_p, Eigen::Ref<MatrixXd> J) const 
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

void Tetrahedron::compute_D_shape(const Mesh& mesh, const Integration_Point& i_p, Eigen::Ref<MatrixXd> d_shape) const
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



