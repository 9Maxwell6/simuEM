#include "entity/mesh/mesh.h"
#include "entity/mesh/e_edge.h"

using namespace simu;

Geometry Edge::get_geometry() const { return Geometry::EDGE; }

Edge::Edge(const size_t *node_idx, size_t id, size_t property_id, int o) : Element(id, property_id, o)
{
   property_id_ = property_id;
   for (int i = 0; i < 2; i++)
   {
      node_idx_[i] = node_idx[i];
   }
}


Edge::Edge(std::vector<std::size_t> node_idx, size_t id, size_t property_id, int o) : Element(id, property_id, o)
{
   property_id_ = property_id;
   for (int i = 0; i < 2; i++)
   {
      node_idx_[i] = node_idx[i];
   }
}


Edge::Edge(size_t idx_1, size_t idx_2, size_t id, size_t property_id, int o) : Element(id, property_id, o)
{
    property_id_  = property_id;
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
}

void Edge::set_node_idx(const size_t *idx)
{
   for (int i = 0; i < 2; i++)
   {
      node_idx_[i] = idx[i];
   }
}

void Edge::set_node_idx(size_t idx_1, size_t idx_2, size_t idx_3)
{
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
}

const size_t * Edge::get_node_idx() const
{
   return node_idx_;
}


bool Edge::compute_Jacobian(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> J) const 
{
   Logger::error("Edge::compute_Jacobian - not implemented.");
   return false;
}

void Edge::compute_D_shape(const Mesh& mesh, const Ref_Coord& coord, Eigen::Ref<MatrixXd> d_shape) const
{
   Logger::error("Edge::compute_D_shape - not implemented.");
}

void Edge::compute_dof_transformation_H_curl(const Mesh& mesh, Eigen::Ref<MatrixXd> P) const
{
   switch (o_)
   {
   case 1:
   {
      // P is 1x1 diagonal matrix of +1/-1
      // Edge 0: 0 -> 1
      size_t idx_0 = node_idx_[0];
      size_t idx_1 = node_idx_[1];
      double s0 = (idx_0 < idx_1) ? 1. : -1.;
      P << s0;
      break;
   }
   default:
      Logger::warning("Edge::compute_dof_transformation_H_curl: higher order case not available.");
      break;
   }
}


void Edge::compute_dof_transformation_H_div(const Mesh& mesh, Eigen::Ref<MatrixXd> P) const
{
   // TODO...
   Logger::warning("Edge::compute_dof_transformation_H_div: not implemented.");
}

