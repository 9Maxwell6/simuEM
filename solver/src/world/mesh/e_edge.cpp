#include "world/mesh/mesh.h"
#include "world/mesh/e_edge.h"

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



void Edge::compute_shape(const Ref_Coord& coord, Eigen::Ref<VectorXd> shape) const
{
   Logger::error("Edge::compute_shape - not implemented.");
}



void Edge::compute_D_shape(const Ref_Coord& coord, Eigen::Ref<MatrixXd> d_shape) const
{
   Logger::error("Edge::compute_D_shape - not implemented.");
}


std::vector<Ref_Coord> Edge::edge_map(const std::vector<Integration_Point>& edge_coord, size_t edge_idx) const
{
   std::vector<Ref_Coord> coord(edge_coord.size());
   switch (edge_idx)
   {
   case 0:
      for(const Integration_Point& i_p : edge_coord) coord.push_back({i_p.coord.x, 0, 0});
      break;
   
   default:
      Logger::error("Edge::edge_map - edge_idx exceed limit: {0}.");
      return {};
   }
   return coord;
   
}

std::vector<Ref_Coord> Edge::face_map(const std::vector<Integration_Point>& face_coord, size_t face_idx) const
{
   Logger::error("Edge::face_map - edge element does not have face_map.");
   return {};
}


void Edge::tangent(Eigen::Ref<MatrixXd> t) const
{
   // 1 edges, 1D
   t << 1;
}


void Edge::normal(Eigen::Ref<MatrixXd> n) const
{
   Logger::error("Edge::normal - 1D element cannot have normal.");
}

