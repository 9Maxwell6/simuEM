#include "entity/mesh/e_edge.h"

Type Edge::get_Type() const
{
    return Type::EDGE;
}

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

void Edge::set_NodeIdx(const size_t *idx)
{
   for (int i = 0; i < 2; i++)
   {
      node_idx_[i] = idx[i];
   }
}

void Edge::set_NodeIdx(size_t idx_1, size_t idx_2, size_t idx_3)
{
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
}

const size_t * Edge::get_NodeIdx() const
{
   return node_idx_;
}



