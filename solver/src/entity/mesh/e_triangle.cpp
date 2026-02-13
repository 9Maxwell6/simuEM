#include "entity/mesh/e_triangle.h"

using namespace simu;

Type Triangle::get_Type() const
{
    return Type::TRIANGLE;
}

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

void Triangle::set_nodeIdx(const size_t *idx)
{
   for (int i = 0; i < 3; i++)
   {
      node_idx_[i] = idx[i];
   }
}

void Triangle::set_nodeIdx(size_t idx_1, size_t idx_2, size_t idx_3)
{
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
    node_idx_[2] = idx_3;
}

const size_t * Triangle::get_nodeIdx() const
{
   return node_idx_;
}



