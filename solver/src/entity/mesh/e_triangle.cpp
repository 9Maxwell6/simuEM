#include "entity/mesh/e_triangle.h"

Type Triangle::get_Type() const
{
    return Type::TETRAHEDRON;
}

Triangle::Triangle(const size_t *idx, size_t property_id, int o) : Element(property_id, o)
{
   property_id_ = property_id;
   for (int i = 0; i < 3; i++)
   {
      node_idx_[i] = idx[i];
   }
}

Triangle::Triangle(size_t idx_1, size_t idx_2, size_t idx_3, size_t property_id, int o) : Element(property_id, o)
{
    property_id_  = property_id;
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
    node_idx_[2] = idx_3;
}

void Triangle::set_NodeIdx(const size_t *idx)
{
   for (int i = 0; i < 3; i++)
   {
      node_idx_[i] = idx[i];
   }
}

void Triangle::set_NodeIdx(size_t idx_1, size_t idx_2, size_t idx_3)
{
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
    node_idx_[2] = idx_3;
}

const size_t * Triangle::get_NodeIdx() const
{
   return node_idx_;
}



