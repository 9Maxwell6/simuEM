#include "entity/tetrahedron.h"

Geometry Tetrahedron::get_Geometry() const
{
    return Geometry::TETRAHEDRON;
}

Tetrahedron::Tetrahedron(const size_t *idx, size_t property_id)
{
   property_id_ = property_id;
   for (int i = 0; i < 4; i++)
   {
      node_idx_[i] = idx[i];
   }
}

Tetrahedron::Tetrahedron(size_t idx_1, size_t idx_2, size_t idx_3, size_t idx_4, size_t property_id)
{
    property_id_  = property_id;
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
    node_idx_[2] = idx_3;
    node_idx_[3] = idx_4;
}

void Tetrahedron::set_NodeIdx(const size_t *idx)
{
   for (int i = 0; i < 4; i++)
   {
      node_idx_[i] = idx[i];
   }
}

void Tetrahedron::set_NodeIdx(size_t idx_1, size_t idx_2, size_t idx_3, size_t idx_4)
{
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
    node_idx_[2] = idx_3;
    node_idx_[3] = idx_4;
}

const size_t * Tetrahedron::get_NodeIdx() const
{
   return node_idx_;
}



