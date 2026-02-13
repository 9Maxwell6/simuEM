#include "entity/mesh/e_tetrahedron.h"

using namespace simu;

Type Tetrahedron::get_Type() const
{
    return Type::TETRAHEDRON;
}

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

void Tetrahedron::set_nodeIdx(const size_t *idx)
{
   for (int i = 0; i < 4; i++)
   {
      node_idx_[i] = idx[i];
   }
}

void Tetrahedron::set_nodeIdx(size_t idx_1, size_t idx_2, size_t idx_3, size_t idx_4)
{
    node_idx_[0] = idx_1;
    node_idx_[1] = idx_2;
    node_idx_[2] = idx_3;
    node_idx_[3] = idx_4;
}

const size_t * Tetrahedron::get_nodeIdx() const
{
   return node_idx_;
}


