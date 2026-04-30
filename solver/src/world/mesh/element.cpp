#include "world/mesh/mesh.h"

#include "world/mesh/element.h"




using namespace simu;


void Element::physical_point(const Mesh& mesh, const Ref_Coord& ref_coord, Eigen::Ref<VectorXd> phy_coord) const
{
    phy_coord.setZero();
    int node_num = get_node_num();
    const size_t * node_idx = get_node_idx();

    int geo_node_num = get_geometry_node_num();
    // TODO: higher geometry order not support yet, mesh object should provide extra nodes for higher order elements.
    // possible idea: using Bounding volume hierarchy (BVH) algorithm to quickly locate current element in mesh. then get all the higher order nodes.

    VectorXd shape(geo_node_num);

    compute_shape(ref_coord, shape);

    for(int i=0; i<node_num; ++i)
    {
        const Node& n = mesh.get_node(node_idx[i]);
        phy_coord += shape[i] * Vector<3>{n.x, n.y, n.z};  // for phy_dim < 8,  n.y and n.z can be zero
    }

    for(int j=node_num; j<geo_node_num; ++j ){
        Logger::error("Element_Data::physical_point - higher geometry order not support yet");
    }
}
