#include "mesh.h"

namespace simu {

namespace filter
{

}

/*  use lambda function to create filter

// Simple wrapper functions that return lambdas
auto filter_by_dimension(int dim) {
    return [dim](const Element* e) {
        return e->dimension() == dim;
    };
}

auto filter_by_volume(double min_vol, double max_vol) {
    return [min_vol, max_vol](const Element* e) {
        return e->volume() >= min_vol && e->volume() <= max_vol;
    };
}

auto filter_on_boundary(const Mesh* mesh) {
    return [mesh](const Element* e) {
        return mesh->is_on_boundary(e);
    };
}

// Usage
mark_elements(filter_by_dimension(3));
mark_elements(filter_by_volume(0.1, 1.0));
mark_elements(filter_on_boundary(this));

*/

}