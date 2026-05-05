#pragma once
#include "math/ref_coord.h"
#include "math/fem/shape.h"


#include "utils/logger.h"

#include <vector>


namespace simu {

namespace Integration{

    const std::vector<Integration_Point>& get_rule(Basis_Shape b_shape, int order);


};

}


