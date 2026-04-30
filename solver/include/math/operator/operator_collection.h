#pragma once

#include "math/operator/operator.h"
#include "math/operator/interpolator.h"
#include "math/operator/integrator.h"
#include "math/operator/integrator_a.h"
#include "math/operator/integrator_l.h"


namespace simu {


#include <mutex>

struct Check_Flag
{
    mutable std::array<std::once_flag, Integrator::SIZE> integrator_check;
    mutable std::array<std::once_flag, Interpolator::SIZE> interpolator_check;
};

}



