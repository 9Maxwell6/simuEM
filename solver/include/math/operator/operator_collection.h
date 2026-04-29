#pragma once

#include "math/operator/operator.h"
#include "math/operator/interpolator.h"
#include "math/operator/integrator.h"
#include "math/operator/integrator_a.h"


namespace simu {


#include <mutex>

struct Check_Flag
{
    mutable std::array<std::once_flag, Integrator::SIZE> integrator_check;
    mutable std::array<std::once_flag, Interpolator::SIZE> interpolator_check;
};

}


// integrator_l.h require field_collection.h, but some of the field.h require element_data.h (they need physical coordinates inside the element).
// element_data.h require Check_Flag.
// Check_Flag require size of each operator types, which defined in integrator.h, interpolator.h, etc.
// so Check_Flag must be defined after [integrator.h, interpolator.h, etc.], and before integrator_l.h
// TODO: restructure all field.h to make this cleaner!
#include "math/operator/integrator_l.h"



