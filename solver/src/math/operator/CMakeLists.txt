add_subdirectory(quadrature)

target_sources(simu_core PRIVATE
    integration.cpp
    operator.cpp
    integrator_a.cpp
    integrator_l.cpp
    interpolator.cpp)

