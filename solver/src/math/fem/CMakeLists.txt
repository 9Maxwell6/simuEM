target_sources(${PROJECT_NAME} PRIVATE 
    block_rack.cpp
    fem_space.cpp
    fem_system.cpp
    space_H1.cpp
    space_Hcurl.cpp
    operation.cpp
    integration.cpp
    integrator_a.cpp
    integrator_l.cpp
    assemble_data.cpp)


if(ENABLE_GUI)
    target_sources(${PROJECT_GUI} PRIVATE 
    fem_space.cpp
    space_H1.cpp
    space_Hcurl.cpp
    integration.cpp
)

endif()