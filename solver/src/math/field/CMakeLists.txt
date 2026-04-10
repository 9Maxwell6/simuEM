target_sources(${PROJECT_NAME} PRIVATE 
    field.cpp
    field_constant.cpp
    field_fespace.cpp
    field_function.cpp
    field_visualizer.cpp)

if(ENABLE_GUI)
    target_sources(${PROJECT_GUI} PRIVATE 
        field_visualizer.cpp)
endif()