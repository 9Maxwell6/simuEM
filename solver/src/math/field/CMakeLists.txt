target_sources(${PROJECT_NAME} PRIVATE 
    field.cpp
    field_constant.cpp
    field_fespace.cpp
    field_function.cpp
    field_visualizer.cpp)

target_sources(${PROJECT_GUI} PRIVATE 
    field_visualizer.cpp)