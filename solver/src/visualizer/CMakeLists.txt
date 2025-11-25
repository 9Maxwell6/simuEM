

add_subdirectory(ui)
add_subdirectory(renderer)

set(Visualizer app_window.cpp)

target_sources(${PROJECT_GUI} PRIVATE ${Visualizer})
