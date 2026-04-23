
add_subdirectory(ui)
add_subdirectory(renderer)


target_sources(${PROJECT_GUI} PRIVATE 
    app_window.cpp
)


