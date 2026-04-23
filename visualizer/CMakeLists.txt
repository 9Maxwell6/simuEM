add_subdirectory(src/app)
add_subdirectory(src/renderer)
add_subdirectory(src/ui)



target_include_directories(${PROJECT_GUI} PUBLIC include)