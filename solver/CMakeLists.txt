





# Add sub-projects
add_subdirectory(src/app)
add_subdirectory(src/entity)
add_subdirectory(src/math)
add_subdirectory(src/physics)
add_subdirectory(src/visualizer)




# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE include)
target_include_directories(${PROJECT_GUI} PRIVATE include)