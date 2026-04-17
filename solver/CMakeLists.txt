# Add sub-projects
add_subdirectory(src/app)
add_subdirectory(src/world)
add_subdirectory(src/math)
add_subdirectory(src/physics)
add_subdirectory(src/visualizer)
add_subdirectory(src/utils)




# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE include)

if(ENABLE_GUI)
    target_include_directories(${PROJECT_GUI} PRIVATE include)
endif()