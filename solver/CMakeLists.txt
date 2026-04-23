# Add sub-projects
add_subdirectory(src/app)
add_subdirectory(src/world)
add_subdirectory(src/math)
add_subdirectory(src/physics)
add_subdirectory(src/visualizer)
add_subdirectory(src/utils)






if(ENABLE_GUI)
    target_include_directories(${PROJECT_GUI} PRIVATE include)
endif()

target_include_directories(simu_core PUBLIC include)