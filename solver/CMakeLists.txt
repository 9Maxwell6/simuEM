# Add sub-projects
add_subdirectory(src/app)
add_subdirectory(src/world)
add_subdirectory(src/math)
add_subdirectory(src/physics)
add_subdirectory(src/utils)








target_include_directories(simu_core PUBLIC include)