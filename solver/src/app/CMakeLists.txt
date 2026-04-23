# Create main executable
add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE simu_core)

if(ENABLE_GUI)
    add_executable(${PROJECT_GUI} main_GUI.cpp)
endif()