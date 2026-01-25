# Create main executable
add_executable(${PROJECT_NAME} main.cpp)

if(ENABLE_GUI)
    add_executable(${PROJECT_GUI} main_GUI.cpp)
endif()