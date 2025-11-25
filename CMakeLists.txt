cmake_minimum_required(VERSION 3.22)
project(simuEM VERSION 1.0.0 LANGUAGES CXX C)
set(PROJECT_GUI "${PROJECT_NAME}_gui")


# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build options
option(BUILD_TESTS "Build test suite" ON)
option(BUILD_EXAMPLES "Build example applications" ON)


# Collect source files
add_subdirectory(solver)



# Set target properties
set_target_properties(${PROJECT_NAME} ${PROJECT_GUI} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Cross-platform linking
if(UNIX)
    # On Linux/macOS, link libdl and pthread if needed
    #target_link_libraries(${PROJECT_GUI} PRIVATE glfw vulkan dl pthread X11 Xxf86vm Xrandr Xi)
    target_link_libraries(${PROJECT_GUI} PRIVATE glfw vulkan)


    # link gmsh
    #set(gmsh "${CMAKE_SOURCE_DIR}/vendor/gmsh/lib/libgmsh.so")

    #set_target_properties(gmsh PROPERTIES
    #    IMPORTED_LOCATION "${GMSH_LIB}"
    #    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/vendor/gmsh/include"
    #)

    #target_link_libraries(${PROJECT_NAME} PRIVATE gmsh)

    set(GMSH_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/vendor/gmsh/include")
    set(GMSH_LIB_DIR "${CMAKE_SOURCE_DIR}/vendor/gmsh/lib")
    set(GMSH_LIB_PATH "${GMSH_LIB_DIR}/libgmsh.so.4.15.0")
    add_library(gmsh SHARED IMPORTED)
    set_target_properties(gmsh PROPERTIES
        IMPORTED_LOCATION "${GMSH_LIB_PATH}"
        INTERFACE_INCLUDE_DIRECTORIES "${GMSH_INCLUDE_DIR}"
    )
    target_link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/vendor/gmsh/lib)   #/vendor/gmsh/lib/

    target_link_libraries(${PROJECT_NAME} PRIVATE gmsh)



elseif(WIN32)
    # On Windows, these are not needed
endif()



# Link dependencies (handled in Dependencies.cmake)
#link_dependencies(${PROJECT_NAME})
#link_dependencies(${PROJECT_GUI})





if(CMAKE_SYSTEM MATCHES Windows)
   message(STATUS "Target system is Windows")
endif()

if(CMAKE_HOST_SYSTEM MATCHES Linux)
   message(STATUS "Build host runs Linux")
endif()




# Optional: Add tests and examples
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()



