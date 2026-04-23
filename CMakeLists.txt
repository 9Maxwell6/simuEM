cmake_minimum_required(VERSION 3.22)
project(simuEM VERSION 1.0.0 LANGUAGES CXX C)
set(PROJECT_GUI "${PROJECT_NAME}_gui")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prohibit in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are prohibited.")
endif ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

# Load local configurations
if(EXISTS "${CMAKE_SOURCE_DIR}/LocalConfig.cmake")
    message(STATUS "Loading local configuration from LocalConfig.cmake")
    include("${CMAKE_SOURCE_DIR}/LocalConfig.cmake")
endif()


enable_language(C)
enable_language(CXX)



# ===============================================================================
# check EIGEN_DIR
# ===============================================================================
if(NOT EXISTS "${EIGEN_DIR}/Eigen/Core")
    message(FATAL_ERROR "EIGEN_DIR invalid: cannot find Eigen/Core in ${EIGEN_DIR}")
endif()


# ===============================================================================
# import MFEM
# ===============================================================================
if(USE_MFEM)
    if(MFEM_DIR)
        find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}" "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
    else()
        find_package(mfem REQUIRED NAMES MFEM)
    endif()
    message(STATUS "Found mfem config in: ${MFEM_DIR} (version ${MFEM_VERSION})")
    include_directories(${MFEM_INCLUDE_DIRS})

endif()
# -------------------------------------------------------------------------------




# ===============================================================================
# import Gmsh
# ===============================================================================
if(USE_GMSH)

    find_library(GMSH_LIB NAMES gmsh PATHS ${GMSH_DIR}/lib)
    add_library(gmsh SHARED IMPORTED)
    set_target_properties(gmsh PROPERTIES
        IMPORTED_LOCATION "${GMSH_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${GMSH_DIR}/include"
    )
    message(STATUS "Found Gmsh-sdk")
    
endif()
# ===============================================================================





# ===============================================================================
# import PETSc
# 
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#
#
# ===============================================================================
if(USE_PETSC)

    if(PETSC_DIR AND PETSC_ARCH)
        find_package(PkgConfig REQUIRED)
        set(ENV{PKG_CONFIG_PATH} "${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig")
        pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)
        message(STATUS "PETSC_FOUND: ${PETSC_FOUND}")
        message(STATUS "PETSC_INCLUDE_DIRS: ${PETSC_INCLUDE_DIRS}")
        message(STATUS "PETSC_CFLAGS: ${PETSC_CFLAGS}")
    else()
        FIND_PACKAGE(PETSc REQUIRED)
    endif()
    
endif()
# ===============================================================================





file(GLOB_RECURSE CORE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/solver/*.cpp"
)

list(FILTER CORE_SOURCES EXCLUDE REGEX ".*/main\\.cpp$")
list(FILTER CORE_SOURCES EXCLUDE REGEX ".*/main_GUI\\.cpp$")

add_library(simu_core STATIC ${CORE_SOURCES})

# Collect source files
add_subdirectory(solver)






# Cross-platform linking
if(UNIX)

    find_package(OpenMP REQUIRED)
    target_link_libraries(simu_core PUBLIC OpenMP::OpenMP_CXX)

    set_target_properties(${PROJECT_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )

    

    
    
    if(ENABLE_GUI)
        #target_link_libraries(${PROJECT_GUI} PRIVATE glfw vulkan)

        set_target_properties(${PROJECT_GUI} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
        )
        target_link_libraries(${PROJECT_GUI} PRIVATE glfw vulkan dl pthread X11 Xxf86vm Xrandr Xi)
    endif()


    


    if(USE_MFEM)
        target_link_libraries(simu_core PUBLIC ${MFEM_LIBRARIES})
        set(LOAD_MFEM ON)
    endif()

    if(USE_GMSH)
        target_link_libraries(simu_core PUBLIC gmsh)
        set(LOAD_GMSH ON)
        if(ENABLE_GUI)
            target_link_libraries(${PROJECT_GUI} PRIVATE gmsh)
        endif()
    endif()

    if(USE_PETSC)
        set(PETSC_DATA_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/data/petsc_data" CACHE PATH "Directory for petsc output data")
        file(MAKE_DIRECTORY ${PETSC_DATA_OUTPUT_DIR})
        
        target_link_libraries(simu_core PUBLIC PkgConfig::PETSC)

        target_compile_definitions(simu_core PUBLIC PETSC_DATA_OUTPUT_DIR="${PETSC_DATA_OUTPUT_DIR}")
        target_compile_definitions(simu_core PUBLIC LOAD_PETSC)

        #set(LOAD_PETSC ON)

        # Check if PETSc was built with MPI
        file(STRINGS "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h" _petsc_mpi REGEX "#define PETSC_MPICC_SHOW")
        if(_petsc_mpi)
            if(DEFINED MPI_DIR)
                set(MPI_HOME ${MPI_DIR})
        else()
            # Extract MPI path from PETSC_MPICC_SHOW
            string(REGEX MATCH "-I([^ ]+)" _ "${_petsc_mpi}")
            string(REPLACE "/include" "" _extracted_mpi_dir "${CMAKE_MATCH_1}")
            set(MPI_HOME ${_extracted_mpi_dir})
            message(STATUS "MPI path extracted from petscconf.h: ${MPI_HOME}")
        endif()
            find_package(MPI REQUIRED)
            target_link_libraries(simu_core PUBLIC MPI::MPI_C MPI::MPI_CXX)
            message(STATUS "MPI found: ${MPI_HOME}")
        endif()
    endif()



    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/config/config.h.in
        ${CMAKE_CURRENT_BINARY_DIR}/config/config.h 
    )

    target_include_directories(simu_core PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/config")
    target_include_directories(simu_core PUBLIC ${EIGEN_DIR})

    if(ENABLE_GUI)
        target_include_directories(${PROJECT_GUI} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/config")
        target_include_directories(${PROJECT_GUI} PRIVATE ${EIGEN_DIR})
    endif()

    #target_compile_definitions(${PROJECT_NAME} PRIVATE 
    #    SCRIPT_PATH="${CMAKE_SOURCE_DIR}/scripts/"
    #)

    target_compile_definitions(simu_core PUBLIC RUNTIME_OUTPUT_PATH="${CMAKE_BINARY_DIR}/bin/")
    target_compile_definitions(${PROJECT_GUI} PRIVATE RUNTIME_OUTPUT_PATH="${CMAKE_BINARY_DIR}/bin/")


    target_compile_definitions(simu_core PUBLIC SCRIPT_PATH="${CMAKE_SOURCE_DIR}/scripts/")


    # for debug
    set(DEBUG_DATA_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/data/debug_data" CACHE PATH "Directory for debug output data")
    file(MAKE_DIRECTORY ${DEBUG_DATA_OUTPUT_DIR})

    target_compile_definitions(simu_core PUBLIC DEBUG_DATA_OUTPUT_DIR="${DEBUG_DATA_OUTPUT_DIR}")



elseif(WIN32)

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
if(BUILD_TEST)
    set(TEST_DATA_OUTPUT_DIR "${CMAKE_SOURCE_DIR}/data/test_data" CACHE PATH "Directory for test output data")
    file(MAKE_DIRECTORY ${TEST_DATA_OUTPUT_DIR})

    target_compile_definitions(simu_core PUBLIC TEST_DATA_OUTPUT_DIR="${TEST_DATA_OUTPUT_DIR}")
    add_subdirectory(tests)
endif()

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()



