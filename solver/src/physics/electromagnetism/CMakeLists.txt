add_subdirectory(formulation)

target_sources(${PROJECT_NAME} PRIVATE 
    eddy_current.cpp
    mfem_eddy_current.cpp)