
target_sources(simu_core PRIVATE 
    auxiliary_solver.cpp
    T_Omega.cpp)



add_executable(t_omega_c T_Omega_complex.cpp)
target_sources(t_omega_c PRIVATE 
    auxiliary_solver_complex.cpp)
target_link_libraries(t_omega_c PRIVATE simu_core)
#set_target_properties(t_omega_c PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set_target_properties(t_omega_c PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


