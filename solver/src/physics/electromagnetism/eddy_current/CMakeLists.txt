
target_sources(simu_core PRIVATE 
    auxiliary_solver.cpp
    T_Omega.cpp)



add_executable(t_omega_c T_Omega_complex.cpp)
target_sources(t_omega_c PRIVATE 
    auxiliary_solver_complex.cpp)
target_link_libraries(t_omega_c PRIVATE simu_core)
set_target_properties(t_omega_c PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(t_omega_c_1 T_Omega_c.cpp)
target_sources(t_omega_c_1 PRIVATE 
    auxiliary_solver_c.cpp)
target_link_libraries(t_omega_c_1 PRIVATE simu_core)
set_target_properties(t_omega_c_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(t_omega_c_2 T_Omega_c2.cpp)
target_sources(t_omega_c_2 PRIVATE 
    auxiliary_solver_c2.cpp)
target_link_libraries(t_omega_c_2 PRIVATE simu_core)
set_target_properties(t_omega_c_2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


