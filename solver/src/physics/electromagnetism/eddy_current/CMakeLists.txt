
target_sources(simu_core PRIVATE 
    auxiliary_solver.cpp
    T_Omega.cpp
    T_Omega_complex.cpp)


#add_executable(T_Omega T_Omega.cpp)
#target_link_libraries(T_Omega PRIVATE simu_core)

#add_executable(T_Omega_complex T_Omega_complex.cpp)
#target_link_libraries(T_Omega_complex PRIVATE simu_core)

#target_sources(T_Omega PRIVATE 
#    auxiliary_solver.cpp
#)

#target_sources(T_Omega_complex PRIVATE 
#    auxiliary_solver.cpp
#)