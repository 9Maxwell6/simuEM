

add_executable(test_poisson Poisson.cpp)
target_link_libraries(test_poisson PRIVATE simu_core)
add_test(NAME test_poisson COMMAND test_poisson)