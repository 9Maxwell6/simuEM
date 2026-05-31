

add_executable(poisson Poisson.cpp)
target_link_libraries(poisson PRIVATE simu_core)
add_test(NAME poisson COMMAND test_poisson)