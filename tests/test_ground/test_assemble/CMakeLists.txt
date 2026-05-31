
add_executable(assemble assemble.cpp)
target_link_libraries(assemble PRIVATE simu_core)
add_test(NAME assemble COMMAND test_assemble)

