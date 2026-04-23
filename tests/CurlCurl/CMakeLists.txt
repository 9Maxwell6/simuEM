

add_executable(test_curlcurl CurlCurl.cpp)
target_link_libraries(test_curlcurl PRIVATE simu_core)
add_test(NAME test_curlcurl COMMAND test_curlcurl)