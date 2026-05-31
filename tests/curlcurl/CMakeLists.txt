

add_executable(curlcurl CurlCurl.cpp)
target_link_libraries(curlcurl PRIVATE simu_core)
add_test(NAME curlcurl COMMAND test_curlcurl)