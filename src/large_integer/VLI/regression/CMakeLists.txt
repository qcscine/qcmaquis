#
# Find Package boost
#
add_definitions (-DBOOST_TEST_DYN_LINK)

foreach(test_vli
        vli_monomial_tests
        vli_polynomial_tests
        vli_vector_cpu_gpu_tests
       )
           add_executable(${test_vli} ${test_vli}.cpp)
           target_link_libraries(${test_vli} ${Boost_LIBRARIES} ${GMP_LIBRARY_CXX} ${PROJECT_BINARY_DIR}/src/vli/libvli_cpu.a ${OpenMP_CXX_FLAGS})
           add_test(${test_vli} ${test_vli}) 
endforeach(test_vli)
