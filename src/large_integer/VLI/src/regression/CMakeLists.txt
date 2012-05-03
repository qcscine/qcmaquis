#
# Find Package boost
#
add_definitions (-DBOOST_TEST_DYN_LINK)

add_executable (vli_cpu_tests vli_cpu_tests.cpp)
add_executable ( vli_monomial_tests vli_monomial_tests.cpp)
add_executable ( vli_polynomial_tests vli_polynomial_tests.cpp)
add_executable (vli_vector_cpu_gpu_tests vli_vector_cpu_gpu_tests.cpp)

foreach(test_vli
        vli_cpu_tests
        vli_monomial_tests
        vli_polynomial_tests
        vli_vector_cpu_gpu_tests
       )
           target_link_libraries(${test_vli} ${Boost_LIBRARIES} ${GMP_LIBRARY_CXX} ${PROJECT_BINARY_DIR}/src/vli/libvli_cpu.a ${OpenMP_CXX_FLAGS})
           add_test(${test_vli} ${test_vli}) 
endforeach(test_vli)

if(VLI_COMPILE_GPU)
    add_executable (vli_gpu_tests vli_gpu_tests.cpp)
    target_link_libraries(vli_gpu_tests ${Boost_LIBRARIES} ${GMP_LIBRARY_CXX} ${PROJECT_BINARY_DIR}/src/vli/libvli_gpu.a ${CUDA_LIBRARIES} ${PROJECT_BINARY_DIR}/src/vli/libvli_cpu.a ${OpenMP_CXX_FLAGS})
    add_test(vli_gpu_tests vli_gpu_tests) 
endif(VLI_COMPILE_GPU)
