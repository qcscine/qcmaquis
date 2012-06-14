#
# Find Package boost
#
if(VLI_COMPILE_GPU)
    include_directories(${CUDA_INCLUDE_DIRS})
    add_definitions(-DVLI_USE_GPU -DVLI_SPLIT_PARAM=${VLI_SPLIT_PARAM} )
endif(VLI_COMPILE_GPU)

add_definitions (-DBOOST_TEST_DYN_LINK)

foreach(test_vli
        vli_monomial_tests
        polynomial_tests
#        vli_polynomial_tests
       vli_vector_cpu_gpu_tests
       )
           add_executable(${test_vli} ${test_vli}.cpp)
           target_link_libraries(${test_vli} ${Boost_LIBRARIES} ${GMP_LIBRARY_CXX} ${PROJECT_BINARY_DIR}/src/vli/libvli_cpu.a ${OpenMP_CXX_FLAGS})

           if(VLI_COMPILE_GPU)
               target_link_libraries(${test_vli} ${Boost_LIBRARIES} ${GMP_LIBRARY_CXX} ${PROJECT_BINARY_DIR}/src/vli/libvli_cpu.a ${PROJECT_BINARY_DIR}/src/vli/libvli_gpu.a ${CUDA_CUDART_LIBRARY} ${CUDA_LIBRARIES}  ${OpenMP_CXX_FLAGS})
           endif(VLI_COMPILE_GPU)

           add_test(${test_vli} ${test_vli}) 
endforeach(test_vli)

file(COPY polynomial_tests_reference.txt DESTINATION .)
