
#
# Find Package boost
#

add_definitions (-DBOOST_TEST_DYN_LINK)

#
# serial test cases
#

add_executable (vli_cpu_tests vli_cpu_tests.cpp)
target_link_libraries (vli_cpu_tests ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES})

add_executable (vli_gpu_tests vli_gpu_tests.cpp)
target_link_libraries (vli_gpu_tests kernels_gpu ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES})

add_test (vli_cpu_tests vli_cpu_tests)
add_test (vli_gpu_tests vli_gpu_tests)


#
# vector test cases
#


add_executable ( vli_cpu_vector_tests vli_cpu_vector_tests.cpp)
target_link_libraries (vli_cpu_vector_tests ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES})

add_executable ( vli_gpu_vector_tests vli_gpu_vector_tests.cpp)
target_link_libraries (vli_gpu_vector_tests kernels_gpu ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES})

add_test (vli_cpu_vector_tests vli_cpu_vector_tests)
add_test (vli_gpu_vector_tests vli_gpu_vector_tests)


#
# monome cpu - gpu test
#


add_executable ( vli_monome_tests vli_monome_tests.cpp)
target_link_libraries (vli_monome_tests kernels_gpu ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES})

add_test (vli_monome_tests vli_monome_tests)


