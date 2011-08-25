
#
# Find Package boost
#

add_definitions (-DBOOST_TEST_DYN_LINK)

#
# serial test cases
#

add_executable (vli_cpu_tests vli_cpu_tests.cpp)
target_link_libraries (vli_cpu_tests ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES} ${GMP_LIBRARIES} -lgmp)

add_executable (vli_gpu_tests vli_gpu_tests.cpp)
target_link_libraries (vli_gpu_tests kernels_gpu ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES} -lgmp -lgmpxx)

add_test (vli_cpu_tests vli_cpu_tests)
add_test (vli_gpu_tests vli_gpu_tests)

#
# monomial polynomial cpu - gpu test
#

add_executable ( vli_monomial_tests vli_monomial_tests.cpp)
target_link_libraries (vli_monomial_tests kernels_gpu ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES})

add_executable ( vli_monomial_polynomial_tests vli_monomial_polynomial_tests.cpp)
target_link_libraries (vli_monomial_polynomial_tests kernels_gpu ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES} ${GMP_LIBRARY_C} ${GMP_LIBRARY_CXX})

add_executable ( vli_polynomial_tests vli_polynomial_tests.cpp)
target_link_libraries (vli_polynomial_tests kernels_gpu ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES})


#add_executable ( vli_monomial_polynomial_tests vli_monomial_polynomial_tests.cpp)
#target_link_libraries (vli_monomial_polynomial_tests kernels_gpu ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES})

add_test(vli_monomial_tests vli_monomial_tests)
#add_test(vli_monomial_polynomial_tests vli_monomial_polynomial_tests)
add_test(vli_polynomial_tests vli_polynomial_tests)
#
# vector
#

#add_executable(vli_vector_cpu_gpu_tests vli_vector_cpu_gpu_tests.cpp)
#target_link_libraries( vli_vector_cpu_gpu_tests kernels_gpu ${Boost_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES})

#add_test(vli_vector_cpu_gpu_tests vli_vector_cpu_gpu_tests)
