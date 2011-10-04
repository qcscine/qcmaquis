
#
# Find Package boost
#

add_definitions (-DBOOST_TEST_DYN_LINK)

#
# serial test cases
#

add_executable (vli_cpu_tests vli_cpu_tests.cpp)
target_link_libraries (vli_cpu_tests ${Boost_LIBRARIES} ${GMP_LIBRARY_CXX})

add_test (vli_cpu_tests vli_cpu_tests)

#
# monomial polynomial cpu - gpu test
#

add_executable ( vli_monomial_tests vli_monomial_tests.cpp)
target_link_libraries (vli_monomial_tests vli_gpu ${Boost_LIBRARIES} ${CUDA_LIBRARIES})

add_executable ( vli_polynomial_tests vli_polynomial_tests.cpp)
target_link_libraries (vli_polynomial_tests vli_gpu ${Boost_LIBRARIES} ${CUDA_LIBRARIES})

add_test(vli_monomial_tests vli_monomial_tests)
add_test(vli_polynomial_tests vli_polynomial_tests)

#
# vector
#
add_executable(vli_vector_cpu_gpu_tests vli_vector_cpu_gpu_tests.cpp)
target_link_libraries( vli_vector_cpu_gpu_tests vli_gpu ${Boost_LIBRARIES} ${CUDA_LIBRARIES} ${OPENMP_LIBRARIES})

add_test(vli_vector_cpu_gpu_tests vli_vector_cpu_gpu_tests)
