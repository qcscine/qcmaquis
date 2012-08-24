
add_definitions (-DBOOST_TEST_DYN_LINK -DVLI_POLYNOMIAL_PRINT_ZEROS)

foreach(test_vli
#        vli_monomial_tests
        polynomial_tests
       )
           add_executable(${test_vli} ${test_vli}.cpp)
           target_link_libraries(${test_vli} ${VLI_TEST_LIBRARIES} ${OpenMP_CXX_FLAGS})
           add_test(${test_vli} ${test_vli}) 
endforeach(test_vli)

file(COPY polynomial_tests_reference.txt DESTINATION .)
