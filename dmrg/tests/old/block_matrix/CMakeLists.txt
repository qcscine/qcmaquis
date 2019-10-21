
add_executable (multi_index.test multi_index.cpp)
add_test(multi_index multi_index.test)

add_test(NAME multi_index2
         # CONFIGURATIONS Debug
         COMMAND multi_index.test
         )
         
add_executable (check_multi_index.test check_multi_index.cpp)

