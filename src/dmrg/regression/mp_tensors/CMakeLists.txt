add_definitions(-DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
add_executable (generic_reshape.test ${DMRG_FRAMEWORK_DIR}/dmrg/utils/utils.cpp generic_reshape.cpp)
add_test(NAME generic_reshape
         # CONFIGURATIONS Debug    
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/generic_reshape.test
         )
