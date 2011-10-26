add_definitions(-DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
add_executable (generic_reshape.test generic_reshape.cpp)
add_test(generic_reshape ${CMAKE_CURRENT_BINARY_DIR}/generic_reshape.test)
