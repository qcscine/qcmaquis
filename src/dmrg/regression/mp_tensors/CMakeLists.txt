add_definitions(-DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
add_executable (generic_reshape.test ${CMAKE_SOURCE_DIR}/utils/utils.cpp generic_reshape.cpp)
add_test(generic_reshape ${CMAKE_CURRENT_BINARY_DIR}/generic_reshape.test)
