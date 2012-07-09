
add_definitions(-DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

# *** Targets
add_executable(diag diag.cpp)
target_link_libraries(diag ${DMRG_LIBRARIES})


# *** Install
install(TARGETS diag RUNTIME DESTINATION bin COMPONENT applications)
