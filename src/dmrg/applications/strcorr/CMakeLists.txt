
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})
if(USE_MTM)
  list(INSERT DMRG_APP_LIBRARIES 0 mtm)
endif(USE_MTM)


# *** Targets
add_executable(dmrg_strcorr dmrg_strcorr.cpp)
target_link_libraries(dmrg_strcorr ${DMRG_APP_LIBRARIES})

add_executable(dmrg_cstrcorr dmrg_strcorr.cpp)
target_link_libraries(dmrg_cstrcorr ${DMRG_APP_LIBRARIES})
set_target_properties(dmrg_cstrcorr PROPERTIES COMPILE_FLAGS -DUSE_COMPLEX)


# *** Install
install(TARGETS dmrg_strcorr dmrg_cstrcorr RUNTIME DESTINATION bin COMPONENT applications)
