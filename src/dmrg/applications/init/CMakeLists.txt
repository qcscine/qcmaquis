
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_models dmrg_utils ${DMRG_LIBRARIES})

if(USE_AMBIENT)
  compile_with_ambient()
endif(USE_AMBIENT)  

# *** Configure symmetries
configure_symm_file("symm/dmrg_init_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/symm/dmrg_init_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND DMRGINIT_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/symm/dmrg_init_{SYMM}.cpp")


# *** Targets
add_executable(dmrg_init ${DMRGINIT_SYMM_SOURCES} dmrg_init.cpp)
target_link_libraries(dmrg_init ${DMRG_APP_LIBRARIES})

# *** Install
install(TARGETS dmrg_init RUNTIME DESTINATION bin COMPONENT applications)
