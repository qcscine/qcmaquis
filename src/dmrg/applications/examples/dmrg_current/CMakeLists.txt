
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})

# *** Configure symmetries
configure_symm_file("symm/runsim_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/symm/runsim_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND RUNSIM_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/symm/runsim_{SYMM}.cpp")


# *** Targets
add_executable(dmrg_current ${RUNSIM_SYMM_SOURCES} dmrg_current.cpp)
target_link_libraries(dmrg_current ${DMRG_APP_LIBRARIES})

