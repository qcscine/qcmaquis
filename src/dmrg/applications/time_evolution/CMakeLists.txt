
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})
if(USE_MTM)
  list(INSERT DMRG_APP_LIBRARIES 0 mtm)
endif(USE_MTM)

if(USE_AMBIENT)
  compile_with_ambient()
endif(USE_AMBIENT)  

# *** Configure symmetries
configure_symm_file("tevol_symm/dmrg_tevol_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/tevol_symm/dmrg_tevol_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND DMRG_TEVOL_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/tevol_symm/dmrg_tevol_{SYMM}.cpp")


# *** Targets
add_executable(dmrg_tevol ${DMRG_TEVOL_SYMM_SOURCES} dmrg_tevol.cpp)
target_link_libraries(dmrg_tevol ${DMRG_APP_LIBRARIES})


# *** Install
install(TARGETS dmrg_tevol RUNTIME DESTINATION bin COMPONENT applications)

