
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})
if(USE_MTM)
  list(INSERT DMRG_APP_LIBRARIES 0 mtm)
endif(USE_MTM)

if(USE_AMBIENT)
  compile_with_ambient()
  list(INSERT DMRG_APP_LIBRARIES 0 ambient ${MPI_CXX_LIBRARIES})
endif(USE_AMBIENT)  

# *** Configure symmetries
configure_symm_file("dmrg_symm/dmrg_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/dmrg_symm/dmrg_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND DMRG_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/dmrg_symm/dmrg_{SYMM}.cpp")

configure_symm_file("measure_symm/measure_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/measure_symm/measure_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND MEASURE_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/measure_symm/measure_{SYMM}.cpp")


# *** Targets
add_library(dmrg_symmfactory STATIC symm_factory.cpp)

add_executable(dmrg ${DMRG_SYMM_SOURCES} dmrg.cpp)
target_link_libraries(dmrg dmrg_symmfactory ${DMRG_APP_LIBRARIES})

add_executable(dmrg_meas ${MEASURE_SYMM_SOURCES} measure.cpp)
target_link_libraries(dmrg_meas dmrg_symmfactory ${DMRG_APP_LIBRARIES})

add_executable(dmrg_cmeas ${MEASURE_SYMM_SOURCES} measure.cpp)
target_link_libraries(dmrg_cmeas dmrg_symmfactory ${DMRG_APP_LIBRARIES})
set_target_properties(dmrg_cmeas PROPERTIES COMPILE_FLAGS -DUSE_COMPLEX)

# *** Install
install(TARGETS dmrg dmrg_meas dmrg_cmeas RUNTIME DESTINATION bin COMPONENT applications)
