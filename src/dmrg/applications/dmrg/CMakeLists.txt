
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})
if(USE_MTM)
  list(INSERT DMRG_APP_LIBRARIES 0 mtm)
endif(USE_MTM)


# *** Configure symmetries
configure_symm_file("dmrg_symm/dmrg_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/dmrg_symm/dmrg_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND DMRG_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/dmrg_symm/dmrg_{SYMM}.cpp")

configure_symm_file("measure_symm/measure_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/measure_symm/measure_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND MEASURE_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/measure_symm/measure_{SYMM}.cpp")


# *** Targets
add_executable(dmrg ${DMRG_SYMM_SOURCES} dmrg.cpp)
target_link_libraries(dmrg ${DMRG_APP_LIBRARIES})

add_executable(dmrg_meas ${MEASURE_SYMM_SOURCES} measure.cpp)
target_link_libraries(dmrg_meas ${DMRG_APP_LIBRARIES})

if(BUILD_OLD_TARGETS)
  add_executable(measure ${MEASURE_SYMM_SOURCES} measure.cpp)
  target_link_libraries(measure ${DMRG_APP_LIBRARIES})
endif(BUILD_OLD_TARGETS)


# *** Install
install(TARGETS dmrg dmrg_meas RUNTIME DESTINATION bin COMPONENT applications)

if(BUILD_OLD_TARGETS)
  install(TARGETS measure RUNTIME DESTINATION bin COMPONENT applications)
endif(BUILD_OLD_TARGETS)

