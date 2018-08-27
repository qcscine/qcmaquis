
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_models dmrg_utils ${DMRG_LIBRARIES})

if(USE_AMBIENT)
  compile_with_ambient()
endif(USE_AMBIENT)

# *** Configure symmetries
configure_symm_file("maquis_dmrg_symm/maquis_dmrg_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND DMRG_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp")


# *** Targets
add_library(maquis_dmrg SHARED ${DMRG_SYMM_SOURCES})
set_property(TARGET maquis_dmrg PROPERTY POSITION_INDEPENDENT_CODE TRUE)
target_link_libraries(maquis_dmrg ${DMRG_APP_LIBRARIES})

# *** Install
install(TARGETS maquis_dmrg EXPORT DMRGTargets COMPONENT libraries DESTINATION lib)
export(TARGETS maquis_dmrg APPEND FILE "${PROJECT_BINARY_DIR}/DMRGTargets.cmake")
