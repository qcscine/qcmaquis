add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_models dmrg_utils ${DMRG_LIBRARIES})

if(USE_AMBIENT)
  compile_with_ambient()
endif(USE_AMBIENT)

# *** Configure symmetries
configure_symm_file("maquis_dmrg_symm/maquis_dmrg_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp" CMAKE_SYMM_GROUP)
configure_file("maquis_dmrg_symm/maquis_dmrg_tpl_nu1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_nu1.cpp")
get_symm_files(APPEND DMRG_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp")


# *** Targets
add_library(maquis_dmrg SHARED ${DMRG_SYMM_SOURCES} maquis_dmrg.cpp maquis_cinterface.cpp mpssi_interface.cpp mpssi_cinterface.cpp starting_guess.cpp maquis_dmrg_detail.cpp complex_interface/maquis_cppinterface.cpp)
set_property(TARGET maquis_dmrg PROPERTY POSITION_INDEPENDENT_CODE TRUE)
if(APPLE)
  target_link_libraries(maquis_dmrg ${DMRG_APP_LIBRARIES} "-framework Accelerate")
# else()
#   target_link_libraries(maquis_dmrg ${DMRG_APP_LIBRARIES})
#   target_link_libraries(maquis_dmrg Eigen3::Eigen)
endif()

# *** Install
install(FILES maquis_dmrg_detail.h DESTINATION include/lib/maquis_dmrg COMPONENT headers)
install(FILES maquis_dmrg.h maquis_cinterface.h complex_interface/maquis_cppinterface.h DESTINATION include COMPONENT headers)
install(TARGETS maquis_dmrg EXPORT DMRGTargets COMPONENT libraries DESTINATION lib)
export(TARGETS maquis_dmrg APPEND FILE "${PROJECT_BINARY_DIR}/DMRGTargets.cmake")