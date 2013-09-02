add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils ${DMRG_LIBRARIES})

if(USE_AMBIENT)
  compile_with_ambient()
  list(INSERT DMRG_APP_LIBRARIES 0 ambient ${MPI_CXX_LIBRARIES})
endif(USE_AMBIENT)

add_executable(mps_overlap mps_overlap.cpp)
target_link_libraries(mps_overlap ${DMRG_APP_LIBRARIES})

add_executable(mps_stats mps_stats.cpp)
target_link_libraries(mps_stats ${DMRG_APP_LIBRARIES})
#set_target_properties(mps_stats PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")

add_executable(mps_upgrade_none mps_upgrade.cpp)
target_link_libraries(mps_upgrade_none ${DMRG_APP_LIBRARIES})
set_target_properties(mps_upgrade_none PROPERTIES COMPILE_DEFINITIONS "USE_NONE")

add_executable(mps_upgrade_u1 mps_upgrade.cpp)
target_link_libraries(mps_upgrade_u1 ${DMRG_APP_LIBRARIES})
set_target_properties(mps_upgrade_u1 PROPERTIES COMPILE_DEFINITIONS "USE_U1")

add_executable(mps_upgrade_2u1 mps_upgrade.cpp)
target_link_libraries(mps_upgrade_2u1 ${DMRG_APP_LIBRARIES})
set_target_properties(mps_upgrade_2u1 PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")

if(USE_AMBIENT)
  ## enable only if needed, and with a revision < r3578
  # add_executable(mps_converter mps_converter.cpp)
  # target_link_libraries(mps_converter ${DMRG_APP_LIBRARIES})
  # set_target_properties(mps_converter PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")
endif(USE_AMBIENT)

