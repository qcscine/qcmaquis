add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils ${DMRG_LIBRARIES})

if(USE_AMBIENT)
  compile_with_ambient()
endif(USE_AMBIENT)

add_executable(mps_overlap mps_overlap.cpp)
target_link_libraries(mps_overlap ${DMRG_APP_LIBRARIES})

add_executable(mps_stats mps_stats.cpp)
target_link_libraries(mps_stats ${DMRG_APP_LIBRARIES})

if(DMRG_HAS_TwoU1)
  add_executable(mps_stats_2u1 mps_stats.cpp)
  target_link_libraries(mps_stats_2u1 ${DMRG_APP_LIBRARIES})
  set_target_properties(mps_stats_2u1 PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")
endif(DMRG_HAS_TwoU1)

#add_executable(mpo_stats_none mpo_stats.cpp)
#target_link_libraries(mpo_stats_none dmrg_models dmrg_utils ${DMRG_LIBRARIES})
#set_target_properties(mpo_stats_none PROPERTIES COMPILE_DEFINITIONS "USE_NONE")

#add_executable(mpo_stats_u1 mpo_stats.cpp)
#target_link_libraries(mpo_stats_u1 dmrg_models dmrg_utils ${DMRG_LIBRARIES})
#set_target_properties(mpo_stats_u1 PROPERTIES COMPILE_DEFINITIONS "USE_U1")

if(DMRG_HAS_TwoU1)
  add_executable(mpo_stats_2u1 mpo_stats.cpp)
  target_link_libraries(mpo_stats_2u1 dmrg_models dmrg_utils ${DMRG_LIBRARIES})
  set_target_properties(mpo_stats_2u1 PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")
endif(DMRG_HAS_TwoU1)

if(DMRG_HAS_TrivialGroup)
  add_executable(mps_upgrade_none mps_upgrade.cpp)
  target_link_libraries(mps_upgrade_none ${DMRG_APP_LIBRARIES})
  set_target_properties(mps_upgrade_none PROPERTIES COMPILE_DEFINITIONS "USE_NONE")
endif(DMRG_HAS_TrivialGroup)

if(DMRG_HAS_U1)
  add_executable(mps_upgrade_u1 mps_upgrade.cpp)
  target_link_libraries(mps_upgrade_u1 ${DMRG_APP_LIBRARIES})
  set_target_properties(mps_upgrade_u1 PROPERTIES COMPILE_DEFINITIONS "USE_U1")
endif(DMRG_HAS_U1)

if(DMRG_HAS_TwoU1)
  add_executable(mps_upgrade_2u1 mps_upgrade.cpp)
  target_link_libraries(mps_upgrade_2u1 ${DMRG_APP_LIBRARIES})
  set_target_properties(mps_upgrade_2u1 PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")
endif(DMRG_HAS_TwoU1)

if(USE_AMBIENT)
  ## enable only if needed, and with a revision < r3578
  # add_executable(mps_converter mps_converter.cpp)
  # target_link_libraries(mps_converter ${DMRG_APP_LIBRARIES})
  # set_target_properties(mps_converter PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")
endif(USE_AMBIENT)

