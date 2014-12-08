add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils ${DMRG_LIBRARIES})

if(USE_AMBIENT)
  compile_with_ambient()
endif(USE_AMBIENT)


macro(create_tools_symm_target TNAME SYMM CPPFILE LIBS)
  get_symm_group_name(SYMM_NAME ${SYMM})
  if(DMRG_HAS_${SYMM_NAME})
    add_executable(${TNAME} ${CPPFILE})
    target_link_libraries(${TNAME} ${LIBS})
    set_target_properties(${TNAME} PROPERTIES COMPILE_DEFINITIONS "USE_${SYMM}")
    install(TARGETS ${TNAME} RUNTIME DESTINATION bin COMPONENT applications)
  endif(DMRG_HAS_${SYMM_NAME})
endmacro(create_tools_symm_target)



create_tools_symm_target("mps_overlap_u1"   "U1"     "mps_overlap.cpp"  "${DMRG_APP_LIBRARIES}")
create_tools_symm_target("mps_overlap_2u1"  "TWOU1"  "mps_overlap.cpp"  "${DMRG_APP_LIBRARIES}")
create_tools_symm_target("mps_overlap_none" "NONE"   "mps_overlap.cpp"  "${DMRG_APP_LIBRARIES}")

# create_tools_symm_target("mps_reflect_u1"   "U1"     "mps_reflect.cpp"  "${DMRG_APP_LIBRARIES}")
create_tools_symm_target("mps_reflect_2u1"  "TWOU1"  "mps_reflect.cpp"  "${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mps_reflect_none" "NONE"   "mps_reflect.cpp"  "${DMRG_APP_LIBRARIES}")

create_tools_symm_target("mps_stats_u1"   "U1"     "mps_stats.cpp"  "${DMRG_APP_LIBRARIES}")
create_tools_symm_target("mps_stats_2u1"  "TWOU1"  "mps_stats.cpp"  "${DMRG_APP_LIBRARIES}")

create_tools_symm_target("mpo_stats_u1"      "U1"      "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
create_tools_symm_target("mpo_stats_2u1"     "TWOU1"   "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
create_tools_symm_target("mpo_stats_2u1pg"   "TWOU1PG" "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
create_tools_symm_target("mpo_stats_su2u1"   "SU2U1"   "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
create_tools_symm_target("mpo_stats_su2u1pg" "SU2U1PG" "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")


create_tools_symm_target("mps_upgrade_u1"   "U1"     "mps_upgrade.cpp"  "${DMRG_APP_LIBRARIES}")
create_tools_symm_target("mps_upgrade_2u1"  "TWOU1"  "mps_upgrade.cpp"  "${DMRG_APP_LIBRARIES}")
create_tools_symm_target("mps_upgrade_none" "NONE"   "mps_upgrade.cpp"  "${DMRG_APP_LIBRARIES}")


if(DMRG_HAS_SU2U1)
  set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
  #find_package(GSL REQUIRED)
  add_executable(su2corr su2corr.cpp)
  target_link_libraries(su2corr ${DMRG_APP_LIBRARIES} ${GSL_LIBRARIES})
  set_target_properties(su2corr PROPERTIES COMPILE_DEFINITIONS "USE_SU2U1")
endif(DMRG_HAS_SU2U1)

if(DMRG_HAS_SU2U1PG)
  set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
  #find_package(GSL REQUIRED)
  add_executable(su2corr_pg su2corr.cpp)
  target_link_libraries(su2corr_pg ${DMRG_APP_LIBRARIES} ${GSL_LIBRARIES})
  set_target_properties(su2corr_pg PROPERTIES COMPILE_DEFINITIONS "USE_SU2U1PG")
endif(DMRG_HAS_SU2U1PG)


if(USE_AMBIENT)
  ## enable only if needed, and with a revision < r3578
  # add_executable(mps_converter mps_converter.cpp)
  # target_link_libraries(mps_converter ${DMRG_APP_LIBRARIES})
  # set_target_properties(mps_converter PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")
endif(USE_AMBIENT)

if(DMRG_HAS_TwoU1PG)
  add_executable(nrdmexpect_2u1pg nrdmexpect.cpp)
  target_link_libraries(nrdmexpect_2u1pg ${DMRG_APP_LIBRARIES})
  set_target_properties(nrdmexpect_2u1pg PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1PG")
endif(DMRG_HAS_TwoU1PG)

