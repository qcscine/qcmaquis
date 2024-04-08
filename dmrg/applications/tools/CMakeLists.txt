add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES maquis_dmrg dmrg_models dmrg_utils ${DMRG_LIBRARIES})

# add_library(srcas_utils srcas_utilities.cpp)
# set_property(TARGET srcas_utils PROPERTY POSITION_INDEPENDENT_CODE TRUE)

macro(create_tools_symm_target TNAME SYMM CPPFILE LIBS)
  get_symm_group_name(SYMM_NAME ${SYMM})
  if(DMRG_HAS_${SYMM_NAME})
    add_executable(${TNAME} ${CPPFILE})
    target_link_libraries(${TNAME} ${LIBS})
    set_target_properties(${TNAME} PROPERTIES COMPILE_DEFINITIONS "USE_${SYMM}")
    install(TARGETS ${TNAME} RUNTIME DESTINATION bin COMPONENT applications)
  endif(DMRG_HAS_${SYMM_NAME})
endmacro(create_tools_symm_target)

macro(create_tools_vib_target TNAME CPPFILE LIBS)
  if(BUILD_VIBRATIONAL)
    add_executable(${TNAME} ${CPPFILE})
    target_link_libraries(${TNAME} ${LIBS})
    install(TARGETS ${TNAME} RUNTIME DESTINATION bin COMPONENT applications)
  endif()
endmacro()

macro(create_tools_target TNAME CPPFILE LIBS)
    add_executable(${TNAME} ${CPPFILE})
    target_link_libraries(${TNAME} ${LIBS})
    install(TARGETS ${TNAME} RUNTIME DESTINATION bin COMPONENT applications)
endmacro()

if (BUILD_MPS_OVERLAP)
  # Real case
  create_tools_symm_target("mps_overlap_u1"       "U1"       "mps_overlap.cpp"  "${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_overlap_2u1"      "TWOU1"    "mps_overlap.cpp"  "${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_overlap_2u1pg"    "TWOU1PG"  "mps_overlap.cpp"  "${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_overlap_su2u1"    "SU2U1"    "mps_overlap.cpp"  "${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_overlap_su2u1pg"  "SU2U1PG"  "mps_overlap.cpp"  "${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_overlap_none"     "NONE"     "mps_overlap.cpp"  "${DMRG_APP_LIBRARIES}")
  # Complex case 
  create_tools_symm_target("mps_overlap_complex_u1"       "U1"       "mps_overlap_complex.cpp"  "${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_overlap_complex_2u1"      "TWOU1"    "mps_overlap_complex.cpp"  "${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_overlap_complex_2u1pg"    "TWOU1PG"  "mps_overlap_complex.cpp"  "${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_overlap_complex_su2u1"    "SU2U1"    "mps_overlap_complex.cpp"  "${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_overlap_complex_su2u1pg"  "SU2U1PG"  "mps_overlap_complex.cpp"  "${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_overlap_complex_none"     "NONE"     "mps_overlap_complex.cpp"  "${DMRG_APP_LIBRARIES}")
endif (BUILD_MPS_OVERLAP)

# create_tools_symm_target("mps_reflect_u1"   "U1"     "mps_reflect.cpp"  "${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mps_reflect_2u1"  "TWOU1"  "mps_reflect.cpp"  "${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mps_reflect_none" "None"   "mps_reflect.cpp"  "${DMRG_APP_LIBRARIES}")

# create_tools_symm_target("mps_join_u1"   "U1"     "mps_join.cpp"  "${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mps_join_2u1"  "TWOU1"  "mps_join.cpp"  "${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mps_join_none" "None"   "mps_join.cpp"  "${DMRG_APP_LIBRARIES}")

# create_tools_symm_target("mps_stats_u1"      "U1"      "mps_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mps_stats_2u1"     "TWOU1"   "mps_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mps_stats_2u1pg"   "TWOU1PG" "mps_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mps_stats_u1df"    "U1DG"    "mps_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mps_stats_su2u1"   "SU2U1"   "mps_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mps_stats_su2u1pg" "SU2U1PG" "mps_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")

# create_tools_symm_target("mpo_stats_u1"      "U1"      "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mpo_stats_2u1"     "TWOU1"   "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mpo_stats_2u1pg"   "TWOU1PG" "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mpo_stats_u1dg"    "U1DG"    "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mpo_stats_su2u1"   "SU2U1"   "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("mpo_stats_su2u1pg" "SU2U1PG" "mpo_stats.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")

# create_tools_symm_target("visual_mpo_u1"      "U1"      "visual_mpo.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("visual_mpo_2u1"     "TWOU1"   "visual_mpo.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")

#create_tools_symm_target("mps_upgrade_u1"   "U1"     "mps_upgrade.cpp"  "${DMRG_APP_LIBRARIES}")
#create_tools_symm_target("mps_upgrade_2u1"  "TWOU1"  "mps_upgrade.cpp"  "${DMRG_APP_LIBRARIES}")
#create_tools_symm_target("mps_upgrade_none" "None"   "mps_upgrade.cpp"  "${DMRG_APP_LIBRARIES}")

if (BUILD_ONETDM_DIAG)
  create_tools_symm_target("onetdm_diag_2u1"     "TWOU1"   "onetdm_diag.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("onetdm_diag_2u1pg"  "TWOU1PG"  "onetdm_diag.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
endif(BUILD_ONETDM_DIAG)

# create_tools_symm_target("nrdmexpect_2u1pg"   "TWOU1PG" "nrdmexpect.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
# create_tools_symm_target("nrdmexpect_su2u1pg" "SU2U1PG" "nrdmexpect.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")

if (BUILD_MPS2CI)
  create_tools_symm_target("mps2ci_2u1"           "TWOU1"     "mps2ci.cpp"          "dmrg_models;${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps2ci_2u1pg"         "TWOU1PG"   "mps2ci.cpp"          "dmrg_models;${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps2ci_2u1_complex"   "TWOU1"     "mps2ci_complex.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps2ci_2u1pg_complex" "TWOU1PG"   "mps2ci_complex.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("srcas_2u1pg"          "TWOU1PG"   "srcas.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
  create_tools_vib_target( "mps2ci_vib"           "mps2ci_vibrational.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
  create_tools_vib_target( "srcas_vib"            "srcas_vibrational.cpp"   "dmrg_models;${DMRG_APP_LIBRARIES}")# ;srcas_utils")
endif (BUILD_MPS2CI)

if (BUILD_SRCAS)
  create_tools_target("srcas"  "srcas.cpp"   "dmrg_models;${DMRG_APP_LIBRARIES}")
endif (BUILD_SRCAS)
# DEAS
if (BUILD_DEAS)
  create_tools_symm_target("deas"   "TWOU1"   "deas.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
endif(BUILD_DEAS)

# GUTZWILLER CORRELATOR APPLICATION
if (BUILD_GUTZWILLER)
  create_tools_symm_target("mps_gutzwiller_complex_2u1"     "TWOU1"   "mps_gutzwiller_complex.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
  create_tools_symm_target("mps_gutzwiller_complex_2u1pg"   "TWOU1PG" "mps_gutzwiller_complex.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
endif(BUILD_GUTZWILLER)

# MPS_TRANSFORM
if (BUILD_MPS_TRANSFORM)
  if(DMRG_HAS_TwoU1)
    create_tools_symm_target("mps_transform"    "SU2U1"   "mps_transform.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
  endif(DMRG_HAS_TwoU1)
  if(DMRG_HAS_TwoU1PG)
    create_tools_symm_target("mps_transform_pg" "SU2U1PG" "mps_transform.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
  endif(DMRG_HAS_TwoU1PG)
endif (BUILD_MPS_TRANSFORM)
