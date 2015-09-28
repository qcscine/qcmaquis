add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_models dmrg_utils ${DMRG_LIBRARIES})
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})


macro(create_init_target TNAME SYMM CPPFILE LIBS)
  get_symm_group_name(SYMM_NAME ${SYMM})
  if(DMRG_HAS_${SYMM_NAME})
    add_executable(${TNAME} ${CPPFILE})
    target_link_libraries(${TNAME} ${LIBS})
    set_target_properties(${TNAME} PROPERTIES COMPILE_DEFINITIONS "USE_${SYMM}")
    install(TARGETS ${TNAME} RUNTIME DESTINATION bin COMPONENT applications)
  endif(DMRG_HAS_${SYMM_NAME})
endmacro(create_init_target)


#create_init_target("ci_generator_2u1"     "TWOU1"   "ci_generator.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
#create_init_target("ci_generator_2u1pg"   "TWOU1PG" "ci_generator.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")


create_init_target("det2mps_2u1"     "TWOU1"   "det2mps.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
create_init_target("det2mps_2u1pg"   "TWOU1PG" "det2mps.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
create_init_target("det2mps_su2u1"   "SU2U1"   "det2mps.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")
create_init_target("det2mps_su2u1pg" "SU2U1PG" "det2mps.cpp"  "dmrg_models;${DMRG_APP_LIBRARIES}")

# *** Configure symmetries
configure_symm_file("simulation_symm/simulation_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND DMRGINIT_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp")


# *** Targets
add_executable(dmrginit ${DMRGINIT_SYMM_SOURCES} dmrginit.cpp)
target_link_libraries(dmrginit ${DMRG_APP_LIBRARIES})

# *** Install
install(TARGETS dmrginit RUNTIME DESTINATION bin COMPONENT applications)
