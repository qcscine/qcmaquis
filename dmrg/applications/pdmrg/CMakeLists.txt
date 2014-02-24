add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)


find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
include_directories(${MPI_CXX_INCLUDE_PATH})

set(DMRG_APP_LIBRARIES dmrg_models dmrg_utils ${DMRG_LIBRARIES} ${MPI_LIBRARIES})


macro(create_pdrmg_symm_target TNAME SYMM)
  get_symm_group_name(SYMM_NAME ${SYMM})
  if(DMRG_HAS_${SYMM_NAME})
    add_executable(${TNAME} main.cpp)
    target_link_libraries(${TNAME} ${DMRG_APP_LIBRARIES})
    set_target_properties(${TNAME} PROPERTIES COMPILE_DEFINITIONS "USE_${SYMM}")
    install(TARGETS ${TNAME} RUNTIME DESTINATION bin COMPONENT applications)
  endif(DMRG_HAS_${SYMM_NAME})
endmacro(create_pdrmg_symm_target)


create_pdrmg_symm_target("pdmrg_2u1"  "TWOU1")
create_pdrmg_symm_target("pdmrg_u1"   "U1"   )
create_pdrmg_symm_target("pdmrg_none" "NONE" )

# *** Targets
# if(DMRG_HAS_TwoU1)
#   add_executable(pdmrg_2u1 main.cpp)
#   target_link_libraries(pdmrg_2u1 dmrg_models dmrg_utils ${DMRG_LIBRARIES})
#   set_target_properties(pdmrg_2u1 PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")
#   install(TARGETS pdmrg_2u1 RUNTIME DESTINATION bin COMPONENT applications)
# endif(DMRG_HAS_TwoU1)
# 
# if(DMRG_HAS_TrivialGroup)
#   add_executable(pdmrg_none main.cpp)
#   target_link_libraries(pdmrg_none ${DMRG_APP_LIBRARIES})
#   set_target_properties(pdmrg_none PROPERTIES COMPILE_DEFINITIONS "USE_NONE")
#   install(TARGETS pdmrg_none RUNTIME DESTINATION bin COMPONENT applications)
# endif(DMRG_HAS_TrivialGroup)
# 
# if(DMRG_HAS_U1)
#   add_executable(pdmrg_u1 main.cpp)
#   target_link_libraries(pdmrg_u1 ${DMRG_APP_LIBRARIES})
#   set_target_properties(pdmrg_u1 PROPERTIES COMPILE_DEFINITIONS "USE_U1")
#   install(TARGETS pdmrg_u1 RUNTIME DESTINATION bin COMPONENT applications)
# endif(DMRG_HAS_U1)
