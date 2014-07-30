
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_models dmrg_utils ${DMRG_LIBRARIES})

if(USE_AMBIENT)
  compile_with_ambient()
else(USE_AMBIENT)
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
  include_directories(${MPI_CXX_INCLUDE_PATH})
  list(APPEND DMRG_APP_LIBRARIES ${MPI_LIBRARIES})
endif(USE_AMBIENT)  


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
