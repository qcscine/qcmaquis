
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_SERIAL_LIBRARIES dmrg_models dmrg_utils ${DMRG_LIBRARIES})



macro(create_pcorr_serial_symm_target TNAME SYMM)
  get_symm_group_name(SYMM_NAME ${SYMM})
  if(DMRG_HAS_${SYMM_NAME})
    add_executable(${TNAME} main_serial.cpp)
    target_link_libraries(${TNAME} ${DMRG_SERIAL_LIBRARIES})
    set_target_properties(${TNAME} PROPERTIES COMPILE_DEFINITIONS "USE_${SYMM}")
    install(TARGETS ${TNAME} RUNTIME DESTINATION bin COMPONENT applications)
  endif(DMRG_HAS_${SYMM_NAME})
endmacro(create_pcorr_serial_symm_target)

create_pcorr_serial_symm_target("pcorr_serial_nu1"  "NU1")
create_pcorr_serial_symm_target("pcorr_serial_2u1"  "TWOU1")
create_pcorr_serial_symm_target("pcorr_serial_u1"   "U1"   )
create_pcorr_serial_symm_target("pcorr_serial_none" "NONE" )



if(USE_AMBIENT)
  compile_with_ambient()
else(USE_AMBIENT)
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
  include_directories(${MPI_CXX_INCLUDE_PATH})
  set(DMRG_APP_LIBRARIES ${DMRG_SERIAL_LIBRARIES} ${MPI_LIBRARIES})
endif(USE_AMBIENT)  

macro(create_pcorr_symm_target TNAME SYMM CPPFILE)
  get_symm_group_name(SYMM_NAME ${SYMM})
  if(DMRG_HAS_${SYMM_NAME})
    add_executable(${TNAME} ${CPPFILE})
    target_link_libraries(${TNAME} ${DMRG_APP_LIBRARIES})
    set_target_properties(${TNAME} PROPERTIES COMPILE_DEFINITIONS "USE_${SYMM}")
    install(TARGETS ${TNAME} RUNTIME DESTINATION bin COMPONENT applications)
  endif(DMRG_HAS_${SYMM_NAME})
endmacro(create_pcorr_symm_target)


create_pcorr_symm_target("pcorr_nu1"  "NU1"   "main.cpp")
create_pcorr_symm_target("pcorr_2u1"  "TWOU1" "main.cpp")
create_pcorr_symm_target("pcorr_u1"   "U1"    "main.cpp")
create_pcorr_symm_target("pcorr_none" "NONE"  "main.cpp")


create_pcorr_symm_target("pcorr_nn_2u1"  "TWOU1" "main_nn.cpp")


