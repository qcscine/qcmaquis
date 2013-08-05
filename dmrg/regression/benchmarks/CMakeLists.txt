add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})

if(USE_AMBIENT)
  compile_with_ambient()
  list(INSERT DMRG_APP_LIBRARIES 0 ambient ${MPI_CXX_LIBRARIES})
endif(USE_AMBIENT)  


# *** Targets
add_executable(ss_optim_u1 ss_optim.cpp)
target_link_libraries(ss_optim_u1 ${DMRG_APP_LIBRARIES})

add_executable(ss_optim_2u1 ss_optim.cpp)
target_link_libraries(ss_optim_2u1 ${DMRG_APP_LIBRARIES})
set_target_properties(ss_optim_2u1 PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")


add_executable(ts_optim_u1 ts_optim.cpp)
target_link_libraries(ts_optim_u1 ${DMRG_APP_LIBRARIES})

add_executable(ts_optim_2u1 ts_optim.cpp)
target_link_libraries(ts_optim_2u1 ${DMRG_APP_LIBRARIES})
set_target_properties(ts_optim_2u1 PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")
