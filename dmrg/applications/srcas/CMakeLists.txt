add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils ${DMRG_LIBRARIES})

##if(USE_AMBIENT)
##  compile_with_ambient()
##  list(INSERT DMRG_APP_LIBRARIES 0 ambient ${MPI_CXX_LIBRARIES})
##endif(USE_AMBIENT)

add_executable(mps_inspection mps_inspection.cpp)
target_link_libraries(mps_inspection ${DMRG_APP_LIBRARIES})
#set_target_properties(mps_inspection PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")

add_executable(mps2ci mps2ci.cpp)
target_link_libraries(mps2ci ${DMRG_APP_LIBRARIES})

add_executable(srcas srcas.cpp)
target_link_libraries(srcas ${DMRG_APP_LIBRARIES})
