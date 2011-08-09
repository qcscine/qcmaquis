
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()

set(DMRG_APP_LIB_SOURCES ${PROJECT_SOURCE_DIR}/utils/utils.cpp models.cpp)
add_library(dmrg_app STATIC ${DMRG_APP_LIB_SOURCES})

set(DMRG_APP_LIBRARIES dmrg_app ${DMRG_LIBRARIES})

add_executable(dmrg_u1 dmrg.cpp)
target_link_libraries(dmrg_u1 ${DMRG_APP_LIBRARIES})

add_executable(dmrg_2u1 dmrg.cpp)
target_link_libraries(dmrg_2u1 ${DMRG_APP_LIBRARIES})

if(BUILD_MULTIGRID)
	add_executable(multigrid_dmrg_u1 multigrid.cpp)
	target_link_libraries(multigrid_dmrg_u1 ${DMRG_APP_LIBRARIES})
	add_executable(multigrid_dmrg_2u1 multigrid.cpp)
	target_link_libraries(multigrid_dmrg_2u1 ${DMRG_APP_LIBRARIES})
endif(BUILD_MULTIGRID)

add_executable(measure_u1 dmrg.cpp)
target_link_libraries(measure_u1 ${DMRG_APP_LIBRARIES})

add_executable(measure_2u1 dmrg.cpp)
target_link_libraries(measure_2u1 ${DMRG_APP_LIBRARIES})

if(BUILD_MULTIGRID)
	add_executable(multigrid_measure_u1 ${PROJECT_SOURCE_DIR}/utils/utils.cpp multigrid.cpp)
	target_link_libraries(multigrid_measure_u1 ${DMRG_APP_LIBRARIES})
	add_executable(multigrid_measure_2u1 ${PROJECT_SOURCE_DIR}/utils/utils.cpp multigrid.cpp)
	target_link_libraries(multigrid_measure_2u1 ${DMRG_APP_LIBRARIES})
endif(BUILD_MULTIGRID)

add_executable(te_nn_u1 te_nn.cpp)
target_link_libraries(te_nn_u1 ${DMRG_APP_LIBRARIES})

add_executable(te_nn_2u1 te_nn.cpp)
target_link_libraries(te_nn_2u1 ${DMRG_APP_LIBRARIES})

add_executable(te_u1 te.cpp)
target_link_libraries(te_u1 ${DMRG_APP_LIBRARIES})

add_executable(te_2u1 te.cpp)
target_link_libraries(te_2u1 ${DMRG_APP_LIBRARIES})

set_target_properties(measure_u1 measure_2u1 PROPERTIES COMPILE_DEFINITIONS MEASURE_ONLY)
set_target_properties(dmrg_2u1 measure_2u1 te_nn_2u1 te_2u1 PROPERTIES COMPILE_DEFINITIONS UseTwoU1)

install(TARGETS dmrg_u1 dmrg_2u1 te_nn_u1 te_nn_2u1 te_u1 te_2u1 measure_u1 measure_2u1 RUNTIME DESTINATION bin COMPONENT applications)

if(BUILD_MULTIGRID)
	set_target_properties(multigrid_measure_u1 multigrid_measure_2u1 PROPERTIES COMPILE_DEFINITIONS MEASURE_ONLY)
	set_target_properties(multigrid_dmrg_2u1 multigrid_measure_2u1 PROPERTIES COMPILE_DEFINITIONS UseTwoU1)

	install(TARGETS multigrid_dmrg_u1 multigrid_dmrg_2u1 multigrid_measure_u1 multigrid_measure_2u1 RUNTIME DESTINATION bin COMPONENT applications)
endif(BUILD_MULTIGRID)