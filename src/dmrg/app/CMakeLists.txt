
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()

configure_file(version.h.in ${CMAKE_CURRENT_BINARY_DIR}/version.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# *** Libraries
if(USE_MTM)
  add_library(mtm STATIC ${PROJECT_SOURCE_DIR}/dense_matrix/mt_matrix/mtm_scheduler.cpp)
endif(USE_MTM)

add_library(dmrg_utils STATIC ${PROJECT_SOURCE_DIR}/utils/utils.cpp ${PROJECT_SOURCE_DIR}/utils/stream_storage.cpp)
add_library(dmrg_models STATIC 	models/models_null.cpp alps_model/alps_model_null.cpp
								models/models_u1.cpp alps_model/alps_model_u1.cpp
								models/models_2u1.cpp alps_model/alps_model_2u1.cpp)


# *** Programs
set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})
if(USE_MTM)
  list(APPEND ${DMRG_APP_LIBRARIES} mtm)
endif(USE_MTM)


add_executable(dmrg dmrg/dmrg_null.cpp dmrg/dmrg_u1.cpp dmrg/dmrg_2u1.cpp dmrg.cpp)
target_link_libraries(dmrg ${DMRG_APP_LIBRARIES})

add_executable(measure measure/measure_null.cpp measure/measure_u1.cpp measure/measure_2u1.cpp measure.cpp)
target_link_libraries(measure ${DMRG_APP_LIBRARIES})

if(BUILD_MULTIGRID)
	add_executable(multigrid_dmrg_u1 multigrid.cpp)
	target_link_libraries(multigrid_dmrg_u1 ${DMRG_APP_LIBRARIES})
	add_executable(multigrid_dmrg_2u1 multigrid.cpp)
	target_link_libraries(multigrid_dmrg_2u1 ${DMRG_APP_LIBRARIES})
endif(BUILD_MULTIGRID)

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

set_target_properties(te_nn_2u1 te_2u1 PROPERTIES COMPILE_DEFINITIONS UseTwoU1)

install(TARGETS dmrg te_nn_u1 te_nn_2u1 te_u1 te_2u1 measure RUNTIME DESTINATION bin COMPONENT applications)

if(BUILD_MULTIGRID)
	set_target_properties(multigrid_measure_u1 multigrid_measure_2u1 PROPERTIES COMPILE_DEFINITIONS MEASURE_ONLY)
	set_target_properties(multigrid_dmrg_2u1 multigrid_measure_2u1 PROPERTIES COMPILE_DEFINITIONS UseTwoU1)

	install(TARGETS multigrid_dmrg_u1 multigrid_dmrg_2u1 multigrid_measure_u1 multigrid_measure_2u1 RUNTIME DESTINATION bin COMPONENT applications)
endif(BUILD_MULTIGRID)