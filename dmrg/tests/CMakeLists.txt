add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES maquis_dmrg dmrg_models dmrg_utils ${DMRG_LIBRARIES})

# *** Targets
add_executable(test1 test1.cpp)
target_link_libraries(test1 ${DMRG_APP_LIBRARIES})

add_executable(test2 test2.cpp)
target_link_libraries(test2 ${DMRG_APP_LIBRARIES})

if(BUILD_DMRG_EVOLVE)
    add_executable(test_imaginarytime test_imaginarytime.cpp)
    target_link_libraries(test_imaginarytime ${DMRG_APP_LIBRARIES})
    add_executable(test_realtime test_realtime.cpp)
    target_link_libraries(test_realtime ${DMRG_APP_LIBRARIES})
endif(BUILD_DMRG_EVOLVE)

if(BUILD_PREBO)
    add_executable(prebo_test1 prebo_test1.cpp)
    target_link_libraries(prebo_test1 ${DMRG_APP_LIBRARIES})
    add_executable(prebo_test2 prebo_test2.cpp)
    target_link_libraries(prebo_test2 ${DMRG_APP_LIBRARIES})
    add_executable(prebo_test_model model/model_prebo.cpp)
    target_link_libraries(prebo_test_model ${DMRG_APP_LIBRARIES})
    # Unit tests
    add_executable(test_prebo_symbolic_operator prebo/symbolicoperator_test.cpp)
    target_link_libraries(test_prebo_symbolic_operator ${DMRG_APP_LIBRARIES})
    add_executable(test_prebo_symbolic_jw prebo/SymbolicJordanWignerTest.cpp)
    target_link_libraries(test_prebo_symbolic_jw ${DMRG_APP_LIBRARIES})
endif(BUILD_PREBO)

add_executable(test_integral_map test_integral_map.cpp)
target_link_libraries(test_integral_map ${DMRG_APP_LIBRARIES})

add_executable(test_rel test_rel.cpp)
target_link_libraries(test_rel ${DMRG_APP_LIBRARIES})

add_executable(test_hirdm test_hirdm.cpp)
target_link_libraries(test_hirdm ${DMRG_APP_LIBRARIES})

add_executable(test_mps_mpo_ops test_mps_mpo_ops/test_mps_mpo_ops.cpp)
target_link_libraries(test_mps_mpo_ops ${DMRG_APP_LIBRARIES})

add_executable(test_siteproblem test_mps_mpo_ops/test_siteproblem.cpp)
target_link_libraries(test_siteproblem ${DMRG_APP_LIBRARIES})

add_executable(test_mpsjoin test_mps_mpo_ops/mpsjoin.cpp)
target_link_libraries(test_mpsjoin ${DMRG_APP_LIBRARIES})

add_executable(test_wigner test_wigner.cpp)
target_link_libraries(test_wigner ${DMRG_APP_LIBRARIES})

# -- Unit tests --
add_executable(test_block_matrix block_matrix/block_matrix.cpp)
target_link_libraries(test_block_matrix ${DMRG_APP_LIBRARIES})

if(BUILD_DMRG_EVOLVE)
    add_executable(test_time_evolvers TimeEvolvers/TimeEvolvers.cpp)
    target_link_libraries(test_time_evolvers ${DMRG_APP_LIBRARIES})

    add_executable(test_site_shifter TimeEvolvers/SiteShifter.cpp)
    target_link_libraries(test_site_shifter ${DMRG_APP_LIBRARIES})

    add_executable(test_time_evolution_sweep TimeEvolvers/TimeEvolutionSweep.cpp)
    target_link_libraries(test_time_evolution_sweep ${DMRG_APP_LIBRARIES})
endif(BUILD_DMRG_EVOLVE)