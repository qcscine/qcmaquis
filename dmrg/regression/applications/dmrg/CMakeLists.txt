set(APPTEST_PREFIX dmrg_)
get_target_property(APPTEST_DMRGAPP dmrg LOCATION)
get_target_property(APPTEST_MEASAPP dmrg_meas LOCATION)

message(STATUS "Using ${APPTEST_DMRGAPP} for dmrg exec")
message(STATUS "Using ${APPTEST_MEASAPP} for dmrg_meas exec")

file(GLOB apptests "test_*py")
foreach(tt ${apptests})
    file(RELATIVE_PATH test_name ${CMAKE_CURRENT_SOURCE_DIR} ${tt})
    set(test_name ${APPTEST_PREFIX}${test_name})
    add_test(NAME ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${tt} ${APPTEST_DMRGAPP} ${APPTEST_MEASAPP})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/lib/python:$ENV{PYTHONPATH})
endforeach()
