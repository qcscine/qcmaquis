set(APPTEST_PREFIX tevol_)
get_target_property(APPTEST_TEVOLAPP dmrg_tevol LOCATION)
get_target_property(APPTEST_MEASAPP dmrg_cmeas LOCATION)

message(STATUS "Using ${APPTEST_TEVOLAPP} for dmrg_tevol exec")
message(STATUS "Using ${APPTEST_MEASAPP} for dmrg_cmeas exec")

file(GLOB apptests "test_*py")
foreach(tt ${apptests})
    file(RELATIVE_PATH test_name ${CMAKE_CURRENT_SOURCE_DIR} ${tt})
    set(test_name ${APPTEST_PREFIX}${test_name})
    add_test(NAME ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${tt} ${APPTEST_TEVOLAPP} ${APPTEST_MEASAPP})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_SOURCE_DIR}/lib/python:$ENV{PYTHONPATH})
endforeach()
