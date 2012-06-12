
set(APPTEST_NAMES bh ff_middle fermiladder)
set(APPTEST_PREFIX dmrg)
set(APPTEST_SOURCEDIR ${CMAKE_CURRENT_SOURCE_DIR})
set(APPTEST_RUNDIR ${CMAKE_CURRENT_BINARY_DIR})
set(APPTEST_PYEXEC ${ALPS_ROOT_DIR}/bin/alpspython)
get_target_property(APPTEST_DMRG1EXEC dmrg LOCATION)
get_target_property(APPTEST_DMRG2EXEC dmrg_meas LOCATION)

message(STATUS "Using ${APPTEST_DMRG1EXEC} for dmrg exec")
message(STATUS "Using ${APPTEST_DMRG2EXEC} for dmrg_meas exec")

foreach(APPTEST_T ${APPTEST_NAMES})
  configure_file(${APPTEST_T}.py.in ${APPTEST_RUNDIR}/${APPTEST_T}.py)
  configure_file(apptest_${APPTEST_T}.py.in ${APPTEST_RUNDIR}/apptest_${APPTEST_T}.py)
  configure_file(${APPTEST_T}.testin.xml.in ${APPTEST_RUNDIR}/${APPTEST_T}.testin.xml)
  add_test(NAME "${APPTEST_PREFIX}_${APPTEST_T}" WORKING_DIRECTORY ${APPTEST_RUNDIR} COMMAND ${APPTEST_PYEXEC} apptest_${APPTEST_T}.py)
endforeach()
