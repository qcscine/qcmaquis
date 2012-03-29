## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.

SET(CTEST_PROJECT_NAME "DMRG")
SET(CTEST_NIGHTLY_START_TIME "01:00:00 EST")

IF(NOT DEFINED CTEST_DROP_METHOD)
  SET(CTEST_DROP_METHOD "http")
ENDIF(NOT DEFINED CTEST_DROP_METHOD)

IF(CTEST_DROP_METHOD STREQUAL "http")
  SET(CTEST_DROP_SITE "maquis.ch")
  SET(CTEST_DROP_LOCATION "/cdash/submit.php?project=DMRG")
  SET(CTEST_DROP_SITE_CDASH TRUE)
# SET(CTEST_TRIGGER_SITE "") ## can be useful
ENDIF(CTEST_DROP_METHOD STREQUAL "http")
