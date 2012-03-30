# Client maintainer: alex.kosenkov@gmail.com
set(CTEST_SITE "castor")
set(CTEST_BUILD_NAME "Linux x64-icc-mvapich")
set(CTEST_BUILD_CONFIGURATION Debug)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_SOURCE_DIRECTORY "/users/maquis/maquis2012/src/types")
set(CTEST_BINARY_DIRECTORY "/users/maquis/maquis2012/src/types/build_intel_mvapich2")
#set(CTEST_CHECKOUT_COMMAND "svn up")
set(CTEST_UPDATE_COMMAND "svn")
set(CTEST_PROJECT_SUBPROJECTS types)

set(ENABLE_REGRESSION_VALIDATIONS ON)

set(subproject ${CTEST_PROJECT_SUBPROJECTS})
set_property(GLOBAL PROPERTY SubProject ${subproject})
set_property(GLOBAL PROPERTY Label ${subproject})

include(${CTEST_SCRIPT_DIRECTORY}/cmake_common.cmake)
