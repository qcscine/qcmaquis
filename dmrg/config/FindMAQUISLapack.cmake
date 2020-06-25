
if(NOT DEFINED MAQUISLapack_LIBRARIES)
  set(MAQUISLapack_LIBRARIES)
endif()
if(NOT DEFINED MAQUISLapack_INCLUDE_DIRS)
  set(MAQUISLapack_INCLUDE_DIRS)
endif()
if(NOT DEFINED MAQUISLapack_LIB_DIRS)
  set(MAQUISLapack_LIB_DIRS)
endif()
if(NOT DEFINED MAQUISLapack_CXX_COMPILER_FLAGS)
  set(MAQUISLapack_CXX_COMPILER_FLAGS)
endif()
if(NOT DEFINED MAQUISLapack_LINKER_FLAGS)
  set(MAQUISLapack_LINKER_FLAGS)
endif()


if(${BLAS_LAPACK_SELECTOR} MATCHES "mkl_sequential")
  if(NOT DEFINED ENV{MKLROOT})
    message(FATAL_ERROR "ENV variable MKLROOT is required for BLAS_LAPACK_SELECTOR in MKL mode.")
  endif(NOT DEFINED ENV{MKLROOT})

  if(NOT APPLE)
    set(MAQUISLapack_LIB_DIRS $ENV{MKLROOT}/lib/intel64)
  else()
    set(MAQUISLapack_LIB_DIRS $ENV{MKLROOT}/lib)
  endif()

  set (MKL_LIB "")
  if (LAPACK_64_BIT)
    message(STATUS "Enabling 64 bit integers with Intel MKL")
    set (MKL_LIB "-lmkl_intel_ilp64")
  else ()
    message(STATUS "Disabled 64 bit integers with Intel MKL")
    set (MKL_LIB "-lmkl_intel_lp64")
  endif ()
  set(MAQUISLapack_LIBRARIES "${MKL_LIB} -lmkl_core -lmkl_sequential -lpthread -lm")
  set(MAQUISLapack_INCLUDE_DIRS $ENV{MKLROOT}/include)


elseif(${BLAS_LAPACK_SELECTOR} MATCHES "mkl_parallel")
  if(NOT DEFINED ENV{MKLROOT})
    message(FATAL_ERROR "ENV variable MKLROOT is required for BLAS_LAPACK_SELECTOR in MKL mode.")
  endif(NOT DEFINED ENV{MKLROOT})

  if(NOT APPLE)
    set(MAQUISLapack_LIB_DIRS $ENV{MKLROOT}/lib/intel64)
  else()
    set(MAQUISLapack_LIB_DIRS $ENV{MKLROOT}/lib)
  endif()
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    set(MAQUISLapack_LIBRARIES "${MKL_LIB} -lmkl_core -lmkl_intel_thread -lpthread -lm -openmp")
  else()
    set(MAQUISLapack_LIBRARIES "${MKL_LIB} -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -fopenmp")
  endif()
  set(MAQUISLapack_INCLUDE_DIRS $ENV{MKLROOT}/include)


# elseif(${BLAS_LAPACK_SELECTOR} MATCHES "alps")
#   set(MAQUISLapack_LIBRARIES ${ALPS_LAPACK_LIBRARIES} ${ALPS_LAPACK_LIBRARY} ${ALPS_BLAS_LIBRARIES} ${ALPS_BLAS_LIBRARY})
#   set(MAQUISLapack_CXX_COMPILER_FLAGS ${ALPS_LAPACK_DEFINITIONS})
#   set(MAQUISLapack_LINKER_FLAGS ${ALPS_LAPACK_LINKER_FLAGS})

elseif(${BLAS_LAPACK_SELECTOR} MATCHES "veclib")
  if(APPLE)
    set(MAQUISLapack_CXX_COMPILER_FLAGS "-framework Accelerate")
    set(MAQUISLapack_LIBRARIES "/usr/lib/libpthread.dylib")
    set(MAQUISLapack_LINKER_FLAGS "-framework Accelerate")
  else(APPLE)
    message(FATAL_ERROR "VecLib available only on Mac.")
  endif(APPLE)

elseif(${BLAS_LAPACK_SELECTOR} MATCHES "openblas")
set(MAQUISLapack_LIB_DIRS ${OPENBLASROOT}/lib)
set(MAQUISLapack_INCLUDE_DIRS ${OPENBLASROOT}/include)

if (APPLE)
    set(MAQUISLapack_LIBRARIES "${MAQUISLapack_LIB_DIRS}/libopenblas.dylib")
else(APPLE)
    set(MAQUISLapack_LIBRARIES "${MAQUISLapack_LIB_DIRS}/libopenblas.so")
endif(APPLE)
## For some reason find_library IGNORES NO_DEFAULT_PATH and includes a wrong path, so we will not use it
# find_library(MAQUISLapack_LIBRARIES NAMES openblas PATHS ${MAQUISLapack_LIB_DIRS}
#              NO_DEFAULT_PATH
#              NO_CMAKE_PATH
#              NO_CMAKE_FIND_ROOT_PATH
#              NO_SYSTEM_ENVIRONMENT_PATH
#              NO_CMAKE_SYSTEM_PATH)

  if (OPENBLASROOT STREQUAL "")
    set (OPENBLASROOT $ENV{OPENBLASROOT} CACHE PATH "OpenBLAS root directory." FORCE)
    if (NOT OPENBLASROOT)
        set (OPENBLASROOT "${OPENBLASROOT}")
        message ("-- Default OPENBLASROOT: ${OPENBLASROOT}")
    else()
        # at this point, OPENBLASROOT should be defined and not empty
        message ("-- OPENBLASROOT = ${OPENBLASROOT}")
    endif ()
  endif ()

  message(STATUS "LAPACK 64-bit integers:" ${LAPACK_64_BIT})
  # check the 64-bit version of OpenBLAS
  if (LAPACK_64_BIT)
    # check if the configuration header is available
    set (OPENBLAS_CONFIG "${OPENBLASROOT}/include/openblas_config.h")
    if (NOT EXISTS ${OPENBLAS_CONFIG})
         # for system-wide OpenBLAS installations the config path might be different
         # so try with an alternative path
         set (OPENBLAS_CONFIG "${OPENBLASROOT}/include/openblas/openblas_config.h")
    endif ()
    if (NOT EXISTS ${OPENBLAS_CONFIG})
            message (FATAL_ERROR
                    "Could not find OpenBLAS config header in: ${OPENBLASROOT} "
                    "(tried ${OPENBLASROOT}/include/openblas_config.h and "
                    " ${OPENBLASROOT}/include/openblas/openblas_config.h) "
                    "Please check if the OPENBLASROOT variable points to a "
                    "valid OpenBLAS installation directory."
                    )
    endif ()
    # check if the OpenBLAS installation was configured with 64bit integer support
    message ("-- Checking OpenBLAS for 64-bit integer interface...")
    include(CheckSymbolExists)
    unset (OPENBLAS_WITH_ILP64 CACHE)
    check_symbol_exists("OPENBLAS_USE64BITINT" ${OPENBLAS_CONFIG} OPENBLAS_WITH_ILP64)
    if (ADDRMODE EQUAL 64 AND NOT OPENBLAS_WITH_ILP64)
            message (FATAL_ERROR
                    "OpenBLAS was not configured for 64-bit integer interface, "
                    "please build OpenBLAS with INTERFACE64=1 defined."
                    )
    endif ()
  endif()

elseif(${BLAS_LAPACK_SELECTOR} MATCHES "manual")
  # variables set manually

else() # auto mode
  message ("No linear algebra library provided, trying to guess... Make sure you compile with the correct integer interface.")
  find_package(LAPACK)
  if(LAPACK_FOUND)
    set(MAQUISLapack_LIBRARIES ${LAPACK_LIBRARIES})
  endif(LAPACK_FOUND)

endif()


# include this to handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
# find_package_handle_standard_args(MAQUISLapack FOUND_VAR MAQUISLapack_FOUND
#                                                REQUIRED_VARS MAQUISLapack_LIBRARIES)
find_package_handle_standard_args(MAQUISLapack REQUIRED_VARS MAQUISLapack_LIBRARIES)
