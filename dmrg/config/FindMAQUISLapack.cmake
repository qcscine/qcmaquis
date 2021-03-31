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


if("${BLAS_LAPACK_SELECTOR}" MATCHES "mkl_sequential" OR "${BLAS_LAPACK_SELECTOR}" MATCHES "mkl_parallel")
    if(NOT DEFINED ENV{MKLROOT})
        if(DEFINED MKLROOT)
            set(ENV{MKLROOT} ${MKLROOT})
        else()
            message(FATAL_ERROR "ENV variable MKLROOT is required for BLAS_LAPACK_SELECTOR in MKL mode.")
        endif()
    endif()

    if(NOT APPLE)
        set(MAQUISLapack_LIB_DIRS $ENV{MKLROOT}/lib/intel64)
        # Fix for Ubuntu 20.04 where standard MKL installation does not set MKLROOT and the libraries are
        # found in /usr/lib/x86_64-linux-gnu
        if (NOT EXISTS ${MAQUISLapack_LIB_DIRS})
            set(MAQUISLapack_LIB_DIRS $ENV{MKLROOT}/lib/x86_64-linux-gnu)
            if (NOT EXISTS ${MAQUISLapack_LIB_DIRS})
                message(FATAL_ERROR "Could not determine path to MKL libraries. Did you set MKLROOT environment variable correctly?")
            endif()
        endif()
    else()
        set(MAQUISLapack_LIB_DIRS $ENV{MKLROOT}/lib)
    endif()

    # the code below has been adapted from OpenMOLCAS CMakeLists.txt, (C) Ignacio Fdez. Galvan
    set (MKL_ARCH "")

    if (LAPACK_64_BIT)
    message(STATUS "Enabled 64 bit integers with Intel MKL")
        set (MKL_ARCH "ilp64")
    else ()
        message(STATUS "Disabled 64 bit integers with Intel MKL")
        set (MKL_ARCH "lp64")
    endif ()

    set (MKL_COMPILER "")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(MKL_COMPILER "gf")
        else()
        set(MKL_COMPILER "intel")
        endif()

    # core library
    find_library (LIBMKL_CORE NAMES "mkl_core"
    PATHS ${MAQUISLapack_LIB_DIRS} NO_DEFAULT_PATH)

    # compiler-specific library interface
    find_library (LIBMKL_INTERFACE NAMES "mkl_${MKL_COMPILER}_${MKL_ARCH}"
            PATHS ${MAQUISLapack_LIB_DIRS} NO_DEFAULT_PATH)

    # sequential/compiler-specific threading interface
    find_library (LIBMKL_PARALLEL NAMES "${BLAS_LAPACK_SELECTOR}"
        PATHS ${MAQUISLapack_LIB_DIRS} NO_DEFAULT_PATH)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        find_library (LIBMKL_THREADING NAMES "mkl_gnu_thread"
                PATHS ${MAQUISLapack_LIB_DIRS} NO_DEFAULT_PATH)
    else()
        find_library (LIBMKL_THREADING NAMES "mkl_intel_thread"
            PATHS ${MAQUISLapack_LIB_DIRS} NO_DEFAULT_PATH)
    endif ()

    list (APPEND MKL_LIBRARIES ${LIBMKL_INTERFACE})
    list (APPEND MKL_LIBRARIES ${LIBMKL_CORE})
    list (APPEND MKL_LIBRARIES ${LIBMKL_PARALLEL})
    if (ENABLE_OMP)
        list (APPEND MKL_LIBRARIES ${LIBMKL_THREADING})
    endif()

  set(MAQUISLapack_LIBRARIES "${MKL_LIBRARIES}")
  set(MAQUISLapack_INCLUDE_DIRS $ENV{MKLROOT}/include)

elseif(${BLAS_LAPACK_SELECTOR} MATCHES "veclib")
  if(APPLE)
    set(MAQUISLapack_CXX_COMPILER_FLAGS "-framework Accelerate")
    if(CMAKE_HOST_SYSTEM_VERSION VERSION_LESS "20.0")
      set(MAQUISLapack_LIBRARIES "/usr/lib/libpthread.dylib")
    endif()
    set(MAQUISLapack_LINKER_FLAGS "-framework Accelerate")
  else(APPLE)
    message(FATAL_ERROR "VecLib available only on Mac.")
  endif(APPLE)

elseif(${BLAS_LAPACK_SELECTOR} MATCHES "openblas")

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

  find_library(MAQUISLapack_LIBRARIES NAMES openblas PATHS ${OPENBLASROOT} PATH_SUFFIXES lib lib64
                NO_DEFAULT_PATH)

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
