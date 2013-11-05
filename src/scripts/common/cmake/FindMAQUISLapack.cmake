
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
  set(MAQUISLapack_LIBRARIES "-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm")
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
    set(MAQUISLapack_LIBRARIES "-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -openmp")
  else()
    set(MAQUISLapack_LIBRARIES "-lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -fopenmp")
  endif()
  set(MAQUISLapack_INCLUDE_DIRS $ENV{MKLROOT}/include)
  
  
elseif(${BLAS_LAPACK_SELECTOR} MATCHES "alps")
  set(MAQUISLapack_LIBRARIES "${ALPS_LAPACK_LIBRARIES} ${ALPS_LAPACK_LIBRARY} ${ALPS_BLAS_LIBRARIES} ${ALPS_BLAS_LIBRARY}")
  set(MAQUISLapack_CXX_COMPILER_FLAGS ${ALPS_LAPACK_DEFINITIONS})
  set(MAQUISLapack_LINKER_FLAGS "${ALPS_LAPACK_LINKER_FLAGS}")
 
elseif(${BLAS_LAPACK_SELECTOR} MATCHES "veclib")
  if(APPLE)
    set(MAQUISLapack_CXX_COMPILER_FLAGS "-framework vecLib")
    set(MAQUISLapack_LIBRARIES "/usr/lib/libpthread.dylib")
    set(MAQUISLapack_LINKER_FLAGS "-framework veclib")
  else(APPLE)
    meassage(FATAL_ERROR "VecLib available only on Mac.")
  endif(APPLE)
  
  
elseif(${BLAS_LAPACK_SELECTOR} MATCHES "manual")
  # variables set manually
  
else() # auto mode
  find_package(LAPACK)
  if(LAPACK_FOUND)
    set(MAQUISLapack_LIBRARIES LAPACK_LIBRARIES)
  endif(LAPACK_FOUND)
  
endif()


# include this to handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
# find_package_handle_standard_args(MAQUISLapack FOUND_VAR MAQUISLapack_FOUND
#                                                REQUIRED_VARS MAQUISLapack_LIBRARIES)
find_package_handle_standard_args(MAQUISLapack REQUIRED_VARS MAQUISLapack_LIBRARIES)
