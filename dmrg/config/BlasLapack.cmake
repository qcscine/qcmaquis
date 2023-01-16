# Detect BLAS/LAPACK vendor
macro (get_blas_vendor BLA_DIR)
    set (BLA_VENDOR "Generic")
    string(FIND ${BLA_DIR} "libmkl_intel" _POS)
    if (NOT _POS EQUAL -1)
        set(BLA_VENDOR "MKL")
    endif()
    string(FIND ${BLA_DIR} "openblas" _POS)
    if (NOT _POS EQUAL -1)
        set(BLA_VENDOR "openblas")
    endif()
    # TODO @stknecht: please add Accelerate detection accordingly
endmacro()

# FIXME: This seems to ignore custom BLAS_DIR and LAPACK_DIR and uses the default path
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
get_blas_vendor(${BLAS_LIBRARIES})
message(STATUS "BLAS/LAPACK vendor: ${BLA_VENDOR}")

# if MKL is found and 64-bit integer interface is required, search again for MKL, this time with the correct integer interface
if (BLA_VENDOR STREQUAL "MKL" AND LAPACK_64_BIT)
    cmake_minimum_required(VERSION 3.13)
    message(STATUS "Enabling 64-bit integer interface with Intel MKL")
    set(BLA_VENDOR "Intel10_64ilp_seq") # Warning: this forces linking of sequential mkl
    # retry finding MKL
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()

message(STATUS "BLAS library: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK library: ${LAPACK_LIBRARIES}")
# Autodetect ILP64 or LP64 integer interface
try_run(_IS_32_BIT_INT _COMPILES ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/config/ilp64.cpp
    LINK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
    COPY_FILE ip64test
)

if (NOT _IS_32_BIT_INT)
    if (LAPACK_64_BIT) # 64-bit integers requested but not supported
        message(FATAL_ERROR "Requested 64 bit integer BLAS but your BLAS library does not support it")
    else()
    message(STATUS "BLAS/LAPACK compiled with 64-bit integer interface")
        set(LAPACK_64_BIT ON)
    endif()
endif ()
