# Client maintainer: dolfim@phys.ethz.ch
set(CTEST_SITE "ironhide.ethz.ch")
set(CTEST_BUILD_NAME "Darwin-Clang-MKL")
set(CTEST_BUILD_CONFIGURATION Debug)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_USE_LAUNCHERS 1)
set(CTEST_TEST_ARGS PARALLEL_LEVEL 2)
set(CTEST_BUILD_FLAGS -j6)

set(dashboard_source_name MAQUIS_DMRG-${CTEST_BUILD_NAME})

set(ALPS_ROOT $ENV{HOME}/opt/alps-clang)

set(ENV{CC}  clang)   # C compiler
set(ENV{CXX} clang++)  # C++ compiler
set(ENV{PATH} /opt/intel/lib:${ALPS_ROOT}/bin:$ENV{PATH})
set(ENV{LD_LIBRARY_PATH} /opt/intel/mkl/lib:/opt/intel/lib:${ALPS_ROOT}/lib)
set(ENV{DYLD_LIBRARY_PATH} /opt/intel/mkl/lib:/opt/intel/lib:${ALPS_ROOT}/lib)
set(ENV{MKL_NUM_THREADS} 4)

message(ENV IS: $ENV{PATH})

SET (dashboard_cache "
    ALPS_ROOT_DIR=${ALPS_ROOT}
    BUILD_DMRG=ON
    BUILD_TIME_EVOLUTION=ON
    BUILD_MULTIGRID=ON
    BUILD_CUSTOMAPP=ON
    BUILD_DIAG=ON
    BUILD_REGRESSION=ON
    ENABLE_APPLICATION_TESTS=ON
    ENABLE_OMP=OFF
    BLAS_LAPACK_SELECTOR=manual
    BLAS_LAPACK_MANUAL_LIBS_DIR=/opt/intel/mkl/lib;/opt/intel/lib
    BLAS_LAPACK_MANUAL_LIBS:STRING=-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
    CMAKE_INSTALL_PREFIX=/dev/null
")

include(${CTEST_SCRIPT_DIRECTORY}/dmrg_testing.cmake)

