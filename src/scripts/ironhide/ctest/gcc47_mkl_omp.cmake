# Client maintainer: dolfim@phys.ethz.ch
set(CTEST_SITE "ironhide.ethz.ch")
set(CTEST_BUILD_NAME "Darwin-GCC47-MKL-OMP")
set(CTEST_BUILD_CONFIGURATION Debug)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_TEST_ARGS PARALLEL_LEVEL 2)
set(CTEST_BUILD_FLAGS -j6)

set(dashboard_source_name MAQUIS_DMRG-${CTEST_BUILD_NAME})

set(ALPS_ROOT $ENV{HOME}/opt/alps-gcc47)

set(ENV{CC}  /opt/local/bin/gcc-mp-4.7)   # C compiler
set(ENV{CXX} /opt/local/bin/g++-mp-4.7)  # C++ compiler
set(ENV{PATH} /opt/intel/lib:${ALPS_ROOT}/bin:$ENV{PATH})
set(ENV{LD_LIBRARY_PATH} /opt/intel/mkl/lib:/opt/intel/lib:${ALPS_ROOT}/lib)
set(ENV{DYLD_LIBRARY_PATH} /opt/intel/mkl/lib:/opt/intel/lib:${ALPS_ROOT}/lib)
set(ENV{OMP_NUM_THREADS} 4)

SET (dashboard_cache "
    ALPS_ROOT_DIR=${ALPS_ROOT}
    BUILD_DMRG=ON
    BUILD_TIME_EVOLUTION=ON
    BUILD_MULTIGRID=ON
    BUILD_CUSTOMAPP=ON
    BUILD_DIAG=ON
    BUILD_REGRESSION=ON
    ENABLE_APPLICATION_TESTS=ON
    ENABLE_OMP=ON
    BLAS_LAPACK_SELECTOR=manual
    BLAS_LAPACK_MANUAL_LIBS_DIR=/opt/intel/mkl/lib
    BLAS_LAPACK_MANUAL_LIBS:STRING=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
    CMAKE_INSTALL_PREFIX=/dev/null
")

include(${CTEST_SCRIPT_DIRECTORY}/dmrg_testing.cmake)

