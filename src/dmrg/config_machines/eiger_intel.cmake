
set(CMAKE_C_COMPILER icc)
set(CMAKE_CXX_COMPILER icpc)
#fill up if pb

# This module will set the following variables per language in your project,
# where <lang> is one of C, CXX, or Fortran:
#   MPI_<lang>_FOUND           TRUE if FindMPI found MPI flags for <lang>
#   MPI_<lang>_COMPILER        MPI Compiler wrapper for <lang>
#   MPI_<lang>_COMPILE_FLAGS   Compilation flags for MPI programs
#   MPI_<lang>_INCLUDE_PATH    Include path(s) for MPI header
#   MPI_<lang>_LINK_FLAGS      Linking flags for MPI programs
#   MPI_<lang>_LIBRARIES       All libraries to link MPI programs against
# Additionally, FindMPI sets the following variables for running MPI
# programs from the command line:
#   MPIEXEC                    Executable for running MPI programs
#   MPIEXEC_NUMPROC_FLAG       Flag to pass to MPIEXEC before giving
#                              it the number of processors to run on
#   MPIEXEC_PREFLAGS           Flags to pass to MPIEXEC directly
#                              before the executable to run.
#   MPIEXEC_POSTFLAGS          Flags to pass to MPIEXEC after other flags

# I just set up the CXX FLAGS
set(MPI_CXX_FOUND TRUE)
set(MPI_CXX_INCLUDE_PATH "/apps/eiger/Intel-MPI-4.0/intel64/include")
set(MPI_CXX_LINK_FLAGS "")
set(MPI_CXX_LIBRARIES "-L/apps/eiger/Intel-MPI-4.0/lib64 -lmpi_mt")
set(MPIEXEC "/apps/eiger/Intel-MPI-4.0/intel64/bin/mpirun")
set(MPIEXEC_NUMPROC_FLAG "-np")
set(MPIEXEC_MAX_NUMPROCS "2" CACHE STRING "Maximum number of processors available to run MPI applications.")
set(MPIEXEC_POSTFLAGS )

#BLAS - LAPACK
set(DEFAULT_BLAS_LAPACK manual)
set(BLAS_LAPACK_MANUAL_LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread)
set(BLAS_LAPACK_MANUAL_LIBS_DIR /apps/eiger/Intel-Composer-XE-2011/mkl/lib/intel64)
set(BLAS_LAPACK_MANUAL_INCLUDES /apps/eiger/Intel-Composer-XE-2011/mkl/include)

