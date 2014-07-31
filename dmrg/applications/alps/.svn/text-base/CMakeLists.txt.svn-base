
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# *** Libraries
add_library(mps_pscan STATIC libpscan/scheduler.cpp libpscan/options.cpp)

# *** Targets
set(MYAPP_LIBRARIES mps_pscan dmrg_models dmrg_utils ${DMRG_LIBRARIES})
add_subdirectory(mps_eigenstate_sims)
add_subdirectory(mps_evolve)
