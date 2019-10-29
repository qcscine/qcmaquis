
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)


# *** List all needed sources
## Lattice library
set(DMRG_MPIInterface_SOURCES "mpi_interface.cpp")
add_library(dmrg_mpiinterface STATIC ${DMRG_MPIInterface_SOURCES})

set_property(TARGET dmrg_mpiinterface PROPERTY POSITION_INDEPENDENT_CODE TRUE)

# *** Install

install(TARGETS dmrg_mpiinterface
        EXPORT DMRGTargets COMPONENT libraries DESTINATION lib)
export(TARGETS dmrg_mpiinterface
       FILE "${PROJECT_BINARY_DIR}/DMRGTargets.cmake")
