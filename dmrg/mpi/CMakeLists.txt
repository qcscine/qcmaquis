
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)


# *** List all needed sources
## Lattice library
set(DMRG_MPIInterface_SOURCES "mpi_interface.cpp")
add_library(mpi_interface STATIC ${DMRG_MPIInterface_SOURCES})

set_property(TARGET mpi_interface PROPERTY POSITION_INDEPENDENT_CODE TRUE)

# *** Install

install(TARGETS mpi_interface
        EXPORT DMRGTargets COMPONENT libraries DESTINATION lib)
export(TARGETS mpi_interface
       FILE "${PROJECT_BINARY_DIR}/DMRGTargets.cmake")
