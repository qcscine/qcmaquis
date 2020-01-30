# *** List all needed sources

set(MPI_Interface_SOURCES "mpi_interface.cpp")

add_library(mpi_interface STATIC ${MPI_Interface_SOURCES})
add_library(mpi_interface DYNAMIC ${MPI_Interface_SOURCES})

set_property(TARGET mpi_interface PROPERTY POSITION_INDEPENDENT_CODE TRUE)

# *** Install

install(TARGETS mpi_interface
        EXPORT DMRGTargets COMPONENT libraries DESTINATION lib)

export(TARGETS mpi_interface
       FILE "${PROJECT_BINARY_DIR}/DMRGTargets.cmake")
