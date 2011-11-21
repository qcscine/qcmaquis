

include_directories ("${P_DENSE_MATRIX_DIR}")
include_directories ("${AMBIENT_DIR}")
#for boost numeric it is not an offical part of boost
#so set up by hand !!!!!!!!
#include_directories ("/opt/alps/include/")
include_directories ("/users/ewartt/ALPS2/include/")

if(BLAS_Accelerate_LIBRARY)
    add_definitions(-DAPPLE -DMPI_PARALLEL -framework Accelerate)
    add_executable (main.out main.cpp)
    target_link_libraries (main.out ambient ${BLAS_vecLib_LIBRARY} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY}  )
endif()

if(BLAS_mkl_core_LIBRARY)
#    find_path(BLAS_mkl_INCLUDE_DIR mkl.h ${BLAS_mkl_core_LIBRARY}/../../../include)
#    include_directories("${BLAS_mkL_INCLUDE_DIR}")
    add_definitions(-DMPI_PARALLEL)
    add_executable (main.out main.cpp)
    if(BLAS_MULTITHREADED)
        target_link_libraries (main.out ambient ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_intel_thread_LIBRARY} iomp5  pthread  ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    else()
        find_library(BLAS_mkl_sequential_LIBRARY mkl_sequential ${BLAS_mkl_core_LIBRARY}/../../../lib/em64t/)
        target_link_libraries (main.out ambient ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY}  pthread  ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY} )
    endif()

endif()



