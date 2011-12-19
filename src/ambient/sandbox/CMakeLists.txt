include_directories ("${P_DENSE_MATRIX_DIR}")
include_directories ("${AMBIENT_DIR}")
#for boost numeric it is not an offical part of boost
#so set up by hand !!!!!!!!
#ok by hand but I do not want a new dependency into ambient
include_directories ("/project/h07/ALPS_INTEL_EIGER/include/")
#Apple target
if(BLAS_Accelerate_LIBRARY)
    add_definitions(-DAPPLE -DMPI_PARALLEL -framework Accelerate)
    add_executable (main.out main.cpp)
    target_link_libraries (main.out ambient ${BLAS_vecLib_LIBRARY} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY}  )
endif()
#MKL target
if(BLAS_mkl_core_LIBRARY)
    add_definitions(-DMPI_PARALLEL)
    add_executable(main.out main.cpp)
    if(BLAS_MULTITHREADED)
        target_link_libraries (main.out ambient ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_intel_thread_LIBRARY} iomp5  ${CMAKE_THREAD_LIBS_INIT} ${MPI_CXX_LIBRARIES})
    else()
        target_link_libraries (main.out ambient ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY}  ${MPI_CXX_LIBRARIES})
    endif()
endif()
