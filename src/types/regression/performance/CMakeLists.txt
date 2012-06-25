compile_with_ambient()
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${AMBIENT_SOURCES_DIR}")

add_definitions(-DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
add_definitions(-DAMBIENT)

# Note: we can't mix openmp and pthreads but for performance tests we need two implementations:
# ambient gemm + sequential blas || mkl gemm + multithreaded blas

if(BLAS_mkl_core_LIBRARY) # intel
    include_directories("${MPI_CXX_INCLUDE_PATH}") # cluster path
    set(BLAS_LAPACK_multi     ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_intel_thread_LIBRARY} iomp5 ${CMAKE_THREAD_LIBS_INIT} pthread) 
    set(BLAS_LAPACK_serial    ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY}) 
    set(BLAS_SCALAPACK_serial ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY} mkl_blacs_intelmpi_lp64 mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64 )
endif()

if(BLAS_Accelerate_LIBRARY) # apple
    include_directories("${MPI_INCLUDE_PATH}") # cluster path
    set(MPI_CXX_LIBRARIES ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    add_definitions(-DAPPLE -DAMBIENT -framework Accelerate)
    set(BLAS_LAPACK ${BLAS_vecLib_LIBRARY}) 
endif()

add_executable( gemm.blas.out    gemm.blas.cpp    )
add_executable( gemm.pblas.out   gemm.pblas.cpp   )
add_executable( gemm.ambient.out gemm.ambient.cpp )

target_link_libraries( gemm.blas.out     ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_LAPACK_multi}          )
target_link_libraries( gemm.pblas.out    ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_SCALAPACK_serial}      )
target_link_libraries( gemm.ambient.out  ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_LAPACK_serial} ambient )

add_test( gemm.blas.out    gemm.blas.out    )
add_test( gemm.pblas.out   gemm.pblas.out   )
add_test( gemm.ambient.out gemm.ambient.out )
