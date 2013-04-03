compile_with_ambient()
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${AMBIENT_SOURCES_DIR}")

add_definitions(-DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

# Note: we can't mix openmp and pthreads but for performance tests we need two implementations:
# ambient gemm + sequential blas || mkl gemm + multithreaded blas

if(BLAS_mkl_core_LIBRARY) # intel
    include_directories("${MPI_CXX_INCLUDE_PATH}") # cluster path
    set(BLAS_LAPACK_multi     ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_intel_thread_LIBRARY} iomp5 ${CMAKE_THREAD_LIBS_INIT} pthread) 
    set(BLAS_LAPACK_serial    ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY} ${PLASMA_LIBRARY} ) 
    set(BLAS_SCALAPACK_serial ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY} ${SCALAPACK_LIBRARY} )
endif()

if(BLAS_Accelerate_LIBRARY) # apple
    include_directories("${MPI_INCLUDE_PATH}") # cluster path
    set(MPI_CXX_LIBRARIES ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    add_definitions(-DAPPLE -framework Accelerate)
    set(BLAS_LAPACK ${BLAS_vecLib_LIBRARY}) 
endif()

add_executable( gemm.blas.out          gemm.blas.cpp          )
add_executable( gemm.pblas.out         gemm.pblas.cpp         )
add_executable( gemm.ambient.out       gemm.ambient.cpp       )
add_executable( gemm.tiled.ambient.out gemm.tiled.ambient.cpp )
add_executable( lq.ambient.out         lq.ambient.cpp         )
add_executable( lq.lapack.out          lq.lapack.cpp          )
add_executable( qr.ambient.out         qr.ambient.cpp         )
add_executable( qr.lapack.out          qr.lapack.cpp          )

target_link_libraries( gemm.pblas.out          ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_SCALAPACK_serial} ambient )
target_link_libraries( gemm.blas.out           ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_LAPACK_multi}     ambient )
target_link_libraries( qr.lapack.out           ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_LAPACK_multi}     ambient )
target_link_libraries( lq.lapack.out           ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_LAPACK_multi}     ambient )
target_link_libraries( gemm.ambient.out        ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_LAPACK_serial}    ambient )
target_link_libraries( gemm.tiled.ambient.out  ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_LAPACK_serial}    ambient )
target_link_libraries( lq.ambient.out          ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_LAPACK_serial}    ambient )
target_link_libraries( qr.ambient.out          ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES} ${BLAS_LAPACK_serial}    ambient )

add_test( gemm.blas.out          gemm.blas.out          )
add_test( gemm.ambient.out       gemm.ambient.out       )
add_test( gemm.tiled.ambient.out gemm.tiled.ambient.out )
add_test( lq.ambient.out         lq.ambient.out         )
add_test( lq.lapack.out          lq.lapack.out          )
add_test( qr.ambient.out         qr.ambient.out         )
add_test( qr.lapack.out          qr.lapack.out          )
