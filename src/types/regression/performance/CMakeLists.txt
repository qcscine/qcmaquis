compile_with_ambient()
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${AMBIENT_SOURCES_DIR}")

add_definitions (-DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
add_definitions(-DMPI_PARALLEL)

# we can not mix OpenMP and pthread into Ambient, but for performance tests I need two implementations : gemm_ambient serial blas, gemm_blas mulithread blas (MKL, ACML, ESSLSMP)
#
if(BLAS_mkl_core_LIBRARY) #intel
    include_directories("${MPI_CXX_INCLUDE_PATH}") #cluster path
    set(BLAS_LAPACK_multi  ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_intel_thread_LIBRARY} iomp5 ${CMAKE_THREAD_LIBS_INIT} pthread) 
    set(BLAS_LAPACK_serial ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY}) 
    set(BLAS_SCALAPACK_serial mkl_blacs_intelmpi_lp64 mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64  ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY})
endif()

if(BLAS_Accelerate_LIBRARY) #apple
    include_directories("${MPI_INCLUDE_PATH}") #cluster path
    set(MPI_CXX_LIBRARIES ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    add_definitions(-DAPPLE -DMPI_PARALLEL -framework Accelerate)
    set(BLAS_LAPACK ${BLAS_vecLib_LIBRARY}) 
endif()

add_executable (gemm_blas.out gemm_blas.cpp)
add_executable (gemm_ambient.out gemm_ambient.cpp)
add_executable (pdgemm_scalapack.out pdgemm_scalapack.cpp)

target_link_libraries (gemm_ambient.out ambient ${BLAS_LAPACK_serial} ${MPI_CXX_LIBRARIES} ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY}) 
target_link_libraries (gemm_blas.out ambient ${BLAS_LAPACK_multi} ${MPI_CXX_LIBRARIES} ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY}) 
target_link_libraries (pdgemm_scalapack.out ${BLAS_SCALAPACK_serial} ${MPI_CXX_LIBRARIES} ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY}) 

add_test(gemm_blas.out gemm_blas.out)
add_test(gemm_ambient.out gemm_ambient.out)
add_test(pdgemm_scalapack.out pdgemm_scalapack.out)

if(ENABLE_PLOT)
    add_executable (plot.out plot.cpp)
    target_link_libraries (plot.out ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY}) 
    add_test(plot.out plot.out)
endif(ENABLE_PLOT)
