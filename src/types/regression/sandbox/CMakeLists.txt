compile_with_ambient()
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${AMBIENT_SOURCES_DIR}")

add_definitions(-DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

if(BLAS_mkl_core_LIBRARY) # intel
    include_directories("${MPI_CXX_INCLUDE_PATH}") # cluster path
    set(BLAS_LAPACK ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY} ${PLASMA_LIBRARY}) 
endif()

if(BLAS_Accelerate_LIBRARY) # apple
    include_directories("${MPI_INCLUDE_PATH}") # cluster path
    set(MPI_CXX_LIBRARIES ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    add_definitions(-DAPPLE -framework Accelerate)
    set(BLAS_LAPACK ${BLAS_vecLib_LIBRARY}) 
endif()

#add_executable (test.binary test.source.cpp)

foreach( sandbox_test
         #test.binary
)
    target_link_libraries (${sandbox_test} ambient ${BLAS_LAPACK} ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES}) 
    add_test(${sandbox_test} ${sandbox_test})
endforeach(sandbox_test)
