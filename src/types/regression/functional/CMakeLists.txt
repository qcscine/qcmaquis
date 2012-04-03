compile_with_ambient()
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${AMBIENT_SOURCES_DIR}")

add_definitions (-DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
add_definitions(-DMPI_PARALLEL)

if(BLAS_mkl_core_LIBRARY) #intel
    include_directories("${MPI_CXX_INCLUDE_PATH}") #cluster path
    set(BLAS_LAPACK ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY}) 
endif()

if(BLAS_Accelerate_LIBRARY) #apple
    include_directories("${MPI_INCLUDE_PATH}") #cluster path
    set(MPI_CXX_LIBRARIES ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    add_definitions(-DAPPLE -DMPI_PARALLEL -framework Accelerate)
    set(BLAS_LAPACK ${BLAS_vecLib_LIBRARY}) 
endif()

add_executable (blas3_multiplication.out blas3_multiplication.cpp)
add_executable (blas3_multiplication_pdiag.out blas3_multiplication_pdiag.cpp)
add_executable (blas3_addition.out blas3_addition.cpp)
add_executable (blas3_substraction.out blas3_substraction.cpp)
add_executable (cast.out cast.cpp)
add_executable (transpose.out transpose.cpp)
add_executable (remove_rows.out remove_rows.cpp)
add_executable (remove_cols.out remove_cols.cpp)
add_executable (trace.out trace.cpp)
add_executable (resize.out resize.cpp)
add_executable (svd_lapack.out svd_lapack.cpp)
add_executable (heev_lapack.out heev_lapack.cpp)
add_executable (syev_lapack.out syev_lapack.cpp)
add_executable (sqrt.out sqrt.cpp)
add_executable (exp.out exp.cpp)
add_executable (serial_access.out serial_access.cpp)
#add_executable (identity.out identity.cpp)

foreach(functional_test 
       blas3_multiplication.out 
       blas3_multiplication_pdiag.out 
       blas3_addition.out 
       blas3_substraction.out 
       cast.out 
       transpose.out 
       remove_rows.out 
       remove_cols.out 
       trace.out 
       resize.out 
       svd_lapack.out 
       heev_lapack.out 
       syev_lapack.out
       sqrt.out
       exp.out
       serial_access.out
#      identiry.out
       )
       target_link_libraries (${functional_test} ambient ${BLAS_LAPACK} ${MPI_CXX_LIBRARIES} ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY}) 
       add_test(${functional_test} ${functional_test})
endforeach(functional_test)
