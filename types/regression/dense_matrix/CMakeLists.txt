include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

add_definitions (-DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

if(BLAS_mkl_core_LIBRARY) #intel
    if(BLAS_MULTITHREADED)
        set(BLAS_LAPACK ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_intel_thread_LIBRARY} iomp5 ${CMAKE_THREAD_LIBS_INIT}) 
    else()
        set(BLAS_LAPACK ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY}) 
    endif()
endif()

if(BLAS_Accelerate_LIBRARY) #apple
    add_definitions(-DAPPLE)
    set(BLAS_LAPACK ${BLAS_Accelerate_LIBRARY}) 
endif()

add_executable (dense_matrix.test dense_matrix_unit_tests.cpp)

target_link_libraries (dense_matrix.test ${BLAS_LAPACK} ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY}) 

add_test(dense_matrix_unit_tets dense_matrix.test)
