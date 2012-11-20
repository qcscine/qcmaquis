compile_with_ambient()
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${AMBIENT_SOURCES_DIR}")

add_definitions(-DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
add_definitions(-DAMBIENT)

if(BLAS_mkl_core_LIBRARY) # intel
    include_directories("${MPI_CXX_INCLUDE_PATH}") # cluster path
    set(BLAS_LAPACK ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLAS_mkl_sequential_LIBRARY} ${PLASMA_LIBRARY}) 
endif()

if(BLAS_Accelerate_LIBRARY) # apple
    include_directories("${MPI_INCLUDE_PATH}") # cluster path
    set(MPI_CXX_LIBRARIES ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    add_definitions(-DAPPLE -DAMBIENT -framework Accelerate)
    set(BLAS_LAPACK ${BLAS_vecLib_LIBRARY}) 
endif()

add_executable( general.transpose.out   general.transpose.cpp )
add_executable( general.rm_rows.out     general.rm_rows.cpp   )
add_executable( general.rm_cols.out     general.rm_cols.cpp   )
add_executable( general.resize.out      general.resize.cpp    )
add_executable( general.elements.out    general.elements.cpp  )
add_executable( general.identity.out    general.identity.cpp  )
add_executable( general.trace.out       general.trace.cpp     )
add_executable( blas.add.out            blas.add.cpp          )
add_executable( blas.sub.out            blas.sub.cpp          )
add_executable( blas.gemm.out           blas.gemm.cpp         )
add_executable( blas.gemmd.out          blas.gemmd.cpp        )
add_executable( lapack.qr.out           lapack.qr.cpp         )
add_executable( lapack.lq.out           lapack.lq.cpp         )
add_executable( lapack.svd.out          lapack.svd.cpp        )
add_executable( lapack.heev.out         lapack.heev.cpp       )
add_executable( aux.cast.out            aux.cast.cpp          )
add_executable( aux.sqrt.out            aux.sqrt.cpp          )
add_executable( aux.exp.out             aux.exp.cpp           )

foreach( functional_test
         blas.add.out
         blas.sub.out
         blas.gemm.out
         blas.gemmd.out
         general.resize.out
         general.rm_cols.out
         general.rm_rows.out
         general.transpose.out
         general.elements.out
         general.identity.out
         general.trace.out
         lapack.qr.out
         lapack.lq.out
         lapack.svd.out
         lapack.heev.out
         aux.cast.out
         aux.sqrt.out
         aux.exp.out
 )
    target_link_libraries (${functional_test} ambient ${BLAS_LAPACK} ${MPI_CXX_LIBRARIES} ${MAQUIS_TYPES_LIBRARIES}) 
    add_test(${functional_test} ${functional_test})
endforeach(functional_test)
