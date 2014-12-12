compile_with_ambient()
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${AMBIENT_SOURCES_DIR}")

add_definitions(-DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

add_executable( gemm.blas.out          gemm.blas.cpp          )
add_executable( gemm.ambient.out       gemm.ambient.cpp       )
add_executable( gemm.tiled.ambient.out gemm.tiled.ambient.cpp )
add_executable( lq.ambient.out         lq.ambient.cpp         )
add_executable( lq.lapack.out          lq.lapack.cpp          )
add_executable( qr.ambient.out         qr.ambient.cpp         )
add_executable( qr.lapack.out          qr.lapack.cpp          )

foreach( performance_test
         gemm.blas.out         
         qr.lapack.out         
         lq.lapack.out         
         gemm.ambient.out      
         gemm.tiled.ambient.out
         lq.ambient.out        
         qr.ambient.out        
       )
    target_link_libraries (${performance_test} ${AMBIENT_TESTS_LIBRARIES}) 
    add_test(${performance_test} ${performance_test})
endforeach(performance_test)
