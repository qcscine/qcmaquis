compile_with_ambient()

add_definitions(-DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DAMBIENT_ALPS)

add_executable( blas.gemm   blas.gemm.cpp )
add_executable( lapack.qr   lapack.qr.cpp )
add_executable( lapack.lq   lapack.lq.cpp )

add_executable( alps.transpose   alps.transpose.cpp )
add_executable( alps.gemm        alps.gemm.cpp      )
add_executable( alps.gemmd       alps.gemmd.cpp     )
add_executable( alps.svd         alps.svd.cpp       )
add_executable( alps.heev        alps.heev.cpp      )
add_executable( alps.geev        alps.geev.cpp      )
add_executable( alps.inverse     alps.inverse.cpp   )
add_executable( alps.cast        alps.cast.cpp      )
add_executable( alps.exp         alps.exp.cpp       )
add_executable( alps.expd        alps.expd.cpp      )
add_executable( alps.conj        alps.conj.cpp      )
add_executable( alps.adjoint     alps.adjoint.cpp   )

foreach( third_party_benchmark
         blas.gemm
         lapack.qr
         lapack.lq
       )
  target_link_libraries (${third_party_benchmark} ${AMBIENT_THIRD_PARTY_LIBRARIES})
endforeach(third_party_benchmark)

foreach( third_party_test
         alps.gemm
         alps.gemmd
         alps.transpose
         alps.svd
         alps.heev
         alps.geev
         alps.inverse
         alps.cast
         alps.exp
         alps.expd
         alps.conj
         alps.adjoint
       )
  target_link_libraries (${third_party_test} ${AMBIENT_THIRD_PARTY_LIBRARIES})
  add_test(${third_party_test} ${third_party_test})
endforeach(third_party_test)
