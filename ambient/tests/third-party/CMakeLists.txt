########################################################################
#
# Dependencies
#
########################################################################

if(NOT ALPS_DIR AND DEFINED ENV{ALPS_ROOT})
  set(ALPS_DIR $ENV{ALPS_ROOT}/share/alps)
  set(PREVENT_ALPS_COMPILERS ON)
  find_package(ALPS NO_SYSTEM_ENVIRONMENT_PATH)
endif(NOT ALPS_DIR AND DEFINED ENV{ALPS_ROOT})

include(${PROJECT_SOURCE_DIR}/config/find_boost.cmake)
find_maquis_boost(COMPONENTS)
message(STATUS "Boost include: ${MAQUIS_Boost_INCLUDE_DIRS}")
message(STATUS "Boost libs: ${MAQUIS_Boost_LIBRARIES}")

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AMBIENT_CXX_FLAGS} -DAMBIENT_OMP ${OpenMP_CXX_FLAGS}")

########################################################################
#
# Compilation
#
########################################################################

include_directories(${MAQUIS_Boost_INCLUDE_DIRS})
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
  target_link_libraries(${third_party_benchmark} ${AMBIENT_LIBRARIES} -lboost_system)
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
  target_link_libraries(${third_party_test} ${AMBIENT_LIBRARIES} -lboost_system)
  add_test(${third_party_test} ${third_party_test})
endforeach(third_party_test)
