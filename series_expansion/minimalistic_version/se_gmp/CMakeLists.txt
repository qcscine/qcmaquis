   add_definitions(-DUSE_GMP_INTEGERS -fopenmp)
   add_executable (a_gmp.out ${PROJECT_SOURCE_DIR}/minimal_high_t_expansion.cpp)
   target_link_libraries (a_gmp.out -fopenmp ${GMP_LIBRARY_CXX} ${GMP_LIBRARY_C})
