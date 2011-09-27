   add_definitions(-D_MIX -DUSE_VLI_INTEGERS_CPU -fopenmp)
   add_executable (a_vli_mix.out ${PROJECT_SOURCE_DIR}/minimal_high_t_expansion.cpp)
   target_link_libraries (a_vli_mix.out -fopenmp ${VLI_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES})
