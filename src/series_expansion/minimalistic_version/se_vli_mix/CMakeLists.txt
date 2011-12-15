   include_directories(${CUDA_INCLUDE_DIRS})
   add_definitions(-DUSE_VLI_INTEGERS_CPU -DVLI_USE_GPU ${OPENMP} -DVLI_SPLIT_PARAM=${SE_SPLIT_PARAM})
   add_executable (a_vli_mix.out ${PROJECT_SOURCE_DIR}/minimal_high_t_expansion.cpp)
   target_link_libraries (a_vli_mix.out ${OPENMP} ${VLI_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES})
