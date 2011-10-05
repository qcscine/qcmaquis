   include_directories(${CUDA_INCLUDE_DIRS})
   add_definitions(-DUSE_VLI_INTEGERS_CPU -DVLI_USE_GPU )
   add_executable (a_vli_gpu.out ${PROJECT_SOURCE_DIR}/minimal_high_t_expansion.cpp)
   target_link_libraries (a_vli_gpu.out ${VLI_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES})
