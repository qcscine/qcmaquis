   add_definitions(-DUSE_VLI_INTEGERS_CPU)
   add_executable (a_vli.out ${PROJECT_SOURCE_DIR}/minimal_high_t_expansion.cpp)
   target_link_libraries (a_vli.out  ${VLI_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES})
