   
   include_directories(${CUDA_INCLUDE_DIRS})
   add_definitions(-D_OPEN -DUSE_VLI_INTEGERS_CPU -fopenmp)
   add_executable (a_vli_omp.out ${PROJECT_SOURCE_DIR}/minimal_high_t_expansion.cpp)
   target_link_libraries (a_vli_omp.out -fopenmp ${VLI_LIBRARY}  ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES}) 
