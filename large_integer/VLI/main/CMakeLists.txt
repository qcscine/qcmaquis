# Make sure the linker can find the vli library once it is built. 
link_directories (${PROJECT_BINARY_DIR}/src/vli) 

if(VLI_COMPILE_GPU)
  #  include_directories(${CUDA_INCLUDE_DIRS})
    add_definitions(-DVLI_USE_GPU -DVLI_SPLIT_PARAM=${VLI_SPLIT_PARAM} )
endif(VLI_COMPILE_GPU)

add_executable (a.out main.cpp )
target_link_libraries (a.out  vli)
if(VLI_COMPILE_GPU)
    target_link_libraries (a.out  accelerator_gpu ${CUDA_CUDART_LIBRARY} ${CUDA_LIBRARIES} ${OpenMP_CXX_FLAGS})
else(VLI_COMPILE_GPU)
    target_link_libraries (a.out  rt ${OpenMP_CXX_FLAGS})
endif(VLI_COMPILE_GPU)
