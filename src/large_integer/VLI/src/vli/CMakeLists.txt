
########################################################################
#
# Configure Header File Template
#
########################################################################
configure_file (
        "${PROJECT_SOURCE_DIR}/src/vli/config.hpp.in"
        "${PROJECT_BINARY_DIR}/vli/config.hpp"
        )


########################################################################
#
# Configure Header File Template
#
########################################################################

add_library(vli SHARED 
            detail/cpu/kernels_cpu_add_asm.cpp
            detail/cpu/kernels_cpu_sub_asm.cpp
            detail/cpu/kernels_cpu_mul_asm.cpp
            detail/cpu/kernels_cpu_muladd_asm.cpp
      )

########################################################################
#
# GPU Kernel Compile
#
########################################################################

if(VLI_COMPILE_GPU)
    CUDA_ADD_LIBRARY(vli_gpu detail/gpu/tasklist/single_coefficient_task.cpp detail/gpu/polynomial_operations.cu)
endif(VLI_COMPILE_GPU)
