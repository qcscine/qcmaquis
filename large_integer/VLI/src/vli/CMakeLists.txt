
########################################################################
#
# Configure Header File Template
#
########################################################################
configure_file (
        "${PROJECT_SOURCE_DIR}/src/vli/vli_config.h.in"
        "${PROJECT_BINARY_DIR}/vli/vli_config.h"
        )


########################################################################
#
# Configure Header File Template
#
########################################################################
add_library(vli_cpu 
            detail/kernels_cpu_add_asm.cpp
            detail/kernels_cpu_sub_asm.cpp
            detail/kernels_cpu_mul_asm.cpp
            detail/kernels_cpu_muladd_asm.cpp
      )
     


########################################################################
#
# GPU Kernel Compile
#
########################################################################

if(VLI_COMPILE_GPU)
    CUDA_ADD_LIBRARY(vli_gpu detail/single_coefficient_task.cpp detail/gpu_hardware_carryover_implementation.cu)
endif(VLI_COMPILE_GPU)
