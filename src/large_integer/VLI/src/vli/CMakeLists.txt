
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
# GPU Kernel Compile
#
########################################################################

if(VLI_COMPILE_GPU)
    CUDA_ADD_LIBRARY(vli_gpu detail/kernels_gpu.cu)
endif(VLI_COMPILE_GPU)
