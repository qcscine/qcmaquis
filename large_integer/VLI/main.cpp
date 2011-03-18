#include <iostream>
#include <gmp.h>

#define SIZE_BITS 256

#include <cuda.h>
#include <cublas.h>

#include "GpuManager.h"
#include "vli_matrix.h"
#include "definition.h"

#include "vli_number_gpu.h"
#include "vli_number_cpu.h"

#include "matrix_gpu.h"
#include "matrix_gpu_functions.hpp"

#include "kernel_number.h"

#include "timings.h"
typedef int TYPE; 

#define MAXMAIN 1000000



int main (int argc, char * const argv[]) 
{

	
	gpu::gpu_manager* GPU;
	GPU->instance();
	TYPE FreqGPU = GPU->instance().GetDeviceProperties().clockRate;
	std::cout << " FreqGPU : " << FreqGPU << std::endl;

	vli::vli_cpu<TYPE> A(2);
	vli::vli_cpu<TYPE> B(2);
	vli::vli_cpu<TYPE> C;
	
	vli::addition_classic_cpu(A,B,C);
	
	std::cout << C << std::endl;
	
	
	
	
	GPU->instance().destructor();
	
    return 0;
}
