#include <iostream>

#include <cuda.h>
#include <cublas.h>

#include "GpuManager.h"
#include "vli_matrix.h"
#include "definition.h"

#include "matrix_gpu.h"
#include "matrix_gpu_functions.hpp"
#include "kernels.h"

typedef int TYPE; 

int main (int argc, char * const argv[]) 
{

	
	gpu::gpu_manager* GPU;
	GPU->instance();
	int FreqGPU = GPU->instance().GetDeviceProperties().clockRate;
	std::cout << " FreqGPU : " << FreqGPU << std::endl;

	vli::vli_matrix< int > A_c(8,5);
	vli::vli_matrix< int > B_c(8,5);
	vli::vli_matrix< int > Caddition_cpu(8,5);
	vli::vli_matrix< int > Caddition_gpu(8,5);	
	vli::vli_matrix< int > Cmultiplication_cpu(8,5);
	vli::vli_matrix< int > Cmultiplication_gpu(8,5);
	

	for(int i=0 ; i < 8; i++)
	{
		A_c(0,i) = 255;
		A_c(1,i) = 255;
		A_c(2,i) = 0;
		A_c(3,i) = 0;
		A_c(4,i) = 0;
		
		B_c(0,i) = 255;
		B_c(1,i) = 255;
		B_c(2,i) = 0;
		B_c(3,i) = 0;
		B_c(4,i) = 0;	
	}
	

	
	std::cout << " cpu addition " << std::endl;
	vli::addition_classic_cpu(A_c,B_c,Caddition_cpu);
	std::cout << Caddition_cpu << std::endl;

	std::cout << " gpu addition " << std::endl;
 	gpu::addition_classic_gpu(A_c,B_c,Caddition_gpu);
	
	
	
	std::cout << " cpu multiplication " << std::endl;
	vli::multiplication_classic_kernel_cpu(A_c,B_c,Cmultiplication_cpu);
	std::cout << Cmultiplication_cpu << std::endl;
			
	std::cout << " gpu multiplication " << std::endl;
	gpu::multiplication_classic_gpu(A_c,B_c,Cmultiplication_gpu);
	
	GPU->instance().destructor();
	
    return 0;
}
