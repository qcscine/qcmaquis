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

	vli::vli_matrix< int > A(25,5);
	vli::vli_matrix< int > B(25,5);
	vli::vli_matrix< int > C(25,5);

	vli::vli_matrix< int > A_c(25,5);
	vli::vli_matrix< int > B_c(25,5);
	vli::vli_matrix< int > C_c(25,5);
	
	

	for(int i=0 ; i < 25; i++)
	{
	A(0,i) = -1;
	A(1,i) = 2;
	A(2,i) = 0;
	A(3,i) = 1;
	A(4,i) = -3;
	
	B(0,i) = -2;
	B(1,i) = 1;
	B(2,i) = -1;
	B(3,i) = 0;
	B(4,i) = 2;
		
		
		A_c(0,i) = 2;
		A_c(1,i) = 2;
		A_c(2,i) = 3;
		A_c(3,i) = 3;
		A_c(4,i) = 0;
		
		B_c(0,i) = 1;
		B_c(1,i) = 0;
		B_c(2,i) = 1;
		B_c(3,i) = 0;
		B_c(4,i) = 0;	
		

	}
	

	
	
	std::cout << " cpu avizienis " << std::endl;
	
	vli::addition_Avizienis_kernel_cpu(A,B,C);
	std::cout << C << std::endl;
	
	std::cout << " gpu avizienis " << std::endl;
	gpu::matrix_matrix_addition(A_c,B_c,C_c);
	
	
	std::cout << " cpu classique " << std::endl;
	vli::addition_classic_kernel_cpu(A_c,B_c,C_c);
	std::cout << C_c << std::endl;
		
	GPU->instance().destructor();
	
    return 0;
}
