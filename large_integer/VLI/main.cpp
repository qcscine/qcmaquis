#include <iostream>

#define SIZE_BITS 256

#include <cuda.h>
#include <cublas.h>

#include "definition.h"

#include "GpuManager.h"
#include "GpuManager.hpp"
#include "vli_number_gpu.h"
#include "vli_number_cpu.h"
#include "vli_number_cpu.hpp"
#include "vli_number_gpu.hpp"
#include "vli_vector_cpu.hpp"


#include "timings.h"

typedef int TYPE; 

int main (int argc, char * const argv[]) 
{

	
	gpu::gpu_manager* GPU;
	GPU->instance();
	TYPE FreqGPU = GPU->instance().GetDeviceProperties().clockRate;
	std::cout << " FreqGPU : " << FreqGPU << std::endl;

	vli::vli_cpu<TYPE> A(254);
	vli::vli_cpu<TYPE> B(254);	
	vli::vli_cpu<TYPE> C;

	vli::vli_gpu<TYPE> D(A);
	std::cout << D << std::endl;
	vli::vli_gpu<TYPE> E(B);
	std::cout << E << std::endl;
	vli::vli_gpu<TYPE> F(0);

	
	A*=B;	
	multiply_gpu(D.p(),E.p(),F.p(),1,D.size());
	
//	D*=E;
	
	std::cout << " cpu "<< A << std::endl;
	std::cout << " gpu "<< F << std::endl;	
/*	
	A = E*B;

	C = A;
	std::cout << A << std::endl;	

//	std::cout << B << std::endl;
	
	vli::vli_vector<TYPE> U(4);
	vli::vli_vector<TYPE> V(4);
	vli::vli_vector<TYPE> W(4);
	
	U = V + W;
	U = V * W;
*/
	GPU->instance().destructor();
	
    return 0;
}
