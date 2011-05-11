
//#include "vli_config.h"
#include <iostream>

#define SIZE_BITS 256


#include "GpuManager.h"
#include "GpuManager.hpp"
#include "vli_number_cpu.hpp"
#include "vli_number_gpu.hpp"
#include "vli_vector_cpu.hpp"
#include "vli_vector_gpu.hpp"

#include "timings.h"

typedef int TYPE; 

int main (int argc, char * const argv[]) 
{

	
	gpu::gpu_manager* GPU;
	GPU->instance();
	TYPE FreqGPU = GPU->instance().GetDeviceProperties().clockRate;
	std::cout << " FreqGPU : " << FreqGPU << std::endl;

	vli::vli_cpu<TYPE> A(255);
	vli::vli_cpu<TYPE> B(255);	
	vli::vli_cpu<TYPE> C1;
	vli::vli_cpu<TYPE> C2;

	A[1] = 255;


	vli::vli_gpu<TYPE> D(A);
	std::cout << D << std::endl;
	vli::vli_gpu<TYPE> E(B);
	std::cout << E << std::endl;
	vli::vli_gpu<TYPE> F1(0);
	vli::vli_gpu<TYPE> F2(0);
	

	

	
	C1 = B+A;	
 //  vli::detail::multiply_gpu(D.p(),E.p(),F.p(),1,vli::vli_gpu<TYPE>::size);
	
//	D*=E;

    F1 = D+E;    
	
	std::cout << " cpu + "<< C1 << std::endl;
	std::cout << " gpu + "<< F1 << std::endl;	

	std::cout << " cpu * "<< C2 << std::endl;
	std::cout << " gpu * "<< F2 << std::endl;	
	
	

	C1 = A;
	std::cout << A << std::endl;	

//	std::cout << B << std::endl;

	vli::vli_vector<vli::vli_cpu<TYPE> > U(1000);
	vli::vli_vector<vli::vli_cpu<TYPE> > V(1000);
	vli::vli_vector<vli::vli_cpu<TYPE> > W(1000);

    U[4] = A;
    U[5] = B;
    W[0] = vli::vli_cpu<int>(1);

    U *= entrywise(W);
    W = entrywise_product(U,V);

    vli::vli_vector_gpu<TYPE> X(U);
    vli::vli_vector_gpu<TYPE> Y(V);
    vli::vli_vector_gpu<TYPE> Z(W);

    Z = entrywise_product(X,Y);
//    X[4] = A;
//    X[5] = B;


    Timer gpu_timer("GPU");
    Timer cpu_timer("CPU");

    gpu_timer.begin();
    for(unsigned int i=0; i< 100000; ++i)    
        Z = X + Z;
    gpu_timer.end();


    cpu_timer.begin();
    for(unsigned int i=0; i< 100000; ++i)    
        W = U + W;
    cpu_timer.end();


	GPU->instance().destructor();
	
    return 0;
}
