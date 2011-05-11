
//#include "vli_config.h"
#include <iostream>

#define SIZE_BITS 256

#include <cuda.h>
#include <cublas.h>

#include "GpuManager.h"
#include "GpuManager.hpp"
#include "vli_number_cpu.hpp"
#include "vli_number_gpu.hpp"
#include "vli_vector_cpu.hpp"
#include "vli_vector_gpu.hpp"

#define BOOST_TEST_MODULE gpu

//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp> //faster to compile

#include "timings.h"

typedef int TYPE; 

BOOST_AUTO_TEST_CASE(gpu_manager)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	TYPE FreqGPU = GPU->instance().GetDeviceProperties().clockRate;
    printf(" FreqGPU %d, Hz", FreqGPU);
    GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE(copy_construction)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
	vli::vli_cpu<TYPE> A(255);

	A[1] = 255;
	A[2] = 255;
	A[3] = 255;
	
    vli::vli_gpu<TYPE> B(A);
		
	BOOST_CHECK_EQUAL(A,B);
	
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE(serial_addition)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
	vli::vli_cpu<TYPE> A(255);
	vli::vli_cpu<TYPE> B(255);
	vli::vli_cpu<TYPE> C(0);

	A[1] = 255;
	A[2] = 255;
	A[3] = 255;
	/** this number is the result by hand */
	vli::vli_cpu<TYPE> Res;

	Res[0] = 254;
	Res[1] = 0;
	Res[2] = 0;
	Res[3] = 0;
	Res[4] = 1;
	/** end calculation by hand */
	
    vli::vli_gpu<TYPE> D(A);
	vli::vli_gpu<TYPE> E(B);
	vli::vli_gpu<TYPE> F(0);

	C = A+B;
	F = D+E;
	
	BOOST_CHECK_EQUAL(C,F);
	BOOST_CHECK_EQUAL(C,Res);
	BOOST_CHECK_EQUAL(Res,F);

	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE(serial_multiplication)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
	vli::vli_cpu<TYPE> A(255);
	vli::vli_cpu<TYPE> B(255);
	vli::vli_cpu<TYPE> C(0);
	
	A[1] = 255;
	A[2] = 255;
	A[3] = 255;
	/** this number is the result by hand */
	vli::vli_cpu<TYPE> Res;
	
	Res[0] = 254;
	Res[1] = 0;
	Res[2] = 0;
	Res[3] = 0;
	Res[4] = 1;
	
	
	
    vli::vli_gpu<TYPE> D(A);
	vli::vli_gpu<TYPE> E(B);
	vli::vli_gpu<TYPE> F(0);
	
	C = A*B;
	F = D*E;
	
	BOOST_CHECK_EQUAL(C,F);
//	BOOST_CHECK_EQUAL(C,Res);
//	BOOST_CHECK_EQUAL(Res,F);
	
	GPU->instance().destructor();

}


/*
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
	
	

	A = E*B;

	C = A;
	std::cout << A << std::endl;	

//	std::cout << B << std::endl;

	vli::vli_vector<vli::vli_cpu<TYPE> > U(1000);
	vli::vli_vector<vli::vli_cpu<TYPE> > V(1000);
	vli::vli_vector<vli::vli_cpu<TYPE> > W(1000);

    U[4] = A;
    U[5] = B;
    W[0] = vli::vli_cpu<int>(1);

    vli::vli_vector_gpu<TYPE> X(U);
    vli::vli_vector_gpu<TYPE> Y(V);
    vli::vli_vector_gpu<TYPE> Z(W);

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

*/