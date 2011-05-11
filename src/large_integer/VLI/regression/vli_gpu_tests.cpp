#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>


#define SIZE_BITS 256
#include "GpuManager.h"
#include "GpuManager.hpp"
#include "vli_number_cpu.hpp"
#include "vli_number_gpu.hpp"


#define TYPE int
using vli::vli_cpu;
using vli::vli_gpu;


BOOST_AUTO_TEST_CASE(gpu_manager)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	TYPE FreqGPU = GPU->instance().GetDeviceProperties().clockRate;
    std::cout<<" FreqGPU "<<FreqGPU<<" Hz"<<std::endl;
    GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE( constructors_test )
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    vli_gpu<int> a;
    vli_gpu<int> b(0);
    BOOST_CHECK_EQUAL(a,b);
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
	
	vli::vli_cpu<TYPE> A(254);
	vli::vli_cpu<TYPE> B(254);
	vli::vli_cpu<TYPE> C(0);
	
	A[1] = 255;
	A[2] = 255;
//	A[3] = 255;
	/** this number is the result by hand */
	vli::vli_cpu<TYPE> Res;
	
	Res[0] = 4;
	Res[1] = 254;
	Res[2] = 255;
	Res[3] = 253;
	
    vli::vli_gpu<TYPE> D(A);
	vli::vli_gpu<TYPE> E(B);
	vli::vli_gpu<TYPE> F(0);
	
	C = A*B;
	F = D*E;
	
	BOOST_CHECK_EQUAL(C,F);
	BOOST_CHECK_EQUAL(C,Res);
	BOOST_CHECK_EQUAL(Res,F);
	
	GPU->instance().destructor();

}

BOOST_AUTO_TEST_CASE(baseten_addition)
{
	vli::vli_cpu<TYPE> A(254);
	vli::vli_cpu<TYPE> B(254);
	vli::vli_cpu<TYPE> C(0);
	
	A[1] = 255;
	B[1] = 255;
	
	C=A+B;
	
	
	std::size_t ATen = A.BaseTen();
	std::size_t BTen = B.BaseTen();
	std::size_t CTen = C.BaseTen();
	std::size_t CTenRes = 0;

	CTenRes = ATen + BTen;

	BOOST_CHECK_EQUAL(CTen,CTenRes);	
}

BOOST_AUTO_TEST_CASE(baseten_multiplication)
{
	vli::vli_cpu<TYPE> A(254);
	vli::vli_cpu<TYPE> B(254);
	vli::vli_cpu<TYPE> C(0);
	
	A[1] = 255;
	B[1] = 255;
	
	C=A*B;
	
	std::size_t ATen = A.BaseTen();
	std::size_t BTen = B.BaseTen();
	std::size_t CTen = C.BaseTen();
	std::size_t CTenRes = 0;
	
	CTenRes = ATen * BTen;

	BOOST_CHECK_EQUAL(CTen,CTenRes);
}
