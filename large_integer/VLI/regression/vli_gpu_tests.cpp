#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>


#define SIZE_BITS 256
#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "gmpxx.h"


#define TYPE unsigned int
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

BOOST_AUTO_TEST_CASE(constructors_test)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    vli_gpu<int,8> a;
    vli_gpu<int,8> b(0);

    BOOST_CHECK_EQUAL(a,b);
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE(copy_construction)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
    vli_cpu<TYPE,8> A;

    A[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
	
    vli::vli_gpu<TYPE,8> B(A);
	    
	BOOST_CHECK_EQUAL(A,B);
	
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE(serial_addition)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
		
    vli_cpu<TYPE,8> A;
    vli_cpu<TYPE,8> B;
    
    A[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
        
    B[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
	vli::vli_cpu<TYPE,8> C(0);
	
    vli::vli_gpu<TYPE,8> D(A);
	vli::vli_gpu<TYPE,8> E(B);
	vli::vli_gpu<TYPE,8> F(0);

	C = A+B;
	F = D+E;
	
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
    mpz_class Cgmp(C.get_str());
    
    Cgmp = Agmp+Bgmp;
    
	BOOST_CHECK_EQUAL(C,F);
	BOOST_CHECK_EQUAL(C.get_str(),Cgmp.get_str());
	BOOST_CHECK_EQUAL(F.get_str(),Cgmp.get_str());

	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE(serial_substraction)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<TYPE,8> A;
    vli_cpu<TYPE,8> B;
    
    A[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    B[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
	vli::vli_cpu<TYPE,8> C(0);
	
    vli::vli_gpu<TYPE,8> D(A);
	vli::vli_gpu<TYPE,8> E(B);
	vli::vli_gpu<TYPE,8> F(0);
    
	C = A-B;
	F = D-E;
	
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
    mpz_class Cgmp(C.get_str());
    
    Cgmp = Agmp-Bgmp;
    
	BOOST_CHECK_EQUAL(C,F);
	BOOST_CHECK_EQUAL(C.get_str(),Cgmp.get_str());
	BOOST_CHECK_EQUAL(F.get_str(),Cgmp.get_str());
    
	GPU->instance().destructor();
}





BOOST_AUTO_TEST_CASE(serial_multiplication)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<TYPE,8> A;
    vli_cpu<TYPE,8> B;
    
    A[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    B[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
	vli::vli_cpu<TYPE,8> C(0);
	
    vli::vli_gpu<TYPE,8> D(A);
	vli::vli_gpu<TYPE,8> E(B);
	vli::vli_gpu<TYPE,8> F(0);
    
	C = A*B;
	F = D*E;
	
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
    mpz_class Cgmp(C.get_str());
    
    Cgmp = Agmp*Bgmp;
    
	BOOST_CHECK_EQUAL(C,F);
	BOOST_CHECK_EQUAL(C.get_str(),Cgmp.get_str());
	BOOST_CHECK_EQUAL(F.get_str(),Cgmp.get_str());
    
	GPU->instance().destructor();

}

BOOST_AUTO_TEST_CASE(PlusAssign)
{
    gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<TYPE,8> A;
    vli_cpu<TYPE,8> B;
    
    A[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    B[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
	
    vli::vli_gpu<TYPE,8> D(A);
	vli::vli_gpu<TYPE,8> E(B);
	
	A+=B;
	D+=E;
	
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
       
    Agmp+=Bgmp;
    
	BOOST_CHECK_EQUAL(A,D);
	BOOST_CHECK_EQUAL(A.get_str(),Bgmp.get_str());
    
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE(MinusAssign)
{
    gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<TYPE,8> A;
    vli_cpu<TYPE,8> B;
    
    A[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    B[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
	
    vli::vli_gpu<TYPE,8> D(A);
	vli::vli_gpu<TYPE,8> E(B);
	
	A-=B;
	D-=E;
	
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
    
    Agmp-=Bgmp;
    
	BOOST_CHECK_EQUAL(A,D);
	BOOST_CHECK_EQUAL(A.get_str(),Bgmp.get_str());
    
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE(MultiplyAssign)
{
    gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<TYPE,8> A;
    vli_cpu<TYPE,8> B;
    
    A[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    A[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    B[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    B[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
	
    vli::vli_gpu<TYPE,8> D(A);
	vli::vli_gpu<TYPE,8> E(B);
	
	A*=B;
	D*=E;
	
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
    
    Agmp*=Bgmp;
    
	BOOST_CHECK_EQUAL(A,D);
	BOOST_CHECK_EQUAL(A.get_str(),Bgmp.get_str());
    
	GPU->instance().destructor();
}

