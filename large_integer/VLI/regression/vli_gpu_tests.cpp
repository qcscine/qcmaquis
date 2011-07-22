#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "gmpxx.h"

using vli::vli_cpu;
using vli::max_int_value;
using vli::vli_gpu;

typedef boost::mpl::list<unsigned int, unsigned long int> test_types;

BOOST_AUTO_TEST_CASE(gpu_manager)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    //TODO correct type?
	double FreqGPU = GPU->instance().GetDeviceProperties().clockRate;
    std::cout<<" FreqGPU "<<FreqGPU<<" Hz"<<std::endl;
    GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(constructors_test, T, test_types)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    vli_gpu<T,8> a;
    vli_gpu<T,8> b(0);

    BOOST_CHECK_EQUAL(a,b);
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(copy_construction, T, test_types)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
    vli_cpu<T,8> A;

    A[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
	
    vli::vli_gpu<T,8> B(A);
	    
	BOOST_CHECK_EQUAL(A,B);
	
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(proxy_access, T, test_types)
{

    // Read access
    vli_cpu<T,8> c;
    c[0] = 219837;
    c[1] = 7823;
    c[2] = 937;
    c[3] = 0;
    c[4] = 898123;
    c[5] = 716324;
    c[6] = 123;
    c[7] = 52782;

    vli_gpu<T,8> cg(c);

    for(std::size_t i=0; i < 8; ++i)
        BOOST_CHECK_EQUAL(c[i],cg[i]);
    
    // Write access
    c[0] = 845234;
    c[1] = 4562;
    c[2] = 98972;
    c[3] = 2343;
    c[4] = 10;
    c[5] = 7723;
    c[6] = 173643;
    c[7] = 52873;

    cg[0] = 845234;
    cg[1] = 4562;
    cg[2] = 98972;
    cg[3] = 2343;
    cg[4] = 10;
    cg[5] = 7723;
    cg[6] = 173643;
    cg[7] = 52873;
    
    for(std::size_t i=0; i < 8; ++i)
        BOOST_CHECK_EQUAL(c[i],cg[i]);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(serial_addition, T, test_types)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
		
    vli_cpu<T,8> A;
    vli_cpu<T,8> B;
    
    A[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
        
    B[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
	vli::vli_cpu<T,8> C(0);
	
    vli::vli_gpu<T,8> D(A);
	vli::vli_gpu<T,8> E(B);
	vli::vli_gpu<T,8> F(0);

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

BOOST_AUTO_TEST_CASE_TEMPLATE(serial_substraction, T, test_types)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<T,8> A;
    vli_cpu<T,8> B;
    
    A[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
    B[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
	vli::vli_cpu<T,8> C(0);
	
    vli::vli_gpu<T,8> D(A);
	vli::vli_gpu<T,8> E(B);
	vli::vli_gpu<T,8> F(0);
    
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

BOOST_AUTO_TEST_CASE_TEMPLATE(serial_multiplication, T, test_types)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<T,8> A;
    vli_cpu<T,8> B;
    
    A[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
    B[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
	vli::vli_cpu<T,8> C(0);
	
    vli::vli_gpu<T,8> D(A);
	vli::vli_gpu<T,8> E(B);
	vli::vli_gpu<T,8> F(0);
    
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

BOOST_AUTO_TEST_CASE_TEMPLATE(PlusAssign, T, test_types)
{
    gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<T,8> A;
    vli_cpu<T,8> B;
    
    A[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
    B[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
	
    vli::vli_gpu<T,8> D(A);
	vli::vli_gpu<T,8> E(B);
	
	A+=B;
	D+=E;
	
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
       
    Agmp+=Bgmp;
    
	BOOST_CHECK_EQUAL(A,D);
	BOOST_CHECK_EQUAL(A.get_str(),Bgmp.get_str());
    
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(MinusAssign, T, test_types)
{
    gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<T,8> A;
    vli_cpu<T,8> B;
    
    A[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
    B[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
	
    vli::vli_gpu<T,8> D(A);
	vli::vli_gpu<T,8> E(B);
	
	A-=B;
	D-=E;
	
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
    
    Agmp-=Bgmp;
    
	BOOST_CHECK_EQUAL(A,D);
	BOOST_CHECK_EQUAL(A.get_str(),Bgmp.get_str());
    
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(MultiplyAssign, T, test_types)
{
    gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<T,8> A;
    vli_cpu<T,8> B;
    
    A[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    A[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
    B[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    B[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,8> >::value);
    
	
    vli::vli_gpu<T,8> D(A);
	vli::vli_gpu<T,8> E(B);
	
	A*=B;
	D*=E;
	
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
    
    Agmp*=Bgmp;
    
	BOOST_CHECK_EQUAL(A,D);
	BOOST_CHECK_EQUAL(A.get_str(),Bgmp.get_str());
    
	GPU->instance().destructor();
}

