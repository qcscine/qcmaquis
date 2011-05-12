#define BOOST_TEST_MODULE vli_cpu
#include <iostream>
#include <boost/test/unit_test.hpp>


#define SIZE_BITS 256
#include "GpuManager.h"
#include "GpuManager.hpp"

#include "vli_number_cpu.hpp"
#include "vli_vector_cpu.hpp"
#include "vli_vector_gpu.hpp"

using vli::vli_cpu;
using vli::vli_vector;
using vli::vli_vector_gpu;

BOOST_AUTO_TEST_CASE( constructors_test )
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
    vli_vector_gpu<int> a(10);
    vli_vector_gpu<int> b(a);
    vli_vector_gpu<int> c(10);
	
	c=a+b;
	
	BOOST_CHECK_EQUAL(a,b);

	GPU->instance().destructor();
}


BOOST_AUTO_TEST_CASE( addition_vector )
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
    vli_vector<vli_cpu<int> > a(10);
    vli_vector<vli_cpu<int> > b(a);
    vli_vector<vli_cpu<int> > c(10);
	
	
	for (int i =0; i < 10; i++)
	{
	    a[i][0] = 255;
	    a[i][1] = 255;
	    a[i][2] = 255;
		
		b[i][0] = 255;
	    b[i][1] = 255;
	    b[i][2] = 255;
	}
	
	vli_vector_gpu<int> a_gpu(a);
	vli_vector_gpu<int> b_gpu(b);
	vli_vector_gpu<int> c_gpu(c);
	
	
	c = a+b;
	c_gpu=a_gpu+b_gpu;
	
	std::cout << c << std::endl;	
	std::cout << c_gpu << std::endl;	

	BOOST_CHECK_EQUAL(a,b);
	GPU->instance().destructor();
}
