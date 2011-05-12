#define BOOST_TEST_MODULE vli_cpu
#include <iostream>
#include <boost/test/unit_test.hpp>


#define SIZE_BITS 256
#include "GpuManager.h"
#include "GpuManager.hpp"

#include "vli_number_cpu.hpp"
#include "vli_vector_cpu.hpp"
#include "vli_vector_gpu.hpp"

using vli::vli_vector_gpu;

BOOST_AUTO_TEST_CASE( constructors_test )
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
    vli_vector_gpu<int> a(10);
    vli_vector_gpu<int> b(a);
    vli_vector_gpu<int> c(10);
	
	c=a+b;

	std::cout << a << std::endl;
	std::cout << b << std::endl;
	std::cout << c << std::endl;	
//    BOOST_CHECK_EQUAL(a,b);
	GPU->instance().destructor();
}

/*
BOOST_AUTO_TEST_CASE( addition )
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
    vli_vector<vli_cpu<int> > a(10);
    vli_vector<vli_cpu<int> > b(a);
    vli_vector<vli_cpu<int> > c(10);
	
	
	for (int i =0; i < 10; i++)
	{
		
	}
	
	c=a+b;
	
	std::cout << a << std::endl;
	std::cout << b << std::endl;
	std::cout << c << std::endl;	
	//    BOOST_CHECK_EQUAL(a,b);
	GPU->instance().destructor();
}*/