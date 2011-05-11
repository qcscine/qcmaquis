#define BOOST_TEST_MODULE vli_cpu
#include <iostream>
#include <boost/test/unit_test.hpp>


#define SIZE_BITS 256
#include "vli_vector_gpu.hpp"

using vli::vli_vector_gpu;

BOOST_AUTO_TEST_CASE( constructors_test )
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
    vli_vector_gpu<int> a(10);
    vli_vector_gpu<int> b(a);
	
    BOOST_CHECK_EQUAL(a,b);
	GPU->instance().destructor();

}


