
//#include "vli_config.h"
#include <iostream>

#define SIZE_BITS 256


#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "vli_cpu/vli_vector_cpu.hpp"
#include "vli_gpu/vli_vector_gpu.hpp"

#include "utils/timings.h"

typedef int TYPE; 
using vli::vli_vector_cpu;
using vli::vli_vector_gpu;
using vli::vli_cpu;

int main (int argc, char * const argv[]) 
{
    int SIZE = 4;
 	gpu::gpu_manager* GPU;
	GPU->instance();
	
    vli_vector_cpu<int> a(SIZE);
    vli_vector_cpu<int> b(a);
    vli_vector_cpu<int> c(SIZE);
	
    
	for (int i =0; i < SIZE; i++)
	{
        a.init(i,3,255,255,255);
        b.init(i,3,255,255,255);
	}
    
    c = a+b;
    
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << c << std::endl;

    a*=b;
	
	vli_vector_gpu<int> a_gpu(a);
	vli_vector_gpu<int> b_gpu(b);
	vli_vector_gpu<int> c_gpu(SIZE);

    c_gpu = a_gpu + b_gpu;
    
    std::cout << a_gpu << std::endl;
    std::cout << b_gpu << std::endl;
    std::cout << c_gpu << std::endl;
    
   

/*
	vli_vector_gpu<vli_cpu<int> > b_gpu(b);
    
    //	c = a+b;
    
    //	c_gpu(c);//=a_gpu+b_gpu;
	
	std::cout << a << std::endl;	
	std::cout << a_gpu << std::endl;	
    
    //	BOOST_CHECK_EQUAL(a[0],b[0]);
 */
	GPU->instance().destructor();
}
