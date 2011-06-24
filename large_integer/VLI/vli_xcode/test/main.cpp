
//#include "vli_config.h"
#include <iostream>

#define SIZE_BITS 256


#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "vli_cpu/vli_vector_cpu.hpp"
#include "vli_gpu/vli_vector_gpu.hpp"
#include "monome/monome.h"

#include "utils/timings.h"

typedef int TYPE; 
using vli::vli_vector_cpu;
using vli::vli_vector_gpu;
using vli::vli_cpu;
using vli::monomial;
using vli::polynomial;

int main (int argc, char * const argv[]) 
{
    int SIZE = 2;
 	gpu::gpu_manager* GPU;
	GPU->instance();

    vli_cpu<int> a;
    vli_cpu<int> b(a);
    vli_cpu<int> c;
	
    
	for (int i =0; i < SIZE; i++)
	{
        a[i] = 255;
        b[i] = 255;
	}
    
    c = a*b;
    
    std::size_t ATen = a.BaseTen();
	std::size_t BTen = b.BaseTen();
	std::size_t CTen = c.BaseTen();
	std::size_t CTenRes = 0;
    
	CTenRes = ATen * BTen;
    
    std::cout << a << std::endl;    
    std::cout << c << std::endl;
    std::cout << CTen << std::endl;
    std::cout << CTenRes << std::endl;

    
/*	
	vli_vector_gpu<int> a_gpu(a);
	vli_vector_gpu<int> b_gpu(b);
	vli_vector_gpu<int> c_gpu(SIZE);

    c_gpu = a_gpu + b_gpu;
    
    std::cout << a_gpu << std::endl;
    std::cout << b_gpu << std::endl;
    std::cout << c_gpu << std::endl;

   */
    monomial<vli_cpu<int> > ma;
    *ma.coeff = a;
    polynomial<vli_vector_cpu<int> >  pb;
    polynomial<vli_vector_cpu<int> >  pc;
    
    std::cout <<  *ma.coeff << std::endl;  

	GPU->instance().destructor();
    return 0;
}
