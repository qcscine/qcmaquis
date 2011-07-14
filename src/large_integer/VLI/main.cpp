#include <iostream>
#include <cstdio>

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"

#include "monome/vector_polynomial_cpu.h"
#include "monome/vector_polynomial_gpu.h"
#include "monome/polynome_gpu.h"
#include "monome/polynome_cpu.h"
#include "monome/monome.h"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"

using vli::vli_cpu;
using vli::vli_gpu;
using vli::monomial;
using vli::polynomial_cpu;
using vli::polynomial_gpu;
using vli::vector_polynomial_gpu;
using vli::vector_polynomial_cpu;

#define SIZE 8

int main (int argc, char * const argv[]) 
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<int,SIZE> a;
    vli_cpu<int,SIZE> b;
        
    a[0] = 155;
    a[1] = 155;
    a[2] = 155;
    a[3] = 155;
 
    b[0] = 200;
    b[1] = 200;
    b[2] = 200;
    b[3] = 200;

    
    a-=b;
     
    std::cout << a << std::endl;   
//    std::cout << a.BaseTen() << std::endl;
 
    
	GPU->instance().destructor();
}
