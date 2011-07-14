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

#define SIZE 4

int main (int argc, char * const argv[]) 
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<unsigned int,SIZE> a;
    vli_cpu<unsigned int,SIZE> b;
    a[0]=100;
    a[1]=100;
    a[2]=1;

    b[0]=123;
    b[1]=245;
    b[2]=2;
//    b[4]=245;
//    b[5]=245;
//    b[6]=231;
//    b[7]=12;

//    a.negate();

    int A = a.BaseTen();
    int B = b.BaseTen();


    std::cout<<A<<std::endl;
    std::cout<< A*B <<std::endl;

   std::cout<< a <<std::endl;
   std::cout<< a.get_string() <<std::endl;
   std::cout<< b <<std::endl;
    b*= a;
   std::cout<< b <<std::endl;
   std::cout<<b.BaseTen()<<std::endl;
    
	GPU->instance().destructor();
}
