#include <iostream>
#include <cstdio>

#include "gmp.h"

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
#define TYPE unsigned int 

int main (int argc, char * const argv[]) 
{

	gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<TYPE,SIZE> a;
    vli_cpu<TYPE,SIZE> b;
    vli_cpu<TYPE,SIZE> c;

    a[0]=100;
    a[1]=100;
    a[2]=1;

    b[0]=123;
    b[1]=245;
    b[2]=1;
    b[3]=1;
    b[5]=0;
    b[6]=0;
    b[7]=0;

    TYPE A = a.BaseTen();
    TYPE B = b.BaseTen();


    std::cout<<A<<std::endl;
    std::cout<<B<<std::endl;

    c = b* a;
    std::cout<<c<<std::endl;

    std::cout<<c.get_string()<<std::endl;
   
    
    mpz_t agmp, bgmp;                 	

    mpz_init_set_str (agmp, a.get_char(), 10);
    mpz_init_set_str (bgmp, b.get_char(), 10);
    mpz_mul (bgmp, bgmp, agmp);	

    gmp_printf ("%s is an mpz %Zd\n", "here", bgmp);
    

	GPU->instance().destructor();
    return 0;

}
