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

#define SIZE 8
#define TYPE unsigned int 

using vli::vli_cpu;
using vli::vli_gpu;
using vli::monomial;
using vli::polynomial_cpu;
using vli::polynomial_gpu;
using vli::vector_polynomial_gpu;
using vli::vector_polynomial_cpu;



int main (int argc, char * const argv[]) 
{

    gpu::gpu_manager* GPU;
    GPU->instance();
    
    vli_cpu<TYPE,SIZE> a;
    vli_cpu<TYPE,SIZE> b;
    vli_cpu<TYPE,SIZE> c;
    



    a[0]=100;
    a[1]=99;
    a[2]=12;
    a[2]=136;
    
    b[0]=234;
    b[1]=20;
    b[2]=100;
    b[3]=233;
    b[5]=1;
    b[6]=0;
    b[7]=0;

    vli_gpu<TYPE,SIZE> agpu(a);
    vli_gpu<TYPE,SIZE> bgpu(b);
    vli_gpu<TYPE,SIZE> cgpu(c);
    
    c -= a;
    cgpu -= agpu;
    
    
    std::cout << c << std::endl;
    std::cout << cgpu << std::endl;
    
    c = b* a;
       
    mpz_t agmp, bgmp;                 	


    mpz_init_set_str (agmp, a.get_char(), 10);
    mpz_init_set_str (bgmp, b.get_char(), 10);
    mpz_mul (bgmp, bgmp, agmp);	

    std::cout << "    Dirty a    "<< a.get_str() << std::endl;
    std::cout << "    Dirty b    "<< b.get_str() << std::endl; 
    std::cout << "    Dirty c    "<< c.get_str() << std::endl;


    gmp_printf ("%s is an mpz %Zd\n", "here", bgmp);



	GPU->instance().destructor();
    return 0;

}
