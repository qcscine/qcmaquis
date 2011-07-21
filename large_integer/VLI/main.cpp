#include <iostream>
#include <cstdio>

#include "gmpxx.h"

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"

#include "monome/vector_polynomial_cpu.h"
#include "monome/vector_polynomial_gpu.h"
#include "monome/polynome_gpu.h"
#include "monome/polynome_cpu.h"
#include "monome/monome.h"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "utils/timings.h"
#define SIZE_POLY 20
#define SIZE_VECTOR 16

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
    
    vli_cpu<TYPE,8> A;
    vli_cpu<TYPE,8> B;
    
    A[0]= 10;
      
    B[0]= 10;
    
	vli::vli_cpu<TYPE,8> C(0);
	
    vli::vli_gpu<TYPE,8> D(A);
	vli::vli_gpu<TYPE,8> E(B);
	vli::vli_gpu<TYPE,8> F(0);

    
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
    mpz_class Cgmp(C.get_str());

    Timer TimerA("VLI");
    TimerA.begin();
    for(int i= 0 ; i< 10 ;i++)
     	A*=B;
	TimerA.end();    

    Timer TimerB("GMP");
    TimerB.begin();
    for(int i= 0 ; i< 10 ;i++)
        Agmp*=Bgmp;
    TimerB.end(); 
    
    std::cout << Agmp << std::endl;
    std::cout << A << std::endl;

	GPU->instance().destructor();
    

    return 0;
}
