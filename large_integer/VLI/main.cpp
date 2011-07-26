#include <iostream>
#include <cstdio>

#include "gmpxx.h"

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"

#include "polynomial/vector_polynomial_cpu.hpp"
#include "polynomial/vector_polynomial_gpu.hpp"
#include "polynomial/polynomial_gpu.hpp"
#include "polynomial/polynomial_cpu.hpp"
#include "polynomial/monomial.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "utils/timings.h"
#define SIZE_POLY 20
#define SIZE_VECTOR 16

#define TYPE unsigned long int 

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
    A[0]= 100;
    B[0]= 100;
       
    mpz_class Agmp(A.get_str());
    mpz_class Bgmp(B.get_str());
 
    Timer TimerA("VLI");
    TimerA.begin();
    for(int i= 0 ; i< 250;++i)
     	A*=3;
     	A*=3;
     	A*=3;
     	A*=3;
     	A*=3;
     	A*=3;
     	A*=3;
     	A*=3;
     	A*=3;
     	A*=3;
    TimerA.end();    

    Timer TimerB("GMP");
    TimerB.begin();
    for(int i= 0 ; i< 250;++i)
        Agmp*=3;
        Agmp*=3;
        Agmp*=3;
        Agmp*=3;
        Agmp*=3;
        Agmp*=3;
        Agmp*=3;
        Agmp*=3;
        Agmp*=3;
        Agmp*=3;
    TimerB.end(); 

    std::cout << A << std::endl;
    std::cout << Agmp << std::endl;

	GPU->instance().destructor();
    

    return 0;
}
