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
#define SIZE_POLY 4
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
    
    polynomial_cpu<vli_cpu<TYPE,8>, SIZE_POLY> pa;
    
    for(int i=0; i<3; i++){
        for(int j = 0; j < SIZE_POLY; j++ ){
            for(int k= 0; k< SIZE_POLY; k++){
                pa(j,k)[i]= MAX_VALUE-1;                
            }
        }
    }
    
    polynomial_gpu<vli_gpu<TYPE,8>, SIZE_POLY> pagpu(pa);
    
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,SIZE_POLY> > VaGPU(SIZE_VECTOR);
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,SIZE_POLY> > VbGPU(SIZE_VECTOR);
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,SIZE_POLY> > VcGPU;
    
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,SIZE_POLY> > VaCPU(SIZE_VECTOR);
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,SIZE_POLY> > VbCPU(SIZE_VECTOR);
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,SIZE_POLY> > VcCPU;
    
    for(int i=0;i < SIZE_VECTOR;i++){
        VaCPU[i]=pa;
        VbCPU[i]=pa;
        
        VaGPU[i]=pa;
        VbGPU[i]=pa;
    }
    
    Timer A("CPU");
    A.begin();
    VcCPU =  inner_product(VaCPU,VbCPU);
    A.end();
    
    TimerCuda B("GPU");
    B.begin();
    VcGPU =  inner_product(VaGPU,VbGPU);
    B.end();
    
  //  std::cout << VcCPU << std::endl;
  //  std::cout << VcGPU << std::endl;
 
    if(VcGPU == VcCPU){ printf("ok \n");}else{printf("no ok \n");}
    
    GPU->instance().destructor();
    return 0;
}
