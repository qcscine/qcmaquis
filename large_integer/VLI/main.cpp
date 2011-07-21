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
#define SIZE_POLY 10
#define SIZE_VECTOR 256

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
/*
    polynomial_cpu<vli_cpu<TYPE,8>, SIZE_POLY> pb;
    polynomial_cpu<vli_cpu<TYPE,8>, SIZE_POLY> pc;
  */  
    for(int i=0; i<3; i++){
        for(int j = 0; j < SIZE_POLY; j++ ){
            for(int k= 0; k< SIZE_POLY; k++){
                pa(j,k)[i]= MAX_VALUE-1;                
    //            pb(j,k)[i]= MAX_VALUE-2;                
            }
        }
    }
    /*
    polynomial_gpu<vli_gpu<TYPE,8>, SIZE_POLY> pagpu(pa);
    polynomial_gpu<vli_gpu<TYPE,8>, SIZE_POLY> pbgpu(pb);
    polynomial_gpu<vli_gpu<TYPE,8>, SIZE_POLY> pcgpu(pa);

    Timer A("CPU");
    A.begin();
    pc = pa*pb;
    A.end();
    
    TimerCuda B("GPU");
    B.begin();
    pcgpu = pagpu*pbgpu;
    B.end();


    if(pc == pcgpu){ printf("ok \n");}else{printf("no ok \n");
    std::cout << pc << std::endl;
    std::cout << pcgpu << std::endl;
    }
*/
    
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,SIZE_POLY> > VaGPU(SIZE_VECTOR);
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,SIZE_POLY> > VbGPU(SIZE_VECTOR);
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,SIZE_POLY> > VcGPU(1);
    
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,SIZE_POLY> > VaCPU(SIZE_VECTOR);
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,SIZE_POLY> > VbCPU(SIZE_VECTOR);
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,SIZE_POLY> > VcCPU(1);
    
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
    
 
    if(VcGPU == VcCPU){ 
        printf("ok \n");
    }else{
        std::cout << VcCPU << std::endl;
        std::cout << VcGPU << std::endl;
        printf("no ok \n");
    }
  
    GPU->instance().destructor();
    return 0;
}
