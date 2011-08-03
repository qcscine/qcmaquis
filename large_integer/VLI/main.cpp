#include <iostream>
#include <cstdio>

#include "gmpxx.h"
#include <boost/lexical_cast.hpp>

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"

#include "polynomial/vector_polynomial_cpu.hpp"
#include "polynomial/vector_polynomial_gpu.hpp"
#include "polynomial/polynomial_gpu.hpp"
#include "polynomial/polynomial_cpu.hpp"
#include "polynomial/monomial.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "utils/timings.h"
#define SIZE_POLY 8
#define SIZE_VECTOR 512

#define TYPE unsigned  int 

using vli::vli_cpu;
using vli::max_int_value;
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
    
    
/*
    polynomial_cpu<vli_cpu<TYPE,8>, SIZE_POLY> pa;

    
    for(int i=0; i< 8 ; i++){
        for(int j=0; j < SIZE_POLY; j++){
            for(int k=0; k < SIZE_POLY; k++)
                pa(j,k)[i] = 9;
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

    TimerCuda C("GPU");
    C.begin();     
    VcGPU =  inner_product(VaGPU,VbGPU);
    C.end();

    if(VcGPU == VcCPU){
        std::cout << " ok " << std::endl; 
    }else{
        std::cout << VcCPU << std::endl;
        std::cout << VcGPU << std::endl;
    }
 
 */

    vli_cpu<TYPE,8> a;
    vli_cpu<TYPE,8> b;
    vli_cpu<TYPE,8> c;

    a[0] = 3;
    b[0] = 3;    
    c[0] = -3;
    
    b.negate();

    std::cout << b << std::endl;
    std::cout << c << std::endl;

    a*=b;
    
  //  std::cout << a << std::endl;
    
    
    return 0;
}
