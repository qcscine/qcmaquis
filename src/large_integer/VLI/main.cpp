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
#include "vli_cpu/vli_number_traits.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "utils/timings.h"

#define SIZE_POLY 20
#define SIZE_VECTOR 131072

#define TYPE unsigned long int 

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

    polynomial_cpu<vli_cpu<TYPE,8>, 2> pa;
    
    for(int i=0; i<2; i++){
        pa(0,0)[i] = static_cast<TYPE>(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(0,1)[i] = static_cast<TYPE>(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(1,0)[i] = static_cast<TYPE>(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(1,1)[i] = static_cast<TYPE>(max_int_value<vli_cpu<TYPE,8> >::value);        
    }
 
    polynomial_gpu<vli_gpu<TYPE,8>, 2> pagpu(pa);
 
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,2> > VaGPU(SIZE_VECTOR);
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,2> > VbGPU(SIZE_VECTOR);
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,2> > VcGPU;

    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,2> > VaCPU(SIZE_VECTOR);
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,2> > VbCPU(SIZE_VECTOR);
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,2> > VcCPU;

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


    Timer B("GPU");
    B.begin();     
    VcGPU =  inner_product(VaGPU,VbGPU);
    B.end();

    if(VcGPU == VcCPU){
        std::cout << " ok " << std::endl; 
    }else{
        std::cout << VcCPU << std::endl;
        std::cout << VcGPU << std::endl;
    }


    GPU->instance().destructor();
    return 0;
}
