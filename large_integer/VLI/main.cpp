
#include <iostream>

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"

#include "monome/vector_polynomial_gpu.h"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "monome/monome.h"
#include "monome/polynome_gpu.h"
#include "monome/polynome.h"

using vli::vli_cpu;
using vli::vli_gpu;
using vli::monomial;
using vli::polynomial;
using vli::polynomial_gpu;
using vli::vector_polynomial_gpu;

int main (int argc, char * const argv[]) 
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    polynomial<vli_cpu<int,8>, 2> pa;
    
    for(int i=0; i<2; i++){
        pa(0,0)[i] = 255;
        pa(0,1)[i] = 255;
        pa(1,0)[i] = 255;
        pa(1,1)[i] = 255;        
    }
 
    polynomial_gpu<vli_gpu<int,8>, 2> pagpu(pa);
 
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<int, 8>,2> > V(4);
    
    V[0] = pagpu;
    
     std::cout << V << std::endl;
    
   
	GPU->instance().destructor();
}
