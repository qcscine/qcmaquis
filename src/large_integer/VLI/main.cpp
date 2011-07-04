
//#include "vli_config.h"
#include <iostream>

#define SIZE_BITS 256


#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "monome/monome.h"
#include "monome/polynome.h"
#include "monome/polynome_gpu.h"

#include "utils/timings.h"

typedef int TYPE; 
using vli::vli_cpu;
using vli::vli_gpu;
using vli::monomial;
using vli::polynomial;
using vli::polynomial_gpu;


int main (int argc, char * const argv[]) 
{
 	gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<int,8> a;
    a[0] = 255;
    a[1] = 255;
    
    vli_cpu<int,8> b(a);
    vli_gpu<int,8> c(a);
     
   
    

	a+=a;
    std::cout << a <<std::endl;
    
    vli::monomial<vli_cpu<int, 8> > ma(a);
    ma*=a;

    vli::monomial<vli_gpu<int, 8> > magpu(c);
    
    //    ma[0] = 255;
 
    
    
    polynomial<vli_cpu<int,8>, 2> pa;
    polynomial<vli_cpu<int,8>, 2> pb;
    polynomial<vli_cpu<int,8>, 2> pc;
    
    for(int i=0; i<2; i++){
        pa(0,0)[i] = 255;
        pa(0,1)[i] = 255;
        pa(1,0)[i] = 255;
        pa(1,1)[i] = 255;        
        
        pb(0,0)[i] = 255;
        pb(0,1)[i] = 255;
        pb(1,0)[i] = 255;
        pb(1,1)[i] = 255;                
    }
    
    polynomial_gpu<vli_gpu<int,8>, 2> pgpua(pa);
    polynomial_gpu<vli_gpu<int,8>, 2> pgpub(pb);
    polynomial_gpu<vli_gpu<int,8>, 2> pgpuc(pc);
    
    pa += pa;
    
    pc = pa*ma;
    pc = ma*pa;
    pc *= ma;
    
    pc = pa*a;
    pc = a*pa;
    pgpua = pgpua*magpu;
    
    
    printf(" Pc=Pa*pb \n");
    printf("CPU \n");
    std::cout << pa(0,0) << std::endl; 
    std::cout << pa(0,1) << std::endl; 
    std::cout << pa(1,0) << std::endl; 
    std::cout << pa(1,1) << std::endl; 
    printf("---------------------------\n");
    printf("GPU \n");
    std::cout << pgpua(0,0) << std::endl; 
    std::cout << pgpua(0,1) << std::endl; 
    std::cout << pgpua(1,0) << std::endl; 
    std::cout << pgpua(1,1) << std::endl; 
    
    
	GPU->instance().destructor();
    return 0;
}

