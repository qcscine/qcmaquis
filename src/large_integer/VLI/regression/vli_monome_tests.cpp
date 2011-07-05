#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <cstdio>

#define SIZE_BITS 256

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "monome/monome.h"
#include "monome/polynome_gpu.h"
#include "monome/polynome.h"


#include "utils/timings.h"

typedef int TYPE; 
using vli::vli_cpu;
using vli::vli_gpu;
using vli::polynomial;
using vli::monomial;
#define SIZE 8
#define SIZE_VECTOR 8

BOOST_AUTO_TEST_CASE( constructors_test_site_monome)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
	
    vli::polynomial<vli_cpu<int,8>,2> pa;
    vli::polynomial<vli_cpu<int,8>,2> pb;
    vli::polynomial<vli_cpu<int,8>,2> pc;
    
    vli::polynomial_gpu<vli_gpu<int,8>,2> pagpu;
    vli::polynomial_gpu<vli_gpu<int,8>,2> pcgpu;
    
    vli::polynomial_gpu<vli_gpu<int,8>,2> pdgpu(pc);
    
    for(int i=0; i <2; i++){     
        
        pa(0,0)[i] = 255;        
        pa(0,1)[i] = 255;
        pa(1,0)[i] = 255;
        pa(1,1)[i] = 255;
        
        pb(0,0)[i] = 255;
        pb(0,1)[i] = 255;
        pb(1,0)[i] = 255;
        pb(1,1)[i] = 255;
        
        pagpu(0,0) = pa(0,0);        
        pagpu(0,1) = pa(0,1);
        pagpu(1,0) = pa(1,0);
        pagpu(1,1) = pa(1,1);        
    }
    
    
    vli::polynomial_gpu<vli_gpu<int,8>,2> pbgpu(pb);
    
    
    /** THE ORDER IS VERY IMPORTANT, first CPU second GPU*/
    BOOST_CHECK_EQUAL(pa,pagpu);
    BOOST_CHECK_EQUAL(pb,pbgpu);
    BOOST_CHECK_EQUAL(pc,pcgpu);

	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE( multiplication_polynomial)
{
  	gpu::gpu_manager* GPU;
	GPU->instance();
        
    vli::polynomial<vli_cpu<int,8>,2> pa;
    vli::polynomial<vli_cpu<int,8>,2> pb;
    vli::polynomial<vli_cpu<int,8>,2> pc;
    
    vli::polynomial_gpu<vli_gpu<int,8>,2> pagpu;
    vli::polynomial_gpu<vli_gpu<int,8>,2> pbgpu;
    vli::polynomial_gpu<vli_gpu<int,8>,2> pcgpu;
    
    
    for(int i=0; i <2; i++){     
        
        pa(0,0)[i] = 255;        
        pa(0,1)[i] = 255;
        pa(1,0)[i] = 255;
        pa(1,1)[i] = 255;
        
        pb(0,0)[i] = 255;
        pb(0,1)[i] = 255;
        pb(1,0)[i] = 255;
        pb(1,1)[i] = 255;
        
        pagpu(0,0)[i] = 255;        
        pagpu(0,1)[i] = 255;
        pagpu(1,0)[i] = 255;
        pagpu(1,1)[i] = 255;
        
        pbgpu(0,0) = pb(0,0);
        pbgpu(0,1) = pb(0,1);
        pbgpu(1,0) = pb(1,0);
        pbgpu(1,1) = pb(1,1); 
    }
    
    Timer A("CPU");
    A.begin();
    pc = pa*pb;
    A.end();
    
    Timer B("GPU classic");
    B.begin();
    TimerCuda C("GPU Cuda");
    B.begin();
    C.begin();    
    pcgpu = pagpu*pbgpu;
    B.end();
    C.end();
    
    printf("GPU \n");
    std::cout << pcgpu << std::endl;
    printf("--------------------------- \n");
    printf("CPU \n");
    std::cout << pc << std::endl;
    
    /** THE ORDER IS VERY IMPORTANT, first CPU second GPU*/
    BOOST_CHECK_EQUAL(pc,pcgpu);
    
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE( addition_polynomial)
{
    gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli::polynomial<vli_cpu<int,8>,2> pa;
    vli::polynomial<vli_cpu<int,8>,2> pb;
    vli::polynomial<vli_cpu<int,8>,2> pc;
    
    for(int i=0; i <2; i++){     
        
        pa(0,0)[i] = 255;        
        pa(0,1)[i] = 255;
        pa(1,0)[i] = 255;
        pa(1,1)[i] = 255;
        
        pb(0,0)[i] = 255;
        pb(0,1)[i] = 255;
        pb(1,0)[i] = 255;
        pb(1,1)[i] = 255;
    }
    
    vli::polynomial_gpu<vli_gpu<int,8>,2> pagpu(pa);
    vli::polynomial_gpu<vli_gpu<int,8>,2> pbgpu(pb);
    vli::polynomial_gpu<vli_gpu<int,8>,2> pcgpu(pb);

    pa+=pb;
    pagpu+=pbgpu;
    
    printf("GPU \n");
    std::cout << pagpu << std::endl;
    printf("--------------------------- \n");
    printf("CPU \n");
    std::cout << pa << std::endl;
    /** THE ORDER IS VERY IMPORTANT, first CPU second GPU*/
    BOOST_CHECK_EQUAL(pa,pagpu);
    
	GPU->instance().destructor();
    
}
