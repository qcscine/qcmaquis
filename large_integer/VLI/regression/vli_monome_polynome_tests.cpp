/*
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <cstdio>

#define SIZE_BITS 256


#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "monome/monome.h"
#include "monome/polynome_cpu.h"
#include "monome/polynome_gpu.h"

//#include "utils/timings.h"

typedef int TYPE; 
using vli::vli_cpu;
using vli::vli_gpu;
using vli::monomial;
using vli::polynomial_cpu;
using vli::polynomial_gpu;

BOOST_AUTO_TEST_CASE(monome)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<int,8> a;
    a[0] = 185;
    a[1]  =254;
    vli_cpu<int,8> b(a);
    a[0] = 255;
    a[1] = 255;

    vli_gpu<int,8> agpu(a);
    vli_gpu<int,8> bgpu(b);

    monomial<vli_cpu<int, 8> > ma(a);
    monomial<vli_cpu<int, 8> > mb(b);

    monomial<vli_gpu<int, 8> > magpu(agpu);
    monomial<vli_gpu<int, 8> > mbgpu(bgpu);

    ma*=b;    
    magpu *= bgpu;    

    std::cout << ma << std::endl;
    std::cout << magpu << std::endl;

    BOOST_CHECK_EQUAL(ma.coeff,magpu.coeff);
}

BOOST_AUTO_TEST_CASE(monome_polynome)
{
  	gpu::gpu_manager* GPU;
	GPU->instance();
    //init vli cpu    
    vli_cpu<int,8> a;
    vli_cpu<int,8> b;
    b[0] = 255;
    b[1] = 255;
    a[0] = 185;
    a[1]  =254;
    //init vli gpu
    vli_cpu<int,8> agpu(a);
    vli_cpu<int,8> bgpu(b);
    //inti momone cpu
    monomial<vli_cpu<int, 8> > mb(b);
    monomial<vli_cpu<int, 8> > ma(a);
    //init monome gpu
    monomial<vli_gpu<int, 8> > mbgpu(b);
    monomial<vli_gpu<int, 8> > magpu(a);    
    //init polynome cpu
    polynomial_cpu<vli_cpu<int,8>, 2> pa;
    polynomial_cpu<vli_cpu<int,8>, 2> pb;
    
    for(int i=0; i<2; i++){
        pa(0,0)[i] = 255;
        pa(0,1)[i] = 255;
        pa(1,0)[i] = 255;
        pa(1,1)[i] = 255;        
    }
    //init polynome gpu
    polynomial_gpu<vli_gpu<int,8>, 2> pagpu(pa);
    polynomial_gpu<vli_gpu<int,8>, 2> pbgpu(pb);
    
    std::cout << pa << std::endl;
    std::cout << pagpu << std::endl;
    std::cout << ma << std::endl;
    std::cout << magpu << std::endl;
    
    
    
    pb = pa*ma;        
    pbgpu = pagpu*magpu;
    
    
    std::cout << pb << std::endl;
    std::cout << pbgpu << std::endl;
    
    printf(" Pb=Pa*ma \n");
    printf("CPU \n");
    std::cout << pb(0,0) << std::endl; 
    std::cout << pb(0,1) << std::endl; 
    std::cout << pb(1,0) << std::endl; 
    std::cout << pb(1,1) << std::endl; 
    printf("---------------------------\n");
    printf("GPU \n");
    std::cout << pbgpu(0,0) << std::endl; 
    std::cout << pbgpu(0,1) << std::endl; 
    std::cout << pbgpu(1,0) << std::endl; 
    std::cout << pbgpu(1,1) << std::endl; 
    
    BOOST_CHECK_EQUAL(pb,pbgpu);
	GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE(monome_polynome_deux)
{
    gpu::gpu_manager* GPU;
	GPU->instance();
    //init vli cpu    
    vli_cpu<int,8> a;
    vli_cpu<int,8> b;
    a[0] = 185;
    a[1]  =254;
    b[0] = 255;
    b[1] = 255;
    //init vli gpu
    vli_gpu<int,8> agpu(a);
    vli_gpu<int,8> bgpu(b);
    //init polynome cpu
    polynomial_cpu<vli_cpu<int,8>, 2> pa;
    polynomial_cpu<vli_cpu<int,8>, 2> pb;
    
    for(int i=0; i<2; i++){
        pa(0,0)[i] = 255;
        pa(0,1)[i] = 255;
        pa(1,0)[i] = 255;
        pa(1,1)[i] = 255;        
    }
    //init polynome gpu
    polynomial_gpu<vli_gpu<int,8>, 2> pagpu(pa);
    polynomial_gpu<vli_gpu<int,8>, 2> pbgpu(pb);


    pb = pa*a;
    pbgpu = pagpu*agpu;
    
    
    printf(" Pb=Pa*ma \n");
    printf("CPU \n");
    std::cout << pb(0,0) << std::endl; 
    std::cout << pb(0,1) << std::endl; 
    std::cout << pb(1,0) << std::endl; 
    std::cout << pb(1,1) << std::endl; 
    printf("---------------------------\n");
    printf("GPU \n");
    std::cout << pbgpu(0,0) << std::endl; 
    std::cout << pbgpu(0,1) << std::endl; 
    std::cout << pbgpu(1,0) << std::endl; 
    std::cout << pbgpu(1,1) << std::endl; 
    
    BOOST_CHECK_EQUAL(pb,pbgpu);
	GPU->instance().destructor();
    
}

