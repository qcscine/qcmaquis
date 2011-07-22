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
#include "vli_cpu/vli_number_traits.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "monome/monome.h"
#include "monome/polynome_cpu.h"
#include "monome/polynome_gpu.h"

//#include "utils/timings.h"

typedef unsigned int TYPE; 
using vli::vli_cpu;
using vli::max_int_value;
using vli::vli_gpu;
using vli::monomial;
using vli::polynomial_cpu;
using vli::polynomial_gpu;

BOOST_AUTO_TEST_CASE(monome)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    vli_cpu<TYPE,8> a;
    a[0] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    a[1]  =static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    vli_cpu<TYPE,8> b(a);
    a[0] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    a[1] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);

    vli_gpu<TYPE,8> agpu(a);
    vli_gpu<TYPE,8> bgpu(b);

    monomial<vli_cpu<TYPE, 8> > ma(a);
    monomial<vli_cpu<TYPE, 8> > mb(b);

    monomial<vli_gpu<TYPE, 8> > magpu(agpu);
    monomial<vli_gpu<TYPE, 8> > mbgpu(bgpu);

    ma*=b;    
    magpu *= bgpu;    

    std::cout << ma << std::endl;
    std::cout << magpu << std::endl;

    BOOST_CHECK_EQUAL(ma.coeff_,magpu.coeff_);
	
    GPU->instance().destructor();
}
BOOST_AUTO_TEST_CASE(monome_polynome)
{
  	enum {vli_size = 8};
    gpu::gpu_manager* GPU;
	GPU->instance();
    
    //init vli cpu    
    vli_cpu<TYPE,vli_size> a;
    vli_cpu<TYPE,vli_size> b;

    b[0] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    b[1] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    a[0] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    a[1]  =static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    //init vli gpu
    vli_cpu<TYPE,vli_size> agpu(a);
    vli_cpu<TYPE,vli_size> bgpu(b);
    //TYPEi momone cpu
    monomial<vli_cpu<TYPE, vli_size> > mb(b);
    monomial<vli_cpu<TYPE, vli_size> > ma(a);
    //init monome gpu
    monomial<vli_gpu<TYPE, vli_size> > mbgpu(b);
    monomial<vli_gpu<TYPE, vli_size> > magpu(a);    
    //init polynome cpu
    polynomial_cpu<vli_cpu<TYPE,vli_size>, 20> pa;
    polynomial_cpu<vli_cpu<TYPE,vli_size>, 20> pb;
    
    for(int i=0; i<vli_size; i++){
        pa(0,0)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(0,1)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(1,0)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(1,1)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);        
    }
    //init polynome gpu
    polynomial_gpu<vli_gpu<TYPE,vli_size>, 20> pagpu(pa);
    polynomial_gpu<vli_gpu<TYPE,vli_size>, 20> pbgpu(pb);
    
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
    vli_cpu<TYPE,8> a;
    vli_cpu<TYPE,8> b;
    a[0] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    a[1]  =static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    b[0] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    b[1] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
    //init vli gpu
    vli_gpu<TYPE,8> agpu(a);
    vli_gpu<TYPE,8> bgpu(b);
    //init polynome cpu
    polynomial_cpu<vli_cpu<TYPE,8>, 2> pa;
    polynomial_cpu<vli_cpu<TYPE,8>, 2> pb;
    
    for(int i=0; i<2; i++){
        pa(0,0)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(0,1)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(1,0)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(1,1)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);        
    }
    //init polynome gpu
    polynomial_gpu<vli_gpu<TYPE,8>, 2> pagpu(pa);
    polynomial_gpu<vli_gpu<TYPE,8>, 2> pbgpu(pb);


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
