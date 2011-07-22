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
#include "monome/vector_polynomial_cpu.h"
#include "monome/vector_polynomial_gpu.h"
#include "monome/polynome_gpu.h"
#include "monome/polynome_cpu.h"
#include "monome/monome.h"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "vli_gpu/vli_number_gpu.hpp"

using vli::vli_cpu;
using vli::max_int_value;
using vli::vli_gpu;
using vli::monomial;
using vli::polynomial_cpu;
using vli::polynomial_gpu;
using vli::vector_polynomial_gpu;
using vli::vector_polynomial_cpu;

#define SIZE 4

typedef unsigned int TYPE;

BOOST_AUTO_TEST_CASE(vector_inner_product)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    polynomial_cpu<vli_cpu<TYPE,8>, 2> pa;
    
    for(int i=0; i<2; i++){
        pa(0,0)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(0,1)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(1,0)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);
        pa(1,1)[i] = static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,8> >::value);        
    }
 
    polynomial_gpu<vli_gpu<TYPE,8>, 2> pagpu(pa);
 
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,2> > VaGPU(SIZE);
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,2> > VbGPU(SIZE);
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,2> > VcGPU;

    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,2> > VaCPU(SIZE);
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,2> > VbCPU(SIZE);
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,2> > VcCPU;

    for(int i=0;i < SIZE;i++){
        VaCPU[i]=pa;
        VbCPU[i]=pa;

        VaGPU[i]=pa;
        VbGPU[i]=pa;
 
    }
     
    VcCPU =  inner_product(VaCPU,VbCPU);
    VcGPU =  inner_product(VaGPU,VbGPU);
    
    BOOST_CHECK_EQUAL(VcGPU,VcCPU);
    GPU->instance().destructor();
}


