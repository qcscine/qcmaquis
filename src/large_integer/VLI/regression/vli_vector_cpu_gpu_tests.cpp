/*
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <cstdio>

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

#include "regression/common_test_functions.hpp"

using vli::vli_cpu;
using vli::max_int_value;
using vli::vli_gpu;
using vli::monomial;
using vli::polynomial_cpu;
using vli::polynomial_gpu;
using vli::vector_polynomial_gpu;
using vli::vector_polynomial_cpu;

#define SIZE 4
#define SIZE_VECTOR 64

typedef unsigned int TYPE;

BOOST_AUTO_TEST_CASE(vector_inner_product)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    
    polynomial_cpu<vli_cpu<TYPE,8>, 2> pa;
    vli::test::fill_poly_random(pa);
    
    polynomial_gpu<vli_gpu<TYPE,8>, 2> pagpu(pa);
 
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,2> > VaGPU(SIZE_VECTOR);
    vector_polynomial_gpu< polynomial_gpu<vli_gpu<TYPE, 8>,2> > VbGPU(SIZE_VECTOR);
    polynomial_gpu<vli_gpu<TYPE, 8>,2> pcGPU;

    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,2> > VaCPU(SIZE_VECTOR);
    vector_polynomial_cpu< polynomial_cpu<vli_cpu<TYPE, 8>,2> > VbCPU(SIZE_VECTOR);
    polynomial_cpu<vli_cpu<TYPE, 8>,2>  pcCPU;

    for(int i=0;i < SIZE;i++){
        VaCPU[i]=pa;
        VbCPU[i]=pa;

        VaGPU[i]=pa;
        VbGPU[i]=pa;
    }
     
    pcCPU =  inner_product(VaCPU,VbCPU);
    pcGPU =  inner_product(VaGPU,VbGPU);
    
     BOOST_CHECK_EQUAL(pcCPU,pcGPU);
     GPU->instance().destructor();
}


