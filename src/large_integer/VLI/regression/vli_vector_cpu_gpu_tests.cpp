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

using vli::test::fill_random;
using vli::test::fill_poly_random;
using vli::test::fill_vector_random;


typedef vli_cpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value> vli_type_cpu;
typedef vli_gpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value> vli_type_gpu;

typedef vli::monomial<vli_cpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value> > monomial_type_cpu;
typedef vli::monomial<vli_gpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value> > monomial_type_gpu;

typedef vli::polynomial_cpu<vli_cpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value>, vli::detail::size_poly_vli::value > polynomial_type_cpu;
typedef vli::polynomial_gpu<vli_gpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value>, vli::detail::size_poly_vli::value > polynomial_type_gpu;

typedef vli::vector_polynomial_gpu<polynomial_type_gpu> vector_type_gpu;
typedef vli::vector_polynomial_cpu<polynomial_type_cpu> vector_type_cpu;


BOOST_AUTO_TEST_CASE(vector_copy)
{
/*
    bool test(false);        
    vector_type_cpu VaCPU( vli::detail::size_vector_vli::value); 
    
    fill_vector_random(VaCPU);
    vector_type_gpu VaGPU(VaCPU); ;
    

    if(VaGPU == VaCPU)
        test = true;

    BOOST_CHECK_EQUAL(test,true); //BUG BOOST CUDA MALLOC IF TEST ON THE VECTOR, DON'T KNOW WHY
*/
}


BOOST_AUTO_TEST_CASE(vector_inner_product)
{
    vector_type_cpu VaCPU( vli::detail::size_vector_vli::value),VbCPU(vli::detail::size_vector_vli::value); 
    
    polynomial_type_cpu pcCPU;
    polynomial_type_gpu pcGPU;
    
    fill_vector_random(VaCPU);
    fill_vector_random(VbCPU);

    vector_type_gpu VaGPU(VaCPU); 
    vector_type_gpu VbGPU(VbCPU); 
            
    pcCPU =  inner_product(VaCPU,VbCPU);
    pcGPU =  inner_product(VaGPU,VbGPU);
      
    BOOST_CHECK_EQUAL(pcCPU,pcGPU); 
}




