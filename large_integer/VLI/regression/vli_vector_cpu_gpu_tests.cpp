/*
 *
 *  Created by Tim Ewart on 18.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#define BOOST_TEST_MODULE vli_vector_cpu_gpu_tests
#include <boost/test/unit_test.hpp>
#include <boost/mpl/front.hpp>
#include <cstdio>

#include "vli/utils/gpu_manager.h"
#include "vli/utils/gpu_manager.hpp"
#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/vector_polynomial_gpu.hpp"
#include "vli/polynomial/polynomial_gpu.hpp"
#include "vli/polynomial/polynomial_cpu.hpp"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.hpp"
#include "vli/vli_gpu.hpp"
#include "vli/vli_traits.hpp"

#include "regression/vli_test.hpp"

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

template <typename Vli>
struct vli_gpu_from_vli_cpu
{
    typedef vli_gpu<typename Vli::value_type, Vli::size> type;
};

typedef boost::mpl::front<vli::test::vli_cpu_type_list>::type vli_type_cpu;
typedef vli_gpu_from_vli_cpu<vli_type_cpu>::type vli_type_gpu;

typedef vli::monomial<vli_type_cpu> monomial_type_cpu;
typedef vli::monomial<vli_type_gpu> monomial_type_gpu;

typedef vli::polynomial_cpu<vli_type_cpu, vli::detail::size_poly_vli::value > polynomial_type_cpu;
typedef vli::polynomial_gpu<vli_type_gpu, vli::detail::size_poly_vli::value > polynomial_type_gpu;

typedef vli::vector_polynomial_gpu<polynomial_type_gpu> vector_type_gpu;
typedef vli::vector_polynomial_cpu<polynomial_type_cpu> vector_type_cpu;

enum { vector_size = 100 };

BOOST_AUTO_TEST_CASE(vector_copy)
{
    vector_type_cpu VaCPU(vector_size); 
    
    fill_vector_random(VaCPU);
    vector_type_gpu VaGPU(VaCPU); ;

    BOOST_CHECK_EQUAL(VaGPU,VaCPU); //BUG BOOST CUDA MALLOC IF TEST ON THE VECTOR, DON'T KNOW WHY
}

BOOST_AUTO_TEST_CASE(vector_polys)
{
    vector_type_cpu VaCPU(vector_size); 
    
    fill_vector_random(VaCPU);

    monomial_type_cpu MaCPU;   
     




    vector_type_gpu VaGPU(VaCPU); ;

    BOOST_CHECK_EQUAL(VaGPU,VaCPU); //BUG BOOST CUDA MALLOC IF TEST ON THE VECTOR, DON'T KNOW WHY
}

BOOST_AUTO_TEST_CASE(vector_inner_product)
{
    vector_type_cpu VaCPU(vector_size),VbCPU(vector_size);
    
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




