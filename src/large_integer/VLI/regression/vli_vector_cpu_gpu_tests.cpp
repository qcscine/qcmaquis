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

#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial_cpu.hpp"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.hpp"
#include "vli/vli_traits.hpp"

#include "regression/vli_test.hpp"

using vli::vli_cpu;
using vli::max_int_value;
using vli::monomial;
using vli::polynomial_cpu;
using vli::vector_polynomial_cpu;

using vli::test::fill_random;
using vli::test::fill_poly_random;
using vli::test::fill_vector_random;

typedef vli::test::vli_cpu_type_list vli_types;

enum { vector_size = 100 };

BOOST_AUTO_TEST_CASE_TEMPLATE(vector_inner_product, Vli, vli_types)
{
    typedef vli::polynomial_cpu<Vli, 13 > polynomial_type_cpu;
    typedef vli::polynomial_cpu<Vli, 26 > polynomial_result_type_cpu;
    typedef vli::vector_polynomial_cpu<polynomial_type_cpu> vector_type_cpu;

    vector_type_cpu VaCPU(vector_size),VbCPU(vector_size);
    polynomial_result_type_cpu pcCPU0;
    fill_vector_random(VaCPU);
    fill_vector_random(VbCPU);
    pcCPU0 = inner_product(VaCPU,VbCPU);

    polynomial_result_type_cpu pcCPU1;
    pcCPU1 = vli::detail::inner_product_plain(VaCPU,VbCPU);
    BOOST_CHECK_EQUAL(pcCPU0,pcCPU1);
}




