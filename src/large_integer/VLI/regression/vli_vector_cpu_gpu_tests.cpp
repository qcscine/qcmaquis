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
#include "vli/polynomial/polynomial_cpu.h"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.h"
#include "vli/vli_traits.hpp"

#include "use_gmp_integers.hpp"
#include "minimal_polynomial.hpp"

#include "regression/vli_test.hpp"

using vli::vli_cpu;
using vli::max_int_value;
using vli::monomial;
using vli::polynomial_cpu;
using vli::vector_polynomial_cpu;

using vli::test::fill_random;
using vli::test::fill_poly_random;
using vli::test::fill_vector_random;
using vli::test::fill_vector_negate;

typedef vli::test::vli_cpu_type_list vli_types;
typedef mpz_class large_int; //mpz_class into gmp
typedef hp2c::monomial<large_int> monomial_type;

typedef vli::test::vli_cpu_type_list vli_types;

enum { vector_size = 64 };

BOOST_AUTO_TEST_CASE_TEMPLATE(vector_inner_product_gmp_positive, Vli, vli_types)
{
    // VLI
    typedef vli::polynomial_cpu<Vli, 11 > polynomial_type_cpu;
    typedef vli::vli_cpu<typename Vli::value_type,  2*Vli::size > vli_result_type_cpu;
    typedef vli::polynomial_cpu<vli_result_type_cpu, 22 > polynomial_result_type_cpu;
    typedef vli::vector_polynomial_cpu<polynomial_type_cpu> vector_type_cpu;
    // GMP
    typedef hp2c::polynomial<large_int,11> poly_gmp;
    typedef hp2c::polynomial<large_int,2*11> poly_gmp_double;
    typedef std::vector<poly_gmp> vector_poly_gmp;
    
    poly_gmp_double pgmpd;
    polynomial_result_type_cpu pd; 
   
    vector_type_cpu v1(vector_size);
    vector_type_cpu v2(vector_size);
   
    vector_poly_gmp vgmp1(vector_size);
    vector_poly_gmp vgmp2(vector_size);
   
    fill_vector_random(v1);
    fill_vector_random(v2);
    
    vli::test::InitVecVLItoVecGMP(v1,vgmp1);
    vli::test::InitVecVLItoVecGMP(v2,vgmp2);
    pd = vli::detail::inner_product_plain(v1,v2);
    pgmpd = inner_product(vgmp1,vgmp2);
    
    BOOST_CHECK_EQUAL(vli::test::ValidatePolyVLI_PolyGMP(pd,pgmpd), true );
}


BOOST_AUTO_TEST_CASE_TEMPLATE(vector_inner_product_gmp_positive_negative, Vli, vli_types)
{
    // VLI
    typedef vli::polynomial_cpu<Vli, 11 > polynomial_type_cpu;
    typedef vli::vli_cpu<typename Vli::value_type,  2*Vli::size > vli_result_type_cpu;
    typedef vli::polynomial_cpu<vli_result_type_cpu, 22 > polynomial_result_type_cpu;
    typedef vli::vector_polynomial_cpu<polynomial_type_cpu> vector_type_cpu;
    // GMP
    typedef hp2c::polynomial<large_int,11> poly_gmp;
    typedef hp2c::polynomial<large_int,2*11> poly_gmp_double;
    typedef std::vector<poly_gmp> vector_poly_gmp;
    
    poly_gmp_double pgmpd;
    polynomial_result_type_cpu pd; 
    
    vector_type_cpu v1(vector_size);
    vector_type_cpu v2(vector_size);
    
    vector_poly_gmp vgmp1(vector_size);
    vector_poly_gmp vgmp2(vector_size);
    
    fill_vector_random(v1,2);
    fill_vector_random(v2,3);
    
    fill_vector_negate(v1,2);
    fill_vector_negate(v2,3);
        
    vli::test::InitVecVLItoVecGMP(v1,vgmp1);
    vli::test::InitVecVLItoVecGMP(v2,vgmp2);
    pd = vli::detail::inner_product_plain(v1,v2);
    pgmpd = inner_product(vgmp1,vgmp2);
        
    BOOST_CHECK_EQUAL(vli::test::ValidatePolyVLI_PolyGMP(pd,pgmpd), true );
}
