#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>


#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "gmpxx.h"

#include "regression/common_test_functions.hpp"

using vli::vli_cpu;
using vli::max_int_value;

using vli::test::rnd_valid_int;
using vli::test::fill_random;

typedef boost::mpl::list<
        vli_cpu<unsigned int,2>,
        vli_cpu<unsigned int,4>,
        vli_cpu<unsigned int,8>,
        vli_cpu<unsigned int,16>, 

        vli_cpu<unsigned long int,2>,
        vli_cpu<unsigned long int,4>, 
        vli_cpu<unsigned long int,8>,
        vli_cpu<unsigned long int,16>    
        > vli_types;

//I remove these tests from the general test case because we specified the type  for results

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_gmp, Vli, vli_types )
{
    Vli a;
    Vli b;
    fill_random(a,Vli::size);
    fill_random(b,Vli::size); 
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    vli_cpu<typename Vli::value_type,2*Vli::size> c = a*b;
    mpz_class cgmp = agmp * bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_negative_numbers_gmp, Vli, vli_types )
{
    Vli a;
    Vli b;
    fill_random(a,Vli::size);
    fill_random(b,Vli::size); 
    a.negate();
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    vli_cpu<typename Vli::value_type,2*Vli::size> c = a*b;
    mpz_class cgmp = agmp * bgmp;
    
    b.negate();
    vli_cpu<typename Vli::value_type,2*Vli::size> d = a*b;
    mpz_class dgmp = agmp * (-bgmp);
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
    BOOST_CHECK_EQUAL(d.get_str(),dgmp.get_str());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_reverse_negative_numbers_gmp, Vli, vli_types )
{
    Vli a;
    Vli b;    
    
    fill_random(a,Vli::size);
    fill_random(b,Vli::size); 
    b.negate();
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    vli_cpu<typename Vli::value_type,2*Vli::size> c = a*b;
    mpz_class cgmp = agmp * bgmp;
    
    b.negate();
    vli_cpu<typename Vli::value_type,2*Vli::size> d = a*b;
    mpz_class dgmp = agmp * (-bgmp);
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
    BOOST_CHECK_EQUAL(d.get_str(),dgmp.get_str());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_double_negative_numbers_gmp, Vli, vli_types )
{
    Vli a;
    Vli b;
    fill_random(a,Vli::size);
    fill_random(b,Vli::size); 
    a.negate();
    b.negate();
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    vli_cpu<typename Vli::value_type,2*Vli::size> c = a*b;
    mpz_class cgmp = agmp * bgmp;
    
    b.negate();
    vli_cpu<typename Vli::value_type,2*Vli::size> d = a*b;
    mpz_class dgmp = agmp * (-bgmp);
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
    BOOST_CHECK_EQUAL(d.get_str(),dgmp.get_str());
}

/**
    load all over tests
*/
#include "regression/vli_number_common_tests.hpp"

