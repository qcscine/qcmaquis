#include <regression/vli/test_header.hpp>
#include<gmpxx.h>

using namespace vlilib::test;

VLI_FUZZABLE_TEST( gmp_minus )
{
    vli_type a;
    vli_type b;
    init(a,max_positive);
    init(b,overflow_safe);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    vli_type c = a-b;
    mpz_class cgmp = agmp - bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_minus_int )
{
    vli_type a;
    init(a,max_positive);
    int b;
    init(b);
    
    mpz_class agmp(a.get_str());
    
    vli_type c = a-b;
    mpz_class cgmp = agmp - b;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}
