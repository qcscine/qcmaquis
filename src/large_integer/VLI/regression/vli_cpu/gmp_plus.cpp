#include <regression/vli_cpu/test_header.hpp>
#include <gmpxx.h>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_plus )
{
    vli_type a,b;

    init(a,overflow_safe);
    init(b,overflow_safe);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    vli_type c = a+b;
    mpz_class cgmp = agmp + bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_plus_int )
{
    vli_type a;
    init(a,overflow_safe);

    int b;
    init(b);
    
    mpz_class agmp(a.get_str());
    
    vli_type c = a+b;
    mpz_class cgmp = agmp + b;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}
