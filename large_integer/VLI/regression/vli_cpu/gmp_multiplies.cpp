#include <regression/vli_cpu/test_header.hpp>

#include <gmpxx.h>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_multiplies_negative )
{
    vli_type a(5);
    vli_type b(100);
    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);

    negate_inplace(a);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    vli_type c = a*b;
    vli_type cr = b*a;
    mpz_class cgmp = agmp * bgmp;
    
    negate_inplace(b);
    
    vli_type d = a*b;
    mpz_class dgmp = agmp * (-bgmp);
    
    BOOST_CHECK_EQUAL(cr,c);
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
    BOOST_CHECK_EQUAL(d.get_str(),dgmp.get_str());
}
