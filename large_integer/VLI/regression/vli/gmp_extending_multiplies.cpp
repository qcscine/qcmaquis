#include <regression/vli/test_header.hpp>
#include <gmpxx.h>

using namespace vlilib::test;

VLI_FUZZABLE_TEST( gmp_extending_multiplies )
{
    vli_type a,b;
    init(a,max_positive);
    init(b,max_positive);
    vli_type a_orig(a);
    vli_type b_orig(b);

    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    typename double_sized_vli<vli_type>::type c; 
    mul(c,a,b); // TODO
    mpz_class cgmp = agmp * bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
    BOOST_CHECK_EQUAL(a,a_orig);
    BOOST_CHECK_EQUAL(b,b_orig);
}

