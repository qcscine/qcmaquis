#include <gmpxx.h>
#include <regression/vli_cpu/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_and_str_conversion_mpz )
{
    vli_type a;
    init(a);
    mpz_class agmp(a);
    BOOST_CHECK_EQUAL(a.get_str(),agmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_and_str_conversion_mpq )
{
    vli_type a;
    init(a);
    mpq_class agmp(a);
    BOOST_CHECK_EQUAL(a.get_str(),agmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_and_str_conversion_negative_mpz )
{
    vli_type a;
    init(a);
    negate_inplace(a);
    mpz_class agmp(a);
    BOOST_CHECK_EQUAL(a.get_str(),agmp.get_str());
}


VLI_FUZZABLE_TEST( gmp_and_str_conversionnegative_mpq )
{
    vli_type a;
    init(a);
    negate_inplace(a);
    mpq_class agmp(a);
    BOOST_CHECK_EQUAL(a.get_str(),agmp.get_str());
}
