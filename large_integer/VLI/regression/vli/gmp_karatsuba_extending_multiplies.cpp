#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_karatsuba_extending_multiplies)
{
    integer_type a,b;
    init(a,max_positive);
    init(b,max_positive);
    integer_type a_orig(a);
    integer_type b_orig(b);

    mpz_class agmp(a), bgmp(b);

    typename double_sized_integer<integer_type>::type c;
    vli::detail::karatsuba(c,a,b);
    mpz_class cgmp = agmp * bgmp;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
    BOOST_CHECK_EQUAL(a,a_orig);
    BOOST_CHECK_EQUAL(b,b_orig);
}

