#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_minus )
{
    integer_type a;
    integer_type b;
    init(a,max_positive);
    init(b,overflow_safe);

    mpz_class agmp(a), bgmp(b);

    integer_type c = a-b;
    mpz_class cgmp = agmp - bgmp;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_minus_int )
{
    integer_type a;
    init(a,max_positive);
    int b;
    init(b);

    mpz_class agmp(a);

    integer_type c = a-b;
    mpz_class cgmp = agmp - b;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}
