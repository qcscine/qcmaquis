#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_plus )
{
    integer_type a,b;

    init(a,overflow_safe);
    init(b,overflow_safe);

    mpz_class agmp(a), bgmp(b);

    integer_type c = a+b;
    mpz_class cgmp = agmp + bgmp;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}
/*
VLI_FUZZABLE_TEST( gmp_plus_int )
{
    integer_type a;
    init(a,overflow_safe);

    int b;
    init(b);

    mpz_class agmp(a);

    integer_type c = a+b;
    mpz_class cgmp = agmp + b;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}*/
