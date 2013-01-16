#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_modulo )
{
    integer_type a,b;

    init(a);
    init(b);

    mpz_class agmp(a), bgmp(b);

    integer_type c = a % b;
    mpz_class cgmp = agmp % bgmp;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_modulo_negative )
{
    integer_type a,b;

    init(a);
    init(b);

    negate_inplace(a);

    mpz_class agmp(a), bgmp(b);

    integer_type c = a % b;
    mpz_class cgmp = agmp % bgmp;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}
