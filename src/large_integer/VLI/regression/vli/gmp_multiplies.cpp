#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_multiplies_negative )
{
    integer_type a;
    integer_type b;
    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);

    negate_inplace(a);

    mpz_class agmp(a), bgmp(b);

    integer_type c = a*b;
    integer_type cr = b*a;
    mpz_class cgmp = agmp * bgmp;

    negate_inplace(b);

    integer_type d = a*b;
    mpz_class dgmp = agmp * (-bgmp);

    BOOST_CHECK_EQUAL(cr,c);
    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
    BOOST_CHECK_EQUAL(mpz_class(d),dgmp);
}
