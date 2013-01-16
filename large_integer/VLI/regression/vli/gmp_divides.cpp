#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_divides )
{
    integer_type a;
    integer_type b;
    init(a);
    init(b);
    mpz_class agmp(a), bgmp(b);
    integer_type c = a / b;
    mpz_class cgmp = agmp / bgmp;
    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_divides_negative )
{
    integer_type a;
    integer_type b;
    init(a);
    init(b);
    std::cout<<a<<std::endl;
    std::cout<<b<<std::endl;

    negate_inplace(a);

    mpz_class agmp(a), bgmp(b);

    integer_type c = a / b;
    mpz_class cgmp = agmp / bgmp;

    negate_inplace(b);

    integer_type d = a / b;
    mpz_class dgmp = agmp / (-bgmp);

    negate_inplace(a);

    integer_type e = a / b;
    mpz_class egmp = (-agmp) / (-bgmp);

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
    BOOST_CHECK_EQUAL(mpz_class(d),dgmp);
    BOOST_CHECK_EQUAL(mpz_class(e),egmp);
}
