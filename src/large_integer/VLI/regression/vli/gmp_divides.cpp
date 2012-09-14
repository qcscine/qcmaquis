#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_divides )
{
    vli_type a;
    vli_type b;
    init(a);
    init(b);
    mpz_class agmp(a), bgmp(b);
    vli_type c = a / b;
    mpz_class cgmp = agmp / bgmp;
    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_divides_negative )
{
    vli_type a;
    vli_type b;
    init(a);
    init(b);
    std::cout<<a<<std::endl;
    std::cout<<b<<std::endl;

    negate_inplace(a);

    mpz_class agmp(a), bgmp(b);

    vli_type c = a / b;
    mpz_class cgmp = agmp / bgmp;

    negate_inplace(b);

    vli_type d = a / b;
    mpz_class dgmp = agmp / (-bgmp);

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
    BOOST_CHECK_EQUAL(mpz_class(d),dgmp);
}
