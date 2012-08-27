#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( plus_assign_extend_positive )
{
    vli_type a;
    vli_type b;
    init(a,max_positive);
    init(b,max_positive);
    mpz_class agmp(a), bgmp(b);

    typename extended<vli_type>::type c = plus_extend(a,b);
    mpz_class cgmp = agmp + bgmp;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( plus_assign_extend_negative )
{
    vli_type a,b;
    init(a);
    init(b);
    negate_inplace(a);
    negate_inplace(b);
    mpz_class agmp(a), bgmp(b);

    typename extended<vli_type>::type c = plus_extend(a,b);
    mpz_class cgmp = agmp + bgmp;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( plus_assign_extend_mixed_sign )
{
    vli_type a,b;
    init(a,max_positive);
    init(b);
    negate_inplace(b);
    mpz_class agmp(a), bgmp(b);

    typename extended<vli_type>::type c = plus_extend(a,b);
    typename extended<vli_type>::type d = plus_extend(b,a);
    mpz_class cgmp = agmp + bgmp;

    BOOST_CHECK_EQUAL(c,d);
    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

