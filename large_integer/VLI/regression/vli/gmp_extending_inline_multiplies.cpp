#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_extending_multiplies )
{
    integer_type a,b;
    init(a,max_positive);
    init(b,max_positive);
    integer_type a_orig(a);
    integer_type b_orig(b);

    mpz_class agmp(a), bgmp(b);

    typename double_sized_integer<integer_type>::type c;
    vli::detail::inline_mul_extend<integer_type::numwords>(&c[0],&a[0],&b[0]);
    
    mpz_class cgmp = agmp * bgmp;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
    BOOST_CHECK_EQUAL(a,a_orig);
    BOOST_CHECK_EQUAL(b,b_orig);
}

