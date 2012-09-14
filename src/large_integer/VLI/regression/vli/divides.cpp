#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( multiply_divide )
{
    vli_type a;
    vli_type b;
    init(a);
    init(b);

    vli_type c(a);

    a /= b;

    a *= b;

    BOOST_CHECK_EQUAL( a <= c, true);
    BOOST_CHECK_EQUAL( a > c-b, true);
}


VLI_FUZZABLE_TEST( divide_assign_divide_equivalence )
{
    vli_type a;
    vli_type b;

    init(a);
    init(b);

    vli_type b_orig(b);

    vli_type ab = a / b;
    a /= b;
    BOOST_CHECK_EQUAL(a,ab);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}
