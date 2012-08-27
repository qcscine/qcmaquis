#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( negate )
{
    vli_type a;
    init(a);

    vli_type b(a);
    BOOST_CHECK_EQUAL(a.is_negative(), false);

    negate_inplace(a);
    BOOST_CHECK_EQUAL(a.is_negative(), true);
    BOOST_CHECK_EQUAL(a == b, false);

    negate_inplace(a);
    BOOST_CHECK_EQUAL(a.is_negative(), false);
    BOOST_CHECK_EQUAL(a,b);
}

VLI_FUZZABLE_TEST( negate_unary_minus_operator_equivalence )
{
    vli_type a;
    init(a);
    vli_type b(a);
    vli_type b_orig(b);
    negate_inplace(a);
    vli_type c(-b);
    BOOST_CHECK_EQUAL(a == b, false);
    BOOST_CHECK_EQUAL(c,a);
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_FUZZABLE_TEST( negate_and_construct_from_negative_int )
{
    int ai;
    init(ai);
    vli_type a(-ai);
    vli_type am(ai);
    negate_inplace(a);
    BOOST_CHECK_EQUAL(a,am);
}
