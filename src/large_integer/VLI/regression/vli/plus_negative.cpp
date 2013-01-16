#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( plus_assign_by_negative_number_minus_assign_equivalence )
{
    integer_type a,b;
    init(a);
    init(b);
    integer_type b_orig(b);
    integer_type c(a);

    a -= b;
    negate_inplace(b); negate_inplace(b_orig);
    c += b;

    BOOST_CHECK_EQUAL(a,c);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_FUZZABLE_TEST( plus_assign_by_negative_number_minus_assign_equivalence_int )
{
    integer_type a;
    int b;
    init(a);
    init(b);
    int b_orig(b);
    integer_type c(a);

    a -= b;
    c += (-b);

    BOOST_CHECK_EQUAL(a,c);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}
