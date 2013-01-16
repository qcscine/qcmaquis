
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( minus_assign_minus_equivalence )
{
    integer_type a;
    integer_type b;
    init(a,max_positive);
    init(b,overflow_safe);
    integer_type b_orig(b);
    integer_type ab = a - b;
    integer_type ba = b - a;
    a -= b;
    BOOST_CHECK_EQUAL(a,ab);
    negate_inplace(a);
    BOOST_CHECK_EQUAL(a,ba);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_FUZZABLE_TEST( minus_assign_minus_equivalence_int )
{
    integer_type a;
    init(a);

    int b;
    init(b);
    int b_orig(b);

    integer_type ab = a - b;
    a -= b;
    BOOST_CHECK_EQUAL(a,ab);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);

    //b become negative
    b=-b;
    ab = a - b;
    a -= b;
    BOOST_CHECK_EQUAL(a,ab);

    BOOST_CHECK_EQUAL(b,-b_orig);
}

VLI_STATIC_TEST( minus_assign_borrow )
{
    integer_type a,b,c;

    a[integer_type::numwords-1] = 1;
    b[0] = 1;
    a-=b;

    for(std::size_t i(0); i < integer_type::numwords-1; ++i)
        c[i] = static_cast<integer_type::value_type>(-1);

    BOOST_CHECK_EQUAL(a,c);
}
