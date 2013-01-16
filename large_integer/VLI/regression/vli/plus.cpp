#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( plus_assign_plus_equivalence )
{
    integer_type a;
    integer_type b;

    init(a);
    init(b);

    integer_type b_orig(b);

    integer_type ab = a + b;
    integer_type ba = b + a;
    a += b;
    BOOST_CHECK_EQUAL(a,ab);
    BOOST_CHECK_EQUAL(a,ba);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_FUZZABLE_TEST( plus_assign_plus_equivalence_int )
{
    integer_type a;
    init(a);
    int b;
    init(b);

    int b_orig(b);

    integer_type ab = a + b;
    integer_type ba = b + a;
    a += b;
    BOOST_CHECK_EQUAL(a,ab);
    BOOST_CHECK_EQUAL(a,ba);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_STATIC_TEST( plus_assign_carry )
{
    integer_type a,b,c;
    for(std::size_t i(0); i < integer_type::numwords-1; ++i)
        a[i] = static_cast<integer_type::value_type>(-1);
    b[0] = 1;
    a+=b;
    c[integer_type::numwords-1] = 1;
    BOOST_CHECK_EQUAL(a,c);
}

VLI_STATIC_TEST( plus_assign_overflow )
{
    // TODO this test could also be a fuzz test
    integer_type a;
    integer_type b(1);
    init(a,max_positive);
    integer_type c(a);
    negate_inplace(c);

    a+=b;
    a+=b;

    BOOST_CHECK_EQUAL(a,c);
}

