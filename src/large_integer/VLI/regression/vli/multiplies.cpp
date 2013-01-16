#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( multiplies_assign_multiplies_equivalence )
{
    integer_type a;
    integer_type b;
    init(a);
    init(b);
    integer_type b_orig(b);

    integer_type ab = a*b;
    integer_type ba = b*a;
    a*=b;

    BOOST_CHECK_EQUAL(a,ab);
    BOOST_CHECK_EQUAL(a,ba);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_FUZZABLE_TEST( multiplies_assign_multiplies_equivalence_int )
{
    integer_type a;
    int b;
    init(a);
    init(b);

    int b_orig(b);

    integer_type ab = a*b;
    integer_type ba = b*a;
    a*=b;

    BOOST_CHECK_EQUAL(a,ab);
    BOOST_CHECK_EQUAL(a,ba);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_STATIC_TEST( multiplies_plus_relation_int )
{
    integer_type a;
    init(a,overflow_safe);
    integer_type a_orig(a);

    integer_type b = a+a+a;
    integer_type c = a * 3;

    BOOST_CHECK_EQUAL(c,b);

    //Check that a hasn't changed
    BOOST_CHECK_EQUAL(a,a_orig);
}

VLI_STATIC_TEST( multiplies_plus_relation )
{
    integer_type a;
    init(a,overflow_safe);

    integer_type a_orig(a);

    integer_type b = a+a+a;
    integer_type c = a * integer_type(3);

    BOOST_CHECK_EQUAL(c,b);

    //Check that a hasn't changed
    BOOST_CHECK_EQUAL(a,a_orig);
}

VLI_STATIC_TEST( multiplies_carry_int )
{
    integer_type a,c;
    long int b = 2;

    for(std::size_t i(0); i < integer_type::numwords-1; ++i)
        a[i] = std::numeric_limits<integer_type::value_type>::max();

    a*=b;

    c[0] = std::numeric_limits<integer_type::value_type>::max()-1;

    for(std::size_t i(1); i < integer_type::numwords-1; ++i)
        c[i] = std::numeric_limits<integer_type::value_type>::max();

    c[integer_type::numwords-1] = 1;

    BOOST_CHECK_EQUAL(a,c);
}

VLI_STATIC_TEST( multiplies_carry )
{
    integer_type a,c,b;

    b[0]=2;

    for(std::size_t i(0); i < integer_type::numwords-1; ++i)
        a[i] = std::numeric_limits<integer_type::value_type>::max();

    a*=b;

    c[0] = std::numeric_limits<integer_type::value_type>::max()-1;

    for(std::size_t i(1); i < integer_type::numwords-1; ++i)
        c[i] = std::numeric_limits<integer_type::value_type>::max();

    c[integer_type::numwords-1] = 1;

    BOOST_CHECK_EQUAL(a,c);
}
