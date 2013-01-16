#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_STATIC_TEST( pointer_range_overflow )
{
    integer_type a;
    init(a);
    integer_type b(a);
    integer_type a_orig(a);

    b -= integer_type(1);

    BOOST_CHECK_EQUAL(b,a-integer_type(1));


    integer_type *c = new integer_type[3];
    c[0] = integer_type(0);
    c[1] = a;
    c[2] = integer_type(0);

    c[1] *= b;

    a *= b;

    BOOST_CHECK_EQUAL(c[0],integer_type(0));
    BOOST_CHECK_EQUAL(c[1],a);
    BOOST_CHECK_EQUAL(c[2],integer_type(0));

    delete[] c;
}

VLI_STATIC_TEST( multiplies_by_two_not_equal_minus_assign_one )
{
    integer_type a;
    for(std::size_t i=0; i<integer_type::numwords; ++i)
        a[i] = std::numeric_limits<integer_type::value_type>::max();

    integer_type c(a);
    integer_type b(2);

    a *= b;
    c -= 1;

    BOOST_CHECK_EQUAL((a == c), true);
}
