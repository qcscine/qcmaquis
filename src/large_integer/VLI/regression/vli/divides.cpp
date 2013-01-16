#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_STATIC_TEST( divide_by_zero )
{
    integer_type a;
    integer_type b(0);
    init(a);
    bool ok = false;

    try{
        a /= b;
    } catch (vli::integer_division_by_zero_error& e){
        ok = true;
    }
    BOOST_CHECK_EQUAL(ok, true);

}

VLI_FUZZABLE_TEST( multiply_divide )
{
    integer_type a;
    integer_type b;
    init(a);
    init(b);

    integer_type c(a);

    a /= b;

    a *= b;

    BOOST_CHECK_EQUAL( a <= c, true);
    BOOST_CHECK_EQUAL( a > c-b, true);
}


VLI_FUZZABLE_TEST( divide_assign_divide_equivalence )
{
    integer_type a;
    integer_type b;

    init(a);
    init(b);

    integer_type b_orig(b);

    integer_type ab = a / b;
    a /= b;
    BOOST_CHECK_EQUAL(a,ab);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}
