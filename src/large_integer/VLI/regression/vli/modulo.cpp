#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_STATIC_TEST( divide_by_zero )
{
    integer_type a;
    integer_type b(0);
    init(a);
    bool ok = false;
    try{
        a %= b;
    } catch (vli::integer_division_by_zero_error& e){
        ok = true;
    }
    BOOST_CHECK_EQUAL(ok, true);
}

VLI_FUZZABLE_TEST( modulo_identity )
{
    integer_type a;
    init(a);
    integer_type a_orig(a);
    integer_type b = a;
    b %= a;
    BOOST_CHECK_EQUAL(b, integer_type(0));
    BOOST_CHECK_EQUAL(a, a_orig);
}

VLI_FUZZABLE_TEST( modulo )
{
    integer_type a;
    init(a);

    integer_type b(a);

    b += 1;
    integer_type a_orig(a);
    integer_type b_orig(b);

    integer_type c = b % a;

    BOOST_CHECK_EQUAL(c, integer_type(1));
    BOOST_CHECK_EQUAL(b, b_orig);
    BOOST_CHECK_EQUAL(a, a_orig);
}

VLI_FUZZABLE_TEST( modulo_assign_modulo_equivalence )
{
    integer_type a;
    integer_type b;
    init(a);
    init(b);

    integer_type a_orig(a);
    integer_type b_orig(b);
    integer_type c(b);

    c %= a;
    integer_type d = b % a;

    BOOST_CHECK_EQUAL(c,d);
    BOOST_CHECK_EQUAL(a,a_orig);
    BOOST_CHECK_EQUAL(b,b_orig);
}
