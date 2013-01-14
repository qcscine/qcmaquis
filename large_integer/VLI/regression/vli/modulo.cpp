#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_STATIC_TEST( divide_by_zero )
{
    vli_type a;
    vli_type b(0);
    init(a);
    bool ok = false;
    try{
        a %= b;
    } catch (vli::vli_division_by_zero_error& e){
        ok = true;
    }
    BOOST_CHECK_EQUAL(ok, true);
}

VLI_FUZZABLE_TEST( modulo_identity )
{
    vli_type a;
    init(a);
    vli_type a_orig(a);
    vli_type b = a;
    b %= a;
    BOOST_CHECK_EQUAL(b, vli_type(0));
    BOOST_CHECK_EQUAL(a, a_orig);
}

VLI_FUZZABLE_TEST( modulo )
{
    vli_type a;
    init(a);

    vli_type b(a);

    b += 1;
    vli_type a_orig(a);
    vli_type b_orig(b);

    vli_type c = b % a;

    BOOST_CHECK_EQUAL(c, vli_type(1));
    BOOST_CHECK_EQUAL(b, b_orig);
    BOOST_CHECK_EQUAL(a, a_orig);
}

VLI_FUZZABLE_TEST( modulo_assign_modulo_equivalence )
{
    vli_type a;
    vli_type b;
    init(a);
    init(b);

    vli_type a_orig(a);
    vli_type b_orig(b);
    vli_type c(b);

    c %= a;
    vli_type d = b % a;

    BOOST_CHECK_EQUAL(c,d);
    BOOST_CHECK_EQUAL(a,a_orig);
    BOOST_CHECK_EQUAL(b,b_orig);
}
