#include <regression/vli/test_header.hpp>

using namespace vlilib::test;

VLI_STATIC_TEST( comparison_vli )
{
    BOOST_STATIC_ASSERT(vli_type::numwords > 1);

    vli_type a(0);
    vli_type b(0);
    vli_type a_orig(a);
    vli_type b_orig(b);
    
    BOOST_CHECK_EQUAL(a<b, false);
    BOOST_CHECK_EQUAL(a>b, false);

    BOOST_CHECK_EQUAL(a,a_orig);
    BOOST_CHECK_EQUAL(b,b_orig);

    a[1] = 1; a_orig[1] = 1;
    b[1] = 1; b_orig[1] = 1;

    BOOST_CHECK_EQUAL(a<b, false);
    BOOST_CHECK_EQUAL(a>b, false);
    
    BOOST_CHECK_EQUAL(a,a_orig);
    BOOST_CHECK_EQUAL(b,b_orig);

    b[1] = 2; b_orig[1] = 2;
    a[0] = 1; a_orig[0] = 1;

    BOOST_CHECK_EQUAL(a<b, true);
    BOOST_CHECK_EQUAL(b<a, false);
    BOOST_CHECK_EQUAL(a>b, false);
    BOOST_CHECK_EQUAL(b>a, true);
    
    BOOST_CHECK_EQUAL(a,a_orig);
    BOOST_CHECK_EQUAL(b,b_orig);

    // How about different signs?
    negate_inplace(b); negate_inplace(b_orig);

    BOOST_CHECK_EQUAL(a<b, false);
    BOOST_CHECK_EQUAL(b<a, true);
    BOOST_CHECK_EQUAL(a>b, true);
    BOOST_CHECK_EQUAL(b>a, false);
    
    BOOST_CHECK_EQUAL(a,a_orig);
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_STATIC_TEST( comparison_int )
{
    BOOST_STATIC_ASSERT(vli_type::numwords > 1);
    
    vli_type a(0);
    vli_type a_orig(a);

    int zero = 0;
    int one = 1;
    int minus_one = -1;
    BOOST_CHECK_EQUAL(a<zero, false);
    BOOST_CHECK_EQUAL(a>zero, false);
    BOOST_CHECK_EQUAL(a<minus_one, false);
    BOOST_CHECK_EQUAL(a>minus_one, true);
    BOOST_CHECK_EQUAL(a<one, true);
    BOOST_CHECK_EQUAL(a>one, false);

    BOOST_CHECK_EQUAL(a,a_orig);
    BOOST_CHECK_EQUAL(zero, 0);
    BOOST_CHECK_EQUAL(one, 1);
    BOOST_CHECK_EQUAL(minus_one,-1);


    a += 1;  a_orig += 1;
    BOOST_CHECK_EQUAL(a<zero, false);
    BOOST_CHECK_EQUAL(a>zero, true);
    BOOST_CHECK_EQUAL(a<minus_one, false);
    BOOST_CHECK_EQUAL(a>minus_one, true);
    BOOST_CHECK_EQUAL(a<one, false);
    BOOST_CHECK_EQUAL(a>one, false);
    
    BOOST_CHECK_EQUAL(a,a_orig);

    negate_inplace(a); negate_inplace(a_orig);
    BOOST_CHECK_EQUAL(a<zero, true);
    BOOST_CHECK_EQUAL(a>zero, false);
    BOOST_CHECK_EQUAL(a<minus_one, false);
    BOOST_CHECK_EQUAL(a>minus_one, false);
    BOOST_CHECK_EQUAL(a<one, true);
    BOOST_CHECK_EQUAL(a>one, false);
    
    BOOST_CHECK_EQUAL(a,a_orig);
    

    vli_type b(0);
    b[1] = 1;
    vli_type b_orig(b);
    BOOST_CHECK_EQUAL(b<zero, false);
    BOOST_CHECK_EQUAL(b>zero, true);
    BOOST_CHECK_EQUAL(b<minus_one, false);
    BOOST_CHECK_EQUAL(b>minus_one, true);
    BOOST_CHECK_EQUAL(b<one, false);
    BOOST_CHECK_EQUAL(b>one, true);
    
    BOOST_CHECK_EQUAL(b,b_orig);

    negate_inplace(b); negate_inplace(b_orig);
    BOOST_CHECK_EQUAL(b<zero, true);
    BOOST_CHECK_EQUAL(b>zero, false);
    BOOST_CHECK_EQUAL(b<minus_one, true);
    BOOST_CHECK_EQUAL(b>minus_one, false);
    BOOST_CHECK_EQUAL(b<one, true);
    BOOST_CHECK_EQUAL(b>one, false);
    
    BOOST_CHECK_EQUAL(b,b_orig);
}
