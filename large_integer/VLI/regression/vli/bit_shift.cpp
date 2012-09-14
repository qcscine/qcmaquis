#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_STATIC_TEST( left_shift )
{
    vli_type a;
    vli_type b;
    a[0] = 0xffffffffffffffff;
    b[0] = 0xfffffffffffffffe;
    b[vli_type::numwords-1] = 1;

    for(vli_type::size_type i=1; i < vli_type::numwords-1; ++i) {
        a[i] = 0xffffffffffffffff;
        b[i] = 0xffffffffffffffff;
    }

    a <<=1;

    BOOST_CHECK_EQUAL(a,b);
}

VLI_STATIC_TEST( right_shift )
{
    vli_type a;
    vli_type b;

    for(vli_type::size_type i=0; i < vli_type::numwords-1; ++i) {
        a[i] = 0xffffffffffffffff;
        b[i] = 0xffffffffffffffff;
    }

    b[vli_type::numwords-2] = 0x7fffffffffffffff;

    a >>=1;

    BOOST_CHECK_EQUAL(a,b);
}

VLI_STATIC_TEST( minus_1_right_shift )
{
    vli_type a(-1);

    for(vli_type::size_type i=0; i < 64; ++i) {
        a >>= 1;
        BOOST_CHECK_EQUAL(a,vli_type(-1));
    }
}

VLI_FUZZABLE_TEST(left_shift_mul )
{
    vli_type a;
    init(a,overflow_safe);
    vli_type b(a);

    a <<= 2;
    b *= 4;
    BOOST_CHECK_EQUAL(a,b);
}

VLI_FUZZABLE_TEST(left_right_shift)
{
    vli_type a;
    init(a,overflow_safe);
    vli_type b(a);

    a <<=13;
    a >>=13;

    BOOST_CHECK_EQUAL(a,b);
}
