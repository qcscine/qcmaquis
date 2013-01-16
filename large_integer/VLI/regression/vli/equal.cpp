#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_STATIC_TEST( equal_operator )
{
    integer_type a(0);
    integer_type b;

    for(integer_type::size_type i=0; i < integer_type::numwords; ++i)
    {
        b[i] = 1;
        BOOST_CHECK_EQUAL((a == b),false);
        b[i] = 0;
    }

    BOOST_CHECK_EQUAL(a,b);
}
