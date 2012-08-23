#include <regression/vli/test_header.hpp>

using namespace vlilib::test;

VLI_STATIC_TEST( equal_operator )
{
    vli_type a(0);
    vli_type b;

    for(vli_type::size_type i=0; i < vli_type::numwords; ++i)
    {
        b[i] = 1;
        BOOST_CHECK_EQUAL((a == b),false);
        b[i] = 0;
    }

    BOOST_CHECK_EQUAL(a,b);
}
