#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_STATIC_TEST( constructors_test )
{
    integer_type a;
    integer_type b(0);

    BOOST_CHECK_EQUAL(a,b);
}

VLI_FUZZABLE_TEST( copy_constructor_and_assignment )
{ 
    integer_type a;
    init(a);
    integer_type b(a);
    BOOST_CHECK_EQUAL(a,b);

    integer_type c;
    init(c);

    c = b;
    BOOST_CHECK_EQUAL(c,b);
    
    // Check if c stays the same if we change b
    b[1] = b[1]+1;
    BOOST_CHECK_EQUAL(c == b, false);
}

