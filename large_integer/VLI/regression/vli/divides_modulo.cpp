#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( divide_modulo )
{
    integer_type a;
    integer_type b;
    init(a);
    init(b);

    integer_type c(a);

    integer_type mod = a % b;
    a /= b;

    BOOST_CHECK_EQUAL(a*b+mod,c);
}
