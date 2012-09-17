#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( divide_modulo )
{
    vli_type a;
    vli_type b;
    init(a);
    init(b);

    vli_type c(a);

    vli_type mod = a % b;
    a /= b;

    BOOST_CHECK_EQUAL(a*b+mod,c);
}
