#include <regression/vli_cpu/test_header.hpp>

using namespace vli::test;

VLI_STATIC_TEST( pointer_range_overflow )
{
    vli_type a;
    init(a);
    vli_type b(a);
    vli_type a_orig(a);
    
    b -= vli_type(1);

    BOOST_CHECK_EQUAL(b,a-vli_type(1));
    

    vli_type *c = new vli_type[3];
    c[0] = vli_type(0);
    c[1] = a;
    c[2] = vli_type(0);

    c[1] *= b;

    a *= b;

    BOOST_CHECK_EQUAL(c[0],vli_type(0));
    BOOST_CHECK_EQUAL(c[1],a);
    BOOST_CHECK_EQUAL(c[2],vli_type(0));

    delete[] c;
}

VLI_STATIC_TEST( multiplies_by_two_not_equal_minus_assign_one )
{
    vli_type a;
    for(std::size_t i=0; i<vli_type::size; ++i)
        a[i] = std::numeric_limits<vli_type::value_type>::max();

    vli_type c(a); 
    vli_type b(2);
    
    a *= b;
    c -= 1;
    
    BOOST_CHECK_EQUAL((a == c), true);
}
