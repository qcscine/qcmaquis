
#include <regression/vli_cpu/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( minus_assign_minus_equivalence )
{
    vli_type a;
    vli_type b;
    init(a,max_positive);
    init(b,overflow_safe);
    vli_type b_orig(b);
    vli_type ab = a - b;
    vli_type ba = b - a;
    a -= b;
    BOOST_CHECK_EQUAL(a,ab);
    negate_inplace(a);
    BOOST_CHECK_EQUAL(a,ba);
    
    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_FUZZABLE_TEST( minus_assign_minus_equivalence_int )
{
    vli_type a;
    init(a);
    
    int b;
    init(b);
    int b_orig(b);
   
    vli_type ab = a - b;
    a -= b;
    BOOST_CHECK_EQUAL(a,ab);
    
    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_STATIC_TEST( minus_assign_borrow )
{
    vli_type a,b,c;

    a[vli_type::size-1] = 1;  
    b[0] = 1; 
    a-=b; 
    
    for(std::size_t i(0); i < vli_type::size-1; ++i) 
        c[i] = static_cast<vli_type::value_type>(-1);  
   
    BOOST_CHECK_EQUAL(a,c);
}
