#include <regression/vli_cpu/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( plus_assign_plus_equivalence )
{
    vli_type a;
    vli_type b;
    
    init(a);
    init(b);

    vli_type b_orig(b);
   
    vli_type ab = a + b;
    vli_type ba = b + a;
    a += b;
    BOOST_CHECK_EQUAL(a,ab);
    BOOST_CHECK_EQUAL(a,ba);

    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_FUZZABLE_TEST( plus_assign_plus_equivalence_int )
{
    vli_type a;
    init(a);
    int b;
    init(b);

    int b_orig(b);

    vli_type ab = a + b;
    vli_type ba = b + a;
    a += b;
    BOOST_CHECK_EQUAL(a,ab);
    BOOST_CHECK_EQUAL(a,ba);
    
    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_STATIC_TEST( plus_assign_carry )
{
    vli_type a,b,c;
    for(std::size_t i(0); i < vli_type::size-1; ++i) 
        a[i] = static_cast<vli_type::value_type>(-1);  
    b[0] = 1;
    a+=b; 
    c[vli_type::size-1] = 1;  
    BOOST_CHECK_EQUAL(a,c);
}

VLI_STATIC_TEST( plus_assign_overflow )
{
    // TODO this test could also be a fuzz test
    vli_type a;
    for(std::size_t i(0); i < vli_type::size; ++i) 
        a[i] = static_cast<vli_type::value_type>(-1);  
    
    vli_type c;

    a+=1;

    std::cout << std::hex << a << std::endl;
    std::cout << std::hex <<  c << std::endl;
    
    BOOST_CHECK_EQUAL(a,c);
} 

VLI_STATIC_TEST( plus_assign_overflow_bug )
{
    // TODO this test could also be a fuzz test
    vli_type a;
    vli_type b(1);
    init(a,max_positive);
    vli_type c(a);
    negate_inplace(c);
    
    a+=b;
    
    BOOST_CHECK_EQUAL(a,c);
}   

