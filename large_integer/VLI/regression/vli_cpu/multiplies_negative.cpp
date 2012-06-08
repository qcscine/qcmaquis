#include <regression/vli_cpu/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( multiplies_assign_multiplies_equivalence_negative )
{
    vli_type a;
    vli_type b;
    init(a);
    init(b);
    vli_type b_orig(b);
    
    negate_inplace(a); 
    
    vli_type ab = a*b;
    vli_type ba = b*a;
    a*=b;
    
    BOOST_CHECK_EQUAL(a,ab);
    BOOST_CHECK_EQUAL(a,ba);
    
    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}


VLI_FUZZABLE_TEST( multiplies_assign_multiplies_equivalence_int_negative )
{
    vli_type a;
    init(a);
    int b;
    init(b);
    
    b = -b;
    
    int b_orig(b);
    
    vli_type ab = a*b;
    vli_type ba = b*a;
    a*=b;
    
    BOOST_CHECK_EQUAL(a,ab);
    BOOST_CHECK_EQUAL(a,ba);
    
    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_STATIC_TEST( multiplies_plus_relation_negative )
{
    vli_type a;
    init(a,overflow_safe);
    
    vli_type a_orig(a);
    
    vli_type b = -a-a-a;
    vli_type c = a * -vli_type(3); 
    
    BOOST_CHECK_EQUAL(c,b);
    
    //Check that a hasn't changed
    BOOST_CHECK_EQUAL(a,a_orig);
}

VLI_STATIC_TEST( multiplies_plus_relation_int_negative )
{
    vli_type a;
    init(a,overflow_safe); 
    vli_type a_orig(a);
    
    vli_type b = -a-a-a;
    vli_type c = a * -3;
    
    BOOST_CHECK_EQUAL(c,b);
    
    //Check that a hasn't changed
    BOOST_CHECK_EQUAL(a,a_orig);
}

VLI_STATIC_TEST( multiplies_carry_int_negative )
{
    vli_type a,c;
    long int b = -2;
    long int d = -1;

    for(std::size_t i(0); i < vli_type::size-1; ++i)
        a[i] = std::numeric_limits<vli_type::value_type>::max();

    a*=b;
    
    c[0] = 0x2;
 
    for(std::size_t i(1); i < vli_type::size; ++i) 
        c[i] = 0;
     
    c[vli_type::size-1] = std::numeric_limits<vli_type::value_type>::max()-1;
    BOOST_CHECK_EQUAL(a,c);

    a*=d;

    c[0] = std::numeric_limits<vli_type::value_type>::max()-1;
 
    for(std::size_t i(1); i < vli_type::size-1; ++i)
        c[i] = std::numeric_limits<vli_type::value_type>::max();
     
    c[vli_type::size-1] = 1;

    BOOST_CHECK_EQUAL(a,c);
}

VLI_STATIC_TEST( multiplies_carry_negative )
{
    vli_type a,c;
    vli_type b(2);
    vli_type d(1);
    negate_inplace(b);
    negate_inplace(d);

    for(std::size_t i(0); i < vli_type::size-1; ++i) 
        a[i] = std::numeric_limits<vli_type::value_type>::max();

    a*=b; 
    
    c[0] = 0x2;
 
    for(std::size_t i(1); i < vli_type::size; ++i) 
        c[i] = 0;  
     
    c[vli_type::size-1] = std::numeric_limits<vli_type::value_type>::max()-1;
    BOOST_CHECK_EQUAL(a,c);

    a*=d;

    c[0] = std::numeric_limits<vli_type::value_type>::max()-1;
 
    for(std::size_t i(1); i < vli_type::size-1; ++i) 
        c[i] = std::numeric_limits<vli_type::value_type>::max();  
     
    c[vli_type::size-1] = 1;

    BOOST_CHECK_EQUAL(a,c);
}
