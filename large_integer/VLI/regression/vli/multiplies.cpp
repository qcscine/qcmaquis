#include <regression/vli/test_header.hpp>

using namespace vlilib::test;


VLI_FUZZABLE_TEST( multiplies_assign_multiplies_equivalence )
{
    vli_type a;
    vli_type b;
    init(a);
    init(b);
    vli_type b_orig(b);

    vli_type ab = a*b;
    vli_type ba = b*a;
    a*=b;

    BOOST_CHECK_EQUAL(a,ab);
    BOOST_CHECK_EQUAL(a,ba);
    
    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_FUZZABLE_TEST( multiplies_assign_multiplies_equivalence_int )
{
    vli_type a;
    int b;
    init(a);
    init(b);
    
    int b_orig(b);

    vli_type ab = a*b;
    vli_type ba = b*a;
    a*=b;

    BOOST_CHECK_EQUAL(a,ab);
    BOOST_CHECK_EQUAL(a,ba);
    
    //Check that b hasn't changed
    BOOST_CHECK_EQUAL(b,b_orig);
}

VLI_STATIC_TEST( multiplies_plus_relation_int )
{
    vli_type a;
    init(a,overflow_safe);
    vli_type a_orig(a);
    
    vli_type b = a+a+a;
    vli_type c = a * 3;

    BOOST_CHECK_EQUAL(c,b);
    
    //Check that a hasn't changed
    BOOST_CHECK_EQUAL(a,a_orig);
}

VLI_STATIC_TEST( multiplies_plus_relation )
{
    vli_type a;
    init(a,overflow_safe);

    vli_type a_orig(a);
    
    vli_type b = a+a+a;
    vli_type c = a * vli_type(3); 

    BOOST_CHECK_EQUAL(c,b);
    
    //Check that a hasn't changed
    BOOST_CHECK_EQUAL(a,a_orig);
}

VLI_STATIC_TEST( multiplies_carry_int )
{
    vli_type a,c;
    long int b = 2;

    for(std::size_t i(0); i < vli_type::numwords-1; ++i) 
        a[i] = std::numeric_limits<vli_type::value_type>::max();

    a*=b; 
    
    c[0] = std::numeric_limits<vli_type::value_type>::max()-1;
 
    for(std::size_t i(1); i < vli_type::numwords-1; ++i) 
        c[i] = std::numeric_limits<vli_type::value_type>::max();
     
    c[vli_type::numwords-1] = 1;
   
    BOOST_CHECK_EQUAL(a,c);
}

VLI_STATIC_TEST( multiplies_carry )
{
    vli_type a,c,b;
      
    b[0]=2;    

    for(std::size_t i(0); i < vli_type::numwords-1; ++i) 
        a[i] = std::numeric_limits<vli_type::value_type>::max();

    a*=b; 

    c[0] = std::numeric_limits<vli_type::value_type>::max()-1;
 
    for(std::size_t i(1); i < vli_type::numwords-1; ++i) 
        c[i] = std::numeric_limits<vli_type::value_type>::max();
     
    c[vli_type::numwords-1] = 1;

    BOOST_CHECK_EQUAL(a,c);
}
