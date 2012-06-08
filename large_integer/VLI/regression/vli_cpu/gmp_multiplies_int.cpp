#include <regression/vli_cpu/test_header.hpp>

#include <gmpxx.h>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_multiplies_positive_positive_int )
{
    vli_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);
    
    mpz_class agmp(a.get_str());
    
    vli_type c = a*b;
    mpz_class cgmp = agmp * b;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_multiplies_positive_negative_int )
{
    vli_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);
    
    b = -b;
    
    mpz_class agmp(a.get_str());
    
    vli_type c = a*b;
    mpz_class cgmp = agmp * b;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_multiplies_negative_positive_int )
{
    vli_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);
    
    negate_inplace(a);
    
    mpz_class agmp(a.get_str());
    
    vli_type c = a*b;
    mpz_class cgmp = agmp * b;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_multiplies_negative_negative_int )
{
    vli_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);
    
    b = -b;
    negate_inplace(a);
    
    mpz_class agmp(a.get_str());
    
    vli_type c = a*b;
    mpz_class cgmp = agmp * b;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_multiplies_assign_positive_positive_int )
{
    vli_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);
    
    mpz_class agmp(a.get_str());
    
    a*=b;
    agmp*=b;
    
    BOOST_CHECK_EQUAL(a.get_str(),agmp.get_str());
}


VLI_FUZZABLE_TEST( gmp_multiplies_assign_negative_positive_int )
{
    vli_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);
    
    negate_inplace(a);
    
    mpz_class agmp(a.get_str());
    
    a*=b;
    agmp*=b;
    
    BOOST_CHECK_EQUAL(a.get_str(),agmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_multiplies_assign_positive_negative_int )
{
    vli_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);
    b = -b;
        
    mpz_class agmp(a.get_str());
    
    a*=b;
    agmp*=b;
    
    BOOST_CHECK_EQUAL(a.get_str(),agmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_multiplies_assign_negative_negative_int )
{
    vli_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);
    negate_inplace(a);
    b = -b;
    
    mpz_class agmp(a.get_str());
    
    a*=b;
    agmp*=b;
    
    BOOST_CHECK_EQUAL(a.get_str(),agmp.get_str());
}
