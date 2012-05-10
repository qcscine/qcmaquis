#include <regression/vli_cpu/test_header.hpp>
#include <gmpxx.h>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_extending_multiplie_add )
{
    vli_type a,b;
    typename double_sized_vli<vli_type>::type c;

    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);
    init(c,multiplies_overflow_safe);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str()), cgmp(c.get_str());
    
    muladd(c,a,b);
    cgmp += agmp * bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_extending_multiplie_add_negative_a )
{
    vli_type a,b;
    typename double_sized_vli<vli_type>::type c;
    
    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);
    init(c,multiplies_overflow_safe);
    
    negate_inplace(a);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str()), cgmp(c.get_str());
    
    muladd(c,a,b);
    cgmp += agmp * bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_extending_multiplie_add_negative_b )
{
    vli_type a,b;
    typename double_sized_vli<vli_type>::type c;
    
    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);
    init(c,multiplies_overflow_safe);
    
    negate_inplace(b);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str()), cgmp(c.get_str());
    
    muladd(c,a,b);
    cgmp += agmp * bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_extending_multiplie_add_negative_c )
{
    vli_type a,b;
    typename double_sized_vli<vli_type>::type c;
    
    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);
    init(c,multiplies_overflow_safe);
    
    negate_inplace(c);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str()), cgmp(c.get_str());
    
    muladd(c,a,b);
    cgmp += agmp * bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_extending_multiplie_add_negative_ab )
{
    vli_type a,b;
    typename double_sized_vli<vli_type>::type c;
    
    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);
    init(c,multiplies_overflow_safe);
    
    negate_inplace(a);
    negate_inplace(b);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str()), cgmp(c.get_str());
    
    muladd(c,a,b);
    cgmp += agmp * bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

VLI_FUZZABLE_TEST( gmp_extending_multiplie_add_negative_ac )
{
    vli_type a,b;
    typename double_sized_vli<vli_type>::type c;
    
    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);
    init(c,multiplies_overflow_safe);
    
    negate_inplace(a);
    negate_inplace(b);
    negate_inplace(c);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str()), cgmp(c.get_str());
    
    muladd(c,a,b);
    cgmp += agmp * bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

