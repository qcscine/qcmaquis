#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_plus_inline )
{
    integer_type a,b;

    init(a,overflow_safe);
    init(b,overflow_safe);

    mpz_class agmp(a), bgmp(b);

    vli::detail::helper_inline_sub<integer_type::numwords>::inline_sub(&a[0],&b[0]);
    agmp -= bgmp;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
}

VLI_FUZZABLE_TEST( gmp_plus_inline_constant )
{
    static boost::random::mt19937                                           rng;
    static boost::random::uniform_int_distribution<integer_type::value_type>    integer_value_type_max_rnd;

    integer_type a;
    
    init(a,overflow_safe);

    long  int constant =  integer_value_type_max_rnd(rng);

    mpz_class agmp(a), bgmp(constant);
    
    
    
    vli::detail::helper_inline_sub<integer_type::numwords>::inline_sub(&a[0],constant);
    agmp -= bgmp;
    
    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
}
/*
VLI_STATIC_TEST( plus_assign_inline_carry )
{
    integer_type a,b,c;
    for(std::size_t i(0); i < integer_type::numwords-1; ++i)
        a[i] = 0;
    b[0] = 1;
    vli::detail::helper_inline_sub<integer_type::numwords>::inline_sub(&a[0],&b[0]);
    c[integer_type::numwords-1] = -1;
    BOOST_CHECK_EQUAL(a,c);
}
*/
