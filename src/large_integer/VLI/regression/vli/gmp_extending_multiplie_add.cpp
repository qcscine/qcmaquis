#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_extending_multiplie_add )
{
    vli_type a,b;
    typename double_sized_vli<vli_type>::type c;

    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);
    init(c,multiplies_overflow_safe);

    mpz_class agmp(a), bgmp(b), cgmp(c);

    multiply_add_extend(c,a,b);
    cgmp += agmp * bgmp;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
    BOOST_CHECK_EQUAL(mpz_class(b),bgmp);
    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_extending_multiplie_add_negative_a )
{
    vli_type a,b;
    typename double_sized_vli<vli_type>::type c;

    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);
    init(c,multiplies_overflow_safe);

    negate_inplace(a);

    mpz_class agmp(a), bgmp(b), cgmp(c);

    multiply_add_extend(c,a,b);
    cgmp += agmp * bgmp;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
    BOOST_CHECK_EQUAL(mpz_class(b),bgmp);
    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_extending_multiplie_add_negative_b )
{
    vli_type a,b;
    typename double_sized_vli<vli_type>::type c;

    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);
    init(c,multiplies_overflow_safe);

    negate_inplace(b);

    mpz_class agmp(a), bgmp(b), cgmp(c);

    multiply_add_extend(c,a,b);
    cgmp += agmp * bgmp;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
    BOOST_CHECK_EQUAL(mpz_class(b),bgmp);
    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_extending_multiplie_add_negative_c )
{
    vli_type a,b;
    typename double_sized_vli<vli_type>::type c;

    init(a,multiplies_overflow_safe);
    init(b,multiplies_overflow_safe);
    init(c,multiplies_overflow_safe);

    negate_inplace(c);

    mpz_class agmp(a), bgmp(b), cgmp(c);

    multiply_add_extend(c,a,b);
    cgmp += agmp * bgmp;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
    BOOST_CHECK_EQUAL(mpz_class(b),bgmp);
    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
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

    mpz_class agmp(a), bgmp(b), cgmp(c);

    multiply_add_extend(c,a,b);
    cgmp += agmp * bgmp;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
    BOOST_CHECK_EQUAL(mpz_class(b),bgmp);
    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
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

    mpz_class agmp(a), bgmp(b), cgmp(c);

    multiply_add_extend(c,a,b);
    cgmp += agmp * bgmp;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
    BOOST_CHECK_EQUAL(mpz_class(b),bgmp);
    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

