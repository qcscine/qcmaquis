#include <gmpxx.h>
#include <regression/vli/test_header.hpp>

using namespace vli::test;

VLI_FUZZABLE_TEST( gmp_multiplies_positive_positive_int )
{
    integer_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);

    mpz_class agmp(a);

    integer_type c = a*b;
    mpz_class cgmp = agmp * b;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_multiplies_positive_negative_int )
{
    integer_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);

    b = -b;

    mpz_class agmp(a);

    integer_type c = a*b;
    mpz_class cgmp = agmp * b;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_multiplies_negative_positive_int )
{
    integer_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);

    negate_inplace(a);

    mpz_class agmp(a);

    integer_type c = a*b;
    mpz_class cgmp = agmp * b;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_multiplies_negative_negative_int )
{
    integer_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);

    b = -b;
    negate_inplace(a);

    mpz_class agmp(a);

    integer_type c = a*b;
    mpz_class cgmp = agmp * b;

    BOOST_CHECK_EQUAL(mpz_class(c),cgmp);
}

VLI_FUZZABLE_TEST( gmp_multiplies_assign_positive_positive_int )
{
    integer_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);

    mpz_class agmp(a);

    a*=b;
    agmp*=b;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
}


VLI_FUZZABLE_TEST( gmp_multiplies_assign_negative_positive_int )
{
    integer_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);

    negate_inplace(a);

    mpz_class agmp(a);

    a*=b;
    agmp*=b;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
}

VLI_FUZZABLE_TEST( gmp_multiplies_assign_positive_negative_int )
{
    integer_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);
    b = -b;

    mpz_class agmp(a);

    a*=b;
    agmp*=b;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
}

VLI_FUZZABLE_TEST( gmp_multiplies_assign_negative_negative_int )
{
    integer_type a;
    int b;
    init(a,overflow_safe); // TODO is this safe enough?
    init(b);
    negate_inplace(a);
    b = -b;

    mpz_class agmp(a);

    a*=b;
    agmp*=b;

    BOOST_CHECK_EQUAL(mpz_class(a),agmp);
}
