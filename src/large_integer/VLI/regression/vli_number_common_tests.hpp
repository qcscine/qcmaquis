#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>

boost::mt11213b rng;

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors_test, Vli, vli_types )
{
    Vli a;
    Vli b(0);

    BOOST_CHECK_EQUAL(a,b);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( equal_operator, Vli, vli_types )
{
    Vli a(0);
    Vli b;

    for(typename Vli::size_type i=0; i < Vli::size; ++i)
    {
        b[i] = 1;
        BOOST_CHECK_EQUAL(false,(a == b));
        b[i] = 0;
    }

    BOOST_CHECK_EQUAL(a,b);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( copy_constructor_and_assignment, Vli, vli_types )
{ 
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);
    
    Vli a;
    
    for(typename Vli::size_type i=0; i<Vli::size; ++i)
            a[i] = rnd(rng);
    
    Vli b(a);
    BOOST_CHECK_EQUAL(a,b);

    Vli c;
    
    for(typename Vli::size_type i=0; i<Vli::size; ++i)
            c[i] = rnd(rng);

    c = b;
    BOOST_CHECK_EQUAL(c,b);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( negate, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);
    
    Vli a;
    
    for(typename Vli::size_type i=0; i<Vli::size; ++i)
            a[i] = rnd(rng);

    Vli b(a);
    BOOST_CHECK_EQUAL(a.is_negative(), false);
    
    a.negate();
    BOOST_CHECK_EQUAL(a.is_negative(), true);
    
    a.negate();
    BOOST_CHECK_EQUAL(a.is_negative(), false);
    BOOST_CHECK_EQUAL(a,b);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( negate_and_construct_from_negative_int, Vli, vli_types )
{
    Vli a(2437284);
    Vli am(-2437284);
    a.negate();
    BOOST_CHECK_EQUAL(a,am);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_assign_plus_equivalence, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);

    Vli a;
    Vli b;
    for(typename Vli::size_type i=0; i<Vli::size; ++i)
    {
        a[i] = rnd(rng);
        b[i] = rnd(rng);
    }
   
    Vli ab = a + b;
    b += a;
    BOOST_CHECK_EQUAL(b,ab);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_assign_minus_equivalence, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);

    Vli a;
    Vli b;
    for(typename Vli::size_type i=0; i<Vli::size; ++i)
    {
        a[i] = rnd(rng);
        b[i] = rnd(rng);
    }
   
    Vli ab = a - b;
    a -= b;
    BOOST_CHECK_EQUAL(a,ab);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_assign_multiplies_equivalence, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);

    Vli a;
    Vli b;
    for(typename Vli::size_type i=0; i<Vli::size; ++i)
    {
        a[i] = rnd(rng);
        b[i] = rnd(rng);
    }

    Vli c = a*b;
    a*=b;

    BOOST_CHECK_EQUAL(a,c);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);
    
    Vli a;
    
    for(typename Vli::size_type i=0; i<Vli::size-1; ++i)
        a[i] = rnd(rng);
    
    
    Vli b = a+a+a;
    Vli c = a * Vli(3);

    BOOST_CHECK_EQUAL(b,c);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_gmp, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);
    Vli a;
    Vli b;
    
    for(typename Vli::size_type i=0; i<Vli::size-1; ++i)
    {
        a[i] = rnd(rng);
        b[i] = rnd(rng);
    }
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    Vli c = a+b;
    mpz_class cgmp = agmp + bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_gmp, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);
    Vli a;
    Vli b;
    
    for(typename Vli::size_type i=0; i<Vli::size-1; ++i)
    {
        a[i] = rnd(rng);
        b[i] = rnd(rng);
    }
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    Vli c = a-b;
    mpz_class cgmp = agmp - bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_gmp, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);
    
    Vli a;
    Vli b;
    
    for(typename Vli::size_type i=0; i<Vli::size/2; ++i)
        a[i] = rnd(rng);
    
    for(typename Vli::size_type i=0; i<Vli::size/4; ++i)
        b[i] = rnd(rng);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    Vli c = a*b;
    mpz_class cgmp = agmp * bgmp;
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_negative_numbers_gmp, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);
    
    Vli a;
    Vli b;
    
    for(typename Vli::size_type i=0; i<Vli::size/2; ++i)
        a[i] = rnd(rng);
    a.negate();
    
    for(typename Vli::size_type i=0; i<Vli::size/4; ++i)
        b[i] = rnd(rng);
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str());
    
    Vli c = a*b;
    mpz_class cgmp = agmp * bgmp;
    
    b.negate();
    Vli d = a*b;
    mpz_class dgmp = agmp * (-bgmp);
    
    BOOST_CHECK_EQUAL(c.get_str(),cgmp.get_str());
    BOOST_CHECK_EQUAL(d.get_str(),dgmp.get_str());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( comparison_vli, Vli, vli_types )
{
    BOOST_STATIC_ASSERT(Vli::size > 1);

    Vli a(0);
    Vli b(0);
    
    BOOST_CHECK_EQUAL(a<b, false);
    BOOST_CHECK_EQUAL(a>b, false);
    
    a[1] = 1;
    b[1] = 1;

    BOOST_CHECK_EQUAL(a<b, false);
    BOOST_CHECK_EQUAL(a>b, false);

    b[1] = 2;
    a[0] = 1;

    BOOST_CHECK_EQUAL(a<b, true);
    BOOST_CHECK_EQUAL(b<a, false);
    BOOST_CHECK_EQUAL(a>b, false);
    BOOST_CHECK_EQUAL(b>a, true);

    // How about different signs?
    b.negate();

    BOOST_CHECK_EQUAL(a<b, false);
    BOOST_CHECK_EQUAL(b<a, true);
    BOOST_CHECK_EQUAL(a>b, true);
    BOOST_CHECK_EQUAL(b>a, false);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( comparison_int, Vli, vli_types )
{
    BOOST_STATIC_ASSERT(Vli::size > 1);
    
    Vli a(0);

    int zero = 0;
    int one = 1;
    int minus_one = -1;
    BOOST_CHECK_EQUAL(a<zero, false);
    BOOST_CHECK_EQUAL(a>zero, false);
    BOOST_CHECK_EQUAL(a<minus_one, false);
    BOOST_CHECK_EQUAL(a>minus_one, true);
    BOOST_CHECK_EQUAL(a<one, true);
    BOOST_CHECK_EQUAL(a>one, false);

    a+=1;
    BOOST_CHECK_EQUAL(a<zero, false);
    BOOST_CHECK_EQUAL(a>zero, true);
    BOOST_CHECK_EQUAL(a<minus_one, false);
    BOOST_CHECK_EQUAL(a>minus_one, true);
    BOOST_CHECK_EQUAL(a<one, false);
    BOOST_CHECK_EQUAL(a>one, false);

    a.negate();
    BOOST_CHECK_EQUAL(a<zero, true);
    BOOST_CHECK_EQUAL(a>zero, false);
    BOOST_CHECK_EQUAL(a<minus_one, false);
    BOOST_CHECK_EQUAL(a>minus_one, false);
    BOOST_CHECK_EQUAL(a<one, true);
    BOOST_CHECK_EQUAL(a>one, false);
    


    Vli b(0);
    b[1] = 1;
    BOOST_CHECK_EQUAL(b<zero, false);
    BOOST_CHECK_EQUAL(b>zero, true);
    BOOST_CHECK_EQUAL(b<minus_one, false);
    BOOST_CHECK_EQUAL(b>minus_one, true);
    BOOST_CHECK_EQUAL(b<one, false);
    BOOST_CHECK_EQUAL(b>one, true);

    b.negate();
    BOOST_CHECK_EQUAL(b<zero, true);
    BOOST_CHECK_EQUAL(b>zero, false);
    BOOST_CHECK_EQUAL(b<minus_one, true);
    BOOST_CHECK_EQUAL(b>minus_one, false);
    BOOST_CHECK_EQUAL(b<one, true);
    BOOST_CHECK_EQUAL(b>one, false);
}
