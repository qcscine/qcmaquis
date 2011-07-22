
#define BOOST_TEST_MODULE vli_cpu


#include <boost/test/unit_test.hpp>
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "gmp.h"

using vli::vli_cpu;
using vli::max_int_value;
#define SIZE 8
#define TYPE int


BOOST_AUTO_TEST_CASE( constructors_test )
{
    vli_cpu<TYPE,SIZE> a;
    vli_cpu<TYPE,SIZE> b(0);

    BOOST_CHECK_EQUAL(a,b);

}

BOOST_AUTO_TEST_CASE( equal_operator )
{
    vli_cpu<TYPE,SIZE> a(0);
    vli_cpu<TYPE,SIZE> b;

    for(unsigned TYPE i=0; i < SIZE; ++i)
    {
        b[i] = 1;
        BOOST_CHECK_EQUAL(false,(a == b));
        b[i] = 0;
    }

    BOOST_CHECK_EQUAL(a,b);

}

BOOST_AUTO_TEST_CASE( copy_constructor_and_assignment )
{ 
    vli_cpu<TYPE,SIZE> a;
    a[0]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[1]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[2]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[3]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[4]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[5]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[6]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);

    
    vli_cpu<TYPE,SIZE> b(a);
    BOOST_CHECK_EQUAL(a,b);

    vli_cpu<TYPE,SIZE> c(static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value));
    c = b;
    BOOST_CHECK_EQUAL(c,b);
}

BOOST_AUTO_TEST_CASE( plus )
{
    vli_cpu<TYPE,SIZE> a;
    a[0]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[1]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[2]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[3]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[4]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[5]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[6]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[7]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    
    vli_cpu<TYPE,SIZE> b;
    b[0]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[1]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[2]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[3]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[4]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[5]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[7]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    
    vli_cpu<TYPE,SIZE> ab = a + b;

    b += a;

    BOOST_CHECK_EQUAL(b,ab);
}



BOOST_AUTO_TEST_CASE( multiplies )
{
    vli_cpu<TYPE,SIZE> a;
    a[0]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[1]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[2]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[3]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[4]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[5]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[6]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[7]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    
    vli_cpu<TYPE,SIZE> b = a+a+a;
    vli_cpu<TYPE,SIZE> c = a * vli_cpu<TYPE,SIZE>(3);

    BOOST_CHECK_EQUAL(b,c);
}

BOOST_AUTO_TEST_CASE( addition_gmp )
{
    vli_cpu<TYPE,SIZE> a;
    a[0]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[1]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[2]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[3]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[4]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[5]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[6]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[7]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    
    vli_cpu<TYPE,SIZE> b;
    b[0]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[1]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[2]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[3]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[4]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[5]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[7]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    
    vli_cpu<TYPE,SIZE> c;

    mpz_t agmp, bgmp, cgmp;
    
    mpz_init_set_str (agmp, a.get_char(), 10);
    mpz_init_set_str (bgmp, b.get_char(), 10);
    mpz_init_set_str (cgmp, c.get_char(), 10);

    
    c = a+b;
    mpz_add (cgmp, bgmp, agmp);	
    
    char str[128];
    
    for(int i=0;i<128; i++)
        str[i]=0;
    
    mpz_get_str(str,10,bgmp);
    
    std::string strname(str);
    
    BOOST_CHECK_EQUAL(strname,c.get_str());
    
}



BOOST_AUTO_TEST_CASE( multi_gmp )
{
    vli_cpu<TYPE,SIZE> a;
    a[0]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[1]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[2]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    a[3]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    
    vli_cpu<TYPE,SIZE> b;
    b[0]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    b[1]=static_cast<TYPE>(drand48())%(max_int_value<vli_cpu<TYPE,SIZE> >::value);
    
    
    vli_cpu<TYPE,SIZE> c;
    
    mpz_t agmp, bgmp, cgmp;
    
    mpz_init_set_str (agmp, a.get_char(), 10);
    mpz_init_set_str (bgmp, b.get_char(), 10);
    mpz_init_set_str (cgmp, c.get_char(), 10);
    
    
    c = a*b;
    mpz_mul (cgmp, bgmp, agmp);	
    
    char str[128];
    
    for(int i=0;i<128; i++)
        str[i]=0;
    
    mpz_get_str(str,10,bgmp);
    
    std::string strname(str);
    
    BOOST_CHECK_EQUAL(strname,c.get_str());
    
}

BOOST_AUTO_TEST_CASE( comparison_vli )
{
    vli_cpu<TYPE,SIZE> a(0);
    vli_cpu<TYPE,SIZE> b(0);
    
    BOOST_CHECK_EQUAL(a<b, false);
    BOOST_CHECK_EQUAL(a>b, false);
    
    a[1] = 1;
    b[1] = 1;

    BOOST_CHECK_EQUAL(a<b, false);
    BOOST_CHECK_EQUAL(a>b, false);

    b[2] = 1;
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

BOOST_AUTO_TEST_CASE( comparison_int )
{
    vli_cpu<TYPE,SIZE> a(0);

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
    


    vli_cpu<TYPE,SIZE> b(0);
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
