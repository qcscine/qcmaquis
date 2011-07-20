
#define BOOST_TEST_MODULE vli_cpu


#include <boost/test/unit_test.hpp>
#include "vli_cpu/vli_number_cpu.hpp"
#include "gmp.h"

using vli::vli_cpu;
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
    a[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);

    
    vli_cpu<TYPE,SIZE> b(a);
    BOOST_CHECK_EQUAL(a,b);

    vli_cpu<TYPE,SIZE> c(static_cast<TYPE>(drand48())%(MAX_VALUE));
    c = b;
    BOOST_CHECK_EQUAL(c,b);
}

BOOST_AUTO_TEST_CASE( plus )
{
    vli_cpu<TYPE,SIZE> a;
    a[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[7]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    vli_cpu<TYPE,SIZE> b;
    b[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[7]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    vli_cpu<TYPE,SIZE> ab = a + b;

    b += a;

    BOOST_CHECK_EQUAL(b,ab);
}



BOOST_AUTO_TEST_CASE( multiplies )
{
    vli_cpu<TYPE,SIZE> a;
    a[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[7]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    vli_cpu<TYPE,SIZE> b = a+a+a;
    vli_cpu<TYPE,SIZE> c = a * vli_cpu<TYPE,SIZE>(3);

    BOOST_CHECK_EQUAL(b,c);
}

BOOST_AUTO_TEST_CASE( addition_gmp )
{
    vli_cpu<TYPE,SIZE> a;
    a[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[6]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[7]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    vli_cpu<TYPE,SIZE> b;
    b[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[4]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[5]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[7]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
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
    a[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[2]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    a[3]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    vli_cpu<TYPE,SIZE> b;
    b[0]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    b[1]=static_cast<TYPE>(drand48())%(MAX_VALUE);
    
    
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


