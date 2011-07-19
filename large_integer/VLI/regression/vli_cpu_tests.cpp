
#define BOOST_TEST_MODULE vli_cpu


#include <boost/test/unit_test.hpp>


#include "vli_cpu/vli_number_cpu.hpp"

using vli::vli_cpu;


BOOST_AUTO_TEST_CASE( constructors_test )
{
    vli_cpu<int,8> a;
    vli_cpu<int,8> b(0);

    BOOST_CHECK_EQUAL(a,b);

}

BOOST_AUTO_TEST_CASE( equal_operator )
{
    vli_cpu<int,8> a(0);
    vli_cpu<int,8> b;

    for(unsigned int i=0; i < 8; ++i)
    {
        b[i] = 1;
        BOOST_CHECK_EQUAL(false,(a == b));
        b[i] = 0;
    }

    BOOST_CHECK_EQUAL(a,b);

}

BOOST_AUTO_TEST_CASE( copy_constructor_and_assignment )
{
    vli_cpu<int,8> a;
    a[0]=123;
    a[1]=198;
    a[2]=76;
    a[3]=22;
    a[4]=81;
    a[5]=249;
    a[6]=40;
    a[7]=107;

    vli_cpu<int,8> b(a);
    BOOST_CHECK_EQUAL(a,b);

    vli_cpu<int,8> c(5);
    c = b;
    BOOST_CHECK_EQUAL(c,b);
}

BOOST_AUTO_TEST_CASE( plus_assign )
{
// make a new one
//    BOOST_CHECK_EQUAL(d,cd);
}

BOOST_AUTO_TEST_CASE( plus )
{
    vli_cpu<int,8> a;
    a[0] = 252;
    a[1] = 245;
    a[2] = 44;
    a[3] = 97;
    a[4] = 106;
    a[5] = 219;
    a[6] = 198;
    a[7] = 104;

    vli_cpu<int,8> b;
    b[0] = 255;
    b[1] = 234;
    b[3] = 232;
    b[4] = 192;
    b[5] = 83;
    b[6] = 139;
    b[7] = 29;

    vli_cpu<int,8> ab = a + b;

    b += a;

    BOOST_CHECK_EQUAL(b,ab);
}



BOOST_AUTO_TEST_CASE( multiplies )
{
    vli_cpu<int,8> a;
    a[0] = 252;
    a[1] = 245;
    a[2] = 44;
    a[3] = 97;
    a[4] = 106;
    a[5] = 219;
    a[6] = 198;
    a[7] = 0;
    
    vli_cpu<int,8> b = a+a+a;
    vli_cpu<int,8> c = a * vli_cpu<int,8>(3);

    BOOST_CHECK_EQUAL(b,c);
}



