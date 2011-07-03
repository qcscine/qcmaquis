
#define BOOST_TEST_MODULE vli_cpu


#include <boost/test/unit_test.hpp>


#define SIZE_BITS 256
#include "vli_cpu/vli_number_cpu.hpp"


using vli::vli_cpu;


BOOST_AUTO_TEST_CASE( constructors_test )
{
    vli_cpu<int,8> a;
    vli_cpu<int,8> b(0);

    BOOST_CHECK_EQUAL(a,b);

}

BOOST_AUTO_TEST_CASE( copy_constructor_and_assignment )
{
    vli_cpu<int,8> a(10);
    vli_cpu<int,8> b(a);
    BOOST_CHECK_EQUAL(a,b);

    vli_cpu<int,8> c(5);
    c = b;
    BOOST_CHECK_EQUAL(c,b);
}
