#define BOOST_TEST_MODULE vli_cpu
#include <iostream>
#include <boost/test/unit_test.hpp>


#define SIZE_BITS 256
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_vector_cpu.hpp"


using vli::vli_cpu;
using vli::vli_vector_cpu;

BOOST_AUTO_TEST_CASE( constructors_test )
{
    vli_vector_cpu< vli_cpu<int> > a(10);
    vli_vector_cpu< vli_cpu<int> > b(a);
	
    BOOST_CHECK_EQUAL(a[0],b[0]);
}


BOOST_AUTO_TEST_CASE( copy_constructor_and_assignment )
{
    vli_vector_cpu<vli_cpu<int> > a(10);
    vli_vector_cpu<vli_cpu<int> > b(a);

    BOOST_CHECK_EQUAL(a[0],b[0]);
	
	vli_vector_cpu<vli_cpu<int> > c(5);
    c = b;
    BOOST_CHECK_EQUAL(c[0],b[0]);
}
