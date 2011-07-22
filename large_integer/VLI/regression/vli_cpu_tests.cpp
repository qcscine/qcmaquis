#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "gmpxx.h"

using vli::vli_cpu;
using vli::max_int_value;
#define SIZE 8

typedef boost::mpl::list<unsigned int, unsigned long int> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors_test, T, test_types )
{
    vli_cpu<T,SIZE> a;
    vli_cpu<T,SIZE> b(0);

    BOOST_CHECK_EQUAL(a,b);

}

BOOST_AUTO_TEST_CASE_TEMPLATE( equal_operator, T, test_types )
{
    vli_cpu<T,SIZE> a(0);
    vli_cpu<T,SIZE> b;

    for(typename vli_cpu<T,SIZE>::size_type i=0; i < SIZE; ++i)
    {
        b[i] = 1;
        BOOST_CHECK_EQUAL(false,(a == b));
        b[i] = 0;
    }

    BOOST_CHECK_EQUAL(a,b);

}

BOOST_AUTO_TEST_CASE_TEMPLATE( copy_constructor_and_assignment, T, test_types )
{ 
    vli_cpu<T,SIZE> a;
    a[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);

    
    vli_cpu<T,SIZE> b(a);
    BOOST_CHECK_EQUAL(a,b);

    vli_cpu<T,SIZE> c(static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value));
    c = b;
    BOOST_CHECK_EQUAL(c,b);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus, T, test_types )
{
    vli_cpu<T,SIZE> a;
    a[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[7]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    
    vli_cpu<T,SIZE> b;
    b[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[7]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    
    vli_cpu<T,SIZE> ab = a + b;

    b += a;

    BOOST_CHECK_EQUAL(b,ab);
}



BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies, T, test_types )
{
    vli_cpu<T,SIZE> a;
    a[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[7]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    
    vli_cpu<T,SIZE> b = a+a+a;
    vli_cpu<T,SIZE> c = a * vli_cpu<T,SIZE>(3);

    BOOST_CHECK_EQUAL(b,c);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( addition_gmp, T, test_types )
{
    vli_cpu<T,SIZE> a;
    a[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[6]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[7]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    
    vli_cpu<T,SIZE> b;
    b[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[4]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[5]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[7]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    
    vli_cpu<T,SIZE> c;

    mpz_class agmp(a.get_str()), bgmp(b.get_str()), cgmp(c.get_str());
    
    c = a+b;
    cgmp = agmp + bgmp;
    
    BOOST_CHECK_EQUAL(cgmp.get_str(),c.get_str());
}



BOOST_AUTO_TEST_CASE_TEMPLATE( multi_gmp, T, test_types )
{
    vli_cpu<T,SIZE> a;
    a[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[2]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    a[3]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    
    vli_cpu<T,SIZE> b;
    b[0]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    b[1]=static_cast<T>(drand48())%(max_int_value<vli_cpu<T,SIZE> >::value);
    
    
    vli_cpu<T,SIZE> c;
    
    mpz_class agmp(a.get_str()), bgmp(b.get_str()), cgmp(c.get_str());
    
    c = a*b;
    cgmp = agmp * bgmp;
    
    BOOST_CHECK_EQUAL(cgmp.get_str(),c.get_str());
}

BOOST_AUTO_TEST_CASE_TEMPLATE( comparison_vli, T, test_types )
{
    vli_cpu<T,SIZE> a(0);
    vli_cpu<T,SIZE> b(0);
    
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

BOOST_AUTO_TEST_CASE_TEMPLATE( comparison_int, T, test_types )
{
    vli_cpu<T,SIZE> a(0);

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
    


    vli_cpu<T,SIZE> b(0);
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
