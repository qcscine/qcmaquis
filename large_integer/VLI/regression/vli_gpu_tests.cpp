#define BOOST_TEST_MODULE vli_gpu
#include <boost/test/unit_test.hpp>

#include "vli/vli_cpu.h"
#include "vli/vli_traits.hpp"

#include "use_gmp_integers.hpp"
#include "minimal_polynomial.hpp"

#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial_cpu.h"
#include "vli/polynomial/monomial.hpp"

#include "vli/detail/bit_masks.hpp"


BOOST_AUTO_TEST_CASE(AdditionGPU)
{

    vli::vli_cpu< unsigned long int, 6> a,b(1);
    vli::vli_cpu< unsigned long int, 6> c(0);
    vli::vli_cpu< unsigned long int, 6> d(0);
   
    a[0] = 0xFFFFFFFFFFFFFFFF; 
    a[1] = 0xFFFFFFFFFFFFFFFF; 
    a[2] = 0xFFFFFFFFFFFFFFFF; 
    a[3] = 0xFFFFFFFFFFFFFFFF; 
    a[4] = 0xFFFFFFFFFFFFFFFF; 
    a[5] = 0xFFFFFFFFFFFFFFFF; 

    c = vli::detail::addition_gpu(a,b);

    BOOST_CHECK_EQUAL(c,d);

    a[5] = 0x0; 

    c = vli::detail::addition_gpu(a,b);

    d[5] = 0x1; 

    BOOST_CHECK_EQUAL(c,d);
}

BOOST_AUTO_TEST_CASE(MultiplicationGPU)
{
    typedef vli::vli_cpu<unsigned long int,3> vli_in;
    typedef vli::vli_cpu<unsigned long int,6> vli_out;
    typedef mpz_class large_int;

    large_int x,y,z; 

    vli_in c,b;
    vli_out a(0);
    vli_out d(0);
    vli_out e(0);

    c[0] = 0xFFFFFFFFFFFFFFFF;
    c[1] = 0xFFFFFFFFFFFFFFFF;
    c[2] = 0x0FFFFFFFFFFFFFFF;

    b[0] = 0xFFFFFFFFFFFFFFFF;
    b[1] = 0xFFFFFFFFFFFFFFFF;
    b[2] = 0x0FFFFFFFFFFFFFFF;
    
    x = b.get_str();
    y = c.get_str();
    z = x * y; 
    vli::mul(e,b,c);
    d = vli::detail::multiplication_gpu(a,b,c);
   
    BOOST_CHECK_EQUAL(d.get_str(),e.get_str());
    BOOST_CHECK_EQUAL(z.get_str(),d.get_str());
}
