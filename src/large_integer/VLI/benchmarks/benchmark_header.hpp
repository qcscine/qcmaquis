#ifndef VLI_BENCH_HEADER_HPP
#define VLI_BENCH_HEADER_HPP

#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/preprocessor/cat.hpp>
#include <stdexcept>
#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial.hpp"
#include "vli/polynomial/monomial.hpp"
#include "vli/vli_cpu.h"
#include "vli/vli_traits.hpp"
#include "utils/timings.h"
#include "regression/vli_test.hpp"

using vli::test::fill_random;
using vli::test::fill_poly_random;
using vli::test::fill_vector_random;
using vli::test::fill_vector_negate;

namespace vli {
    namespace bench {

     typedef vli::vli_cpu<unsigned long int, VLI_SIZE> vli_type;
     typedef vli::vli_cpu<unsigned long int, 2*VLI_SIZE> vli_type_double;
    
     typedef vli::polynomial< vli_type, ORDER_POLY > polynomial_type;
     typedef vli::polynomial< vli_type_double, 2*ORDER_POLY > polynomial_type_double;
    
     typedef vli::vector_polynomial<polynomial_type> vector_type;

    }
}

#endif // VLI_BENCH_HEADER_HPP
