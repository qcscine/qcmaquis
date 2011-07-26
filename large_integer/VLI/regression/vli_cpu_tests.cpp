#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>


#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "gmpxx.h"

#include "regression/common_test_functions.hpp"

using vli::vli_cpu;
using vli::max_int_value;

using vli::test::rnd_valid_int;
using vli::test::fill_random;

typedef boost::mpl::list<
        vli_cpu<unsigned int,2>,
        vli_cpu<unsigned int,4>,
        vli_cpu<unsigned int,8>,
        vli_cpu<unsigned int,16>,
        vli_cpu<unsigned long int,2>,
        vli_cpu<unsigned long int,4>,
        vli_cpu<unsigned long int,8>,
        vli_cpu<unsigned long int,16>
        > vli_types;

#include "regression/vli_number_common_tests.hpp"
