#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>

#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "gmpxx.h"

using vli::vli_cpu;
using vli::max_int_value;

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

boost::mt11213b rng;

#include "regression/vli_number_common_tests.hpp"
