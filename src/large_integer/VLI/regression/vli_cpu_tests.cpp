#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>

#include "detail/vli_size_param.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "gmpxx.h"
#include "regression/vli_test.hpp"

using vli::vli_cpu;
using vli::max_int_value;

using vli::test::rnd_valid_int;
using vli::test::fill_random;


typedef vli::test::vli_cpu_type_list vli_types;

//
//  load all tests
//
#include "regression/vli_number_common_tests.hpp"

