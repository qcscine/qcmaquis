#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "detail/vli_size_param.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "gmpxx.h"

#include "regression/common_test_functions.hpp"

using vli::vli_cpu;
using vli::max_int_value;

using vli::test::rnd_valid_int;
using vli::test::fill_random;

// Andreas as the lib is now templated we can only run one test only true for the gpu
typedef boost::mpl::list<
        vli_cpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value>
        > vli_types;

/*
        vli_cpu<unsigned int,2>,
        vli_cpu<unsigned int,4>,
        vli_cpu<unsigned int,8>,
        vli_cpu<unsigned int,16>, 

        vli_cpu<unsigned long int,2>,
        vli_cpu<unsigned long int,4>, 
        vli_cpu<unsigned long int,8>,
        vli_cpu<unsigned long int,16>,  */
       

/**
    load all over tests
*/
#include "regression/vli_number_common_tests.hpp"

