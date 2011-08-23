#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "polynomial/monomial.hpp"
#include "gmpxx.h"

#include "regression/common_test_functions.hpp"

using vli::vli_cpu;
using vli::max_int_value;
using vli::vli_gpu;
using vli::monomial;

using vli::test::fill_random;

typedef boost::mpl::list<
        vli_cpu<unsigned int,8>,
        vli_cpu<unsigned int,16>,
        vli_cpu<unsigned long int,8>,
        vli_cpu<unsigned long int,16>,        
        vli_gpu<unsigned int,8>,
        vli_gpu<unsigned int,16>,
        vli_gpu<unsigned long int,8>,
        vli_gpu<unsigned long int,16> 
        > vli_types;

#include "regression/vli_monomial_commun_tests.hpp" 