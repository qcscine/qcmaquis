#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include "gpu/GpuManager.h"
#include "gpu/GpuManager.hpp"
#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "vli_gpu/vli_number_gpu.hpp"
#include "gmpxx.h"

#include "regression/common_test_functions.hpp"

using vli::vli_cpu;
using vli::max_int_value;
using vli::vli_gpu;

using vli::test::rnd_valid_int;
using vli::test::fill_random;

typedef boost::mpl::list<
        vli_gpu<unsigned int,2>,
        vli_gpu<unsigned int,4>,
        vli_gpu<unsigned int,8>,
        vli_gpu<unsigned int,16>,
        vli_gpu<unsigned long int,2>,
        vli_gpu<unsigned long int,4>,
        vli_gpu<unsigned long int,8>,
        vli_gpu<unsigned long int,16>
        > vli_types;


//TODO is the gpu_manager needed for the execution of the other vli_gpu tests?
BOOST_AUTO_TEST_CASE(gpu_manager)
{
	gpu::gpu_manager* GPU;
	GPU->instance();
    //TODO correct type?
	double FreqGPU = GPU->instance().GetDeviceProperties().clockRate;
    std::cout<<" FreqGPU "<<FreqGPU<<" Hz"<<std::endl;
    GPU->instance().destructor();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( conversion_from_cpu_and_equal_test, Vli, vli_types )
{
    vli_cpu<typename Vli::value_type, Vli::size> a;
    Vli ag(a);

    BOOST_CHECK_EQUAL(a,ag);
    BOOST_CHECK_EQUAL(ag,a);

    vli_cpu<typename Vli::value_type, Vli::size> b;
    fill_random(b);

    Vli bg(b);

    BOOST_CHECK_EQUAL(b,bg);
    BOOST_CHECK_EQUAL(bg,b);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( proxy_access, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);

    // Read access
    vli_cpu<typename Vli::value_type,Vli::size> a;
    fill_random(a);

    Vli ag(a);

    for(typename Vli::size_type i=0; i < Vli::size; ++i)
        BOOST_CHECK_EQUAL(ag[i],a[i]);
    
    // Write access
    Vli bg;
    for(typename Vli::size_type i=0; i < Vli::size; ++i)
        bg[i] = a[i];
    
    for(typename Vli::size_type i=0; i < Vli::size; ++i)
        BOOST_CHECK_EQUAL(bg[i],a[i]);
}

// Load all the other tests
#include "regression/vli_number_common_tests.hpp"
