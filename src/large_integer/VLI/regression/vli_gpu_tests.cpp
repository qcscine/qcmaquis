#define BOOST_TEST_MODULE vli_gpu
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

// Andreas as the lib is now templated we can only run one test
typedef boost::mpl::list<
        vli_gpu<vli::detail::type_vli::BaseInt, vli::detail::size_vli::value>
        > vli_types;


//TODO is the gpu_manager needed for the execution of the other vli_gpu tests?
//yeap I am now getting the maximum of threads per block per the manager
//and I always hope have the manager as master chief for the mix mode CPU/GPU
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
    fill_random(a);
    vli_cpu<typename Vli::value_type, Vli::size> a_orig(a);

    Vli ag(a);

    BOOST_CHECK_EQUAL((a == vli_cpu_from(ag)), true);
    BOOST_CHECK_EQUAL((vli_cpu_from(ag) == a), true);
    BOOST_CHECK_EQUAL(a,a_orig);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( proxy_access, Vli, vli_types )
{
    boost::uniform_int<typename Vli::value_type> rnd(0,max_int_value<Vli>::value);

    // Read access
    vli_cpu<typename Vli::value_type,Vli::size> a;
    fill_random(a);
    vli_cpu<typename Vli::value_type,Vli::size> a_orig(a);

    Vli ag(a);

    for(typename Vli::size_type i=0; i < Vli::size; ++i)
        BOOST_CHECK_EQUAL(ag[i],a[i]);
    
    // Write access
    Vli bg;
    for(typename Vli::size_type i=0; i < Vli::size; ++i)
        bg[i] = a[i];
    
    for(typename Vli::size_type i=0; i < Vli::size; ++i)
        BOOST_CHECK_EQUAL(bg[i],a[i]);
    
    BOOST_CHECK_EQUAL(a,a_orig);
}

// Load all the other tests
#include "regression/vli_number_common_tests.hpp"
