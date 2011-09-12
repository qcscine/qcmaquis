#define BOOST_TEST_MODULE vli_gpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/transform.hpp>
#include <gmpxx.h>

#include "vli/utils/gpu_manager.h"
#include "vli/utils/gpu_manager.hpp"
#include "vli/vli_cpu.hpp"
#include "vli/vli_gpu.hpp"
#include "vli/vli_traits.hpp"

#include "regression/vli_test.hpp"

using vli::vli_cpu;
using vli::max_int_value;
using vli::vli_gpu;

using vli::test::rnd_valid_int;
using vli::test::fill_random;




template <typename Vli>
struct vli_gpu_from_vli_cpu
{
    typedef vli_gpu<typename Vli::value_type, Vli::size> type;
};


typedef boost::mpl::transform<
      vli::test::vli_cpu_type_list
    , vli_gpu_from_vli_cpu<boost::mpl::_1>
    >::type vli_types;



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
