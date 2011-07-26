#define BOOST_TEST_MODULE vli_cpu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/list_c.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/push_front.hpp>

#include "vli_cpu/vli_number_cpu.hpp"
#include "vli_cpu/vli_number_traits.hpp"
#include "gmpxx.h"

using vli::vli_cpu;
using vli::max_int_value;

/**
  * Specify which BaseInt types and vli sizes ought to be tested
  */

typedef boost::mpl::list<unsigned int, unsigned long int> test_types;
typedef boost::mpl::list_c<std::size_t,2,4,8,16> test_sizes;


namespace detail
{
/**
  * Generate all combinations from test types and test sizes
  */
template <typename SizeList>
struct all_sizes
{
    template<typename MplList,typename Element>
    struct apply
    {
        typedef typename boost::mpl::fold<
            SizeList,
            MplList,
            boost::mpl::push_front<
                  typename boost::mpl::_1
                , typename boost::mpl::pair<Element,boost::mpl::_2>
                >
            >::type type;
    };
};

typedef boost::mpl::fold<
    test_types,
    boost::mpl::list<>,
    all_sizes<test_sizes>
    >::type type_size_pairs;

/**
  * Convert the Boost MPL Pairs pair<Base_int_type, size> to vli_cpu types.
  */
template <typename Pair>
struct vli_type
{
    typedef vli::vli_cpu<typename boost::mpl::first<Pair>::type,boost::mpl::second<Pair>::type::value> type;
};

typedef boost::mpl::transform<
        detail::type_size_pairs,
        vli_type<boost::mpl::_1>
        >::type vli_types;
} //namespace detail

typedef detail::vli_types vli_types;

#include "regression/vli_number_common_tests.hpp"
