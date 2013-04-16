#define BOOST_TEST_MODULE functional
#define DISTORT_SCALE 128
#define distort(x,y) x*DISTORT_SCALE+(y), x*DISTORT_SCALE-(y)
#include <utils/testing.hpp>

typedef boost::mpl::list< input<distort(5,-7),double,1>,
                          input<distort(5,-7),std::complex<double>,1>
                        > test_types; 
