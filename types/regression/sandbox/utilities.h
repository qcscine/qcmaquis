#include <boost/mpl/push_front.hpp>


#define ValueWG 11
#define tuple1(z, n, unused) BOOST_PP_COMMA_IF(n) size<(BOOST_PP_ADD(n,1)),3,double> 
#define tuple2(z, n, unused) BOOST_PP_COMMA_IF(n) size<(BOOST_PP_ADD(n,1)),-3,double> 


template <int n, int m, typename T> // n # of workgroup, T double or std::complex<double> 
struct size {
   BOOST_STATIC_ASSERT(n>0);
   BOOST_STATIC_ASSERT(n*ValueWG > m);
   typedef T value_type; // To template later
   enum {valuex = n*ValueWG+m};// n is the number or work group, m how we resize
   enum {valuey = n*ValueWG-m};// n is the number or work group, m how we resize
   enum {null = 0};// n is the number or work group, m how we resize
};

//    typedef boost::mpl::list< size<4,0,double>, BOOST_PP_REPEAT(8,tuple1,~), BOOST_PP_REPEAT(8,tuple2,~) > test_types;
      typedef boost::mpl::list< size<1,-3,double> > test_types;

