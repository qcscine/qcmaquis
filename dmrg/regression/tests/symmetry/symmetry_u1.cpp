#define BOOST_TEST_MODULE symmetry

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/geometry/geometries/adapted/boost_array.hpp>
#include "dmrg/block_matrix/symmetry.h"

BOOST_AUTO_TEST_CASE(symmetry_u1_type_test){

  typedef U1::charge value_type;
  BOOST_MPL_ASSERT(( boost::is_same< int, value_type> ));  
  value_type q(0); 
  BOOST_CHECK_EQUAL(q,0);
}

BOOST_AUTO_TEST_CASE(symmetry_u1_fuse_test){

  typedef U1::charge value_type;
  value_type q0(0); 
  value_type q1(1); 
  value_type q2(2); 
  value_type qres; 

  boost::array<value_type,3> array = { { q0,q1,q2 } };

  qres = U1::fuse<3>(array);

  BOOST_CHECK_EQUAL(qres,3);
}
