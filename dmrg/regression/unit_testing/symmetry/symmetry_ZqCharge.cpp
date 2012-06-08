#define BOOST_TEST_MODULE symmetry

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/geometry/geometries/adapted/boost_array.hpp>
#include "dmrg/block_matrix/symmetry.h"

BOOST_AUTO_TEST_CASE(symmetry_ZqCharge_constructor_test){
  ZqCharge<0> q;
  BOOST_CHECK_EQUAL(q,0);
  ZqCharge<0> q0(0);
  BOOST_CHECK_EQUAL(q,q0);
  ZqCharge<0> q1(1);
  ZqCharge<0> q1bis(1);
  BOOST_CHECK_EQUAL(q1,q1bis);
}

BOOST_AUTO_TEST_CASE(symmetry_ZqCharge_operator_test){
  ZqCharge<1> q;
  ZqCharge<1> q0(0);
  ZqCharge<1> q1(1);
  ZqCharge<1> qres;
  BOOST_CHECK_EQUAL(true,q==q0);
  BOOST_CHECK_EQUAL(false,q==q1);
  BOOST_CHECK_EQUAL(false,q!=q0);
  BOOST_CHECK_EQUAL(true,q!=q1);
  BOOST_CHECK_EQUAL(true,q<q1);
  BOOST_CHECK_EQUAL(false,q1<q);
  qres = q0 + q1;
  BOOST_CHECK_EQUAL(q0,qres);
}
