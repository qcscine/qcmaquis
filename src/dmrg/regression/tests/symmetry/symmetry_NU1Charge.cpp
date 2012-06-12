#define BOOST_TEST_MODULE symmetry

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/geometry/geometries/adapted/boost_array.hpp>
#include "dmrg/block_matrix/symmetry.h"

BOOST_AUTO_TEST_CASE(symmetry_NU1Charge_constructor_test){
  NU1Charge<1> q;
  BOOST_CHECK_EQUAL(q[0],0);
  NU1Charge<2> q2(1);
  BOOST_CHECK_EQUAL(q2[0],1);
  BOOST_CHECK_EQUAL(q2[1],1);

  boost::array<int, 2> array = { { 1,1 } };
  NU1Charge<2> q3(array);
  BOOST_CHECK_EQUAL(q3[0],1);
  BOOST_CHECK_EQUAL(q3[1],1);
}

BOOST_AUTO_TEST_CASE(symmetry_NU1Charge_operator_test){

  NU1Charge<2> q1(1);
  NU1Charge<2> q2(1);
  NU1Charge<2> q22(2);

  BOOST_CHECK_EQUAL(true,q1==q2);
  BOOST_CHECK_EQUAL(false,q1==q22);

  BOOST_CHECK_EQUAL(false,q1!=q2);
  BOOST_CHECK_EQUAL(true,q1!=q22);

  BOOST_CHECK_EQUAL(false,q1<q2);
  BOOST_CHECK_EQUAL(true,q1<q22);

  BOOST_CHECK_EQUAL(false,q1>q2);
  BOOST_CHECK_EQUAL(false,q1>q22);
  BOOST_CHECK_EQUAL(true,q22>q1);
}

BOOST_AUTO_TEST_CASE(symmetry_NU1Charge_operator_plus_test){

  NU1Charge<2> qa(1);
  NU1Charge<2> qb(1);
  NU1Charge<2> qres;
  NU1Charge<2> qc(2);

  qres = qa + qb;

  BOOST_CHECK_EQUAL(qres,qc);
}

BOOST_AUTO_TEST_CASE(symmetry_NU1Charge_operator_negate_test){
  NU1Charge<2> qa(1);
  NU1Charge<2> mqa(-1);
  NU1Charge<2> qres;
  qres = -qa;
  BOOST_CHECK_EQUAL(qres,mqa);
}

BOOST_AUTO_TEST_CASE(symmetry_NU1Charge_operator_divide_test){
  NU1Charge<2> qa(4);
  NU1Charge<2> qb(2);
  NU1Charge<2> qres;
  qres = qa/2;
  BOOST_CHECK_EQUAL(qres,qb);
  qres = 2/qa;
  BOOST_CHECK_EQUAL(qres,qb);
}

