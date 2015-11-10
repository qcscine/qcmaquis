/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

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

