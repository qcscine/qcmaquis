/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022 by Alberto Baiardi <abaiardi@ethz.ch>
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

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "dmrg/models/MolecularHamiltonians/IndexTuple.hpp"
#include "dmrg/block_matrix/symmetry.h"

BOOST_AUTO_TEST_CASE(Test_IndexTuple_Constructor_TwoU1) {
#ifdef HAVE_TwoU1
  typename chem::detail::IndexTuple<TwoU1, 4> index1({3, 2, 1, 1});
  BOOST_CHECK_EQUAL(index1.data().size(), 4);
#endif
}

#ifdef HAVE_TwoU1

BOOST_AUTO_TEST_CASE(Test_IndexTuple_Zero_Constructor_TwoU1) {
  typename chem::detail::IndexTuple<TwoU1, 4> index1;
  BOOST_CHECK_EQUAL(index1.data().size(), 4);
  for (int i = 0; i < 4; i++) BOOST_CHECK_EQUAL(index1[i], 0);
}

#endif

BOOST_AUTO_TEST_CASE(Test_IndexTuple_Constructor_SU2U1) {
#ifdef HAVE_SU2U1
  typename chem::detail::IndexTuple<SU2U1, 4> index1({3, 2, 1, 1});
  BOOST_CHECK_EQUAL(index1.data().size(), 4);
#endif
}

#ifdef HAVE_TwoU1

BOOST_AUTO_TEST_CASE(Test_IndexTuple_Merge) {
  typename chem::detail::IndexTuple<TwoU1, 4> index1({1, 2, 3, 4});
  typename chem::detail::IndexTuple<TwoU1, 4> index2({5, 6, 7, 8});
  typename chem::detail::IndexTuple<TwoU1, 8> indexMerged(index1, index2);
  for (int i = 0; i < 8; i++) BOOST_CHECK_EQUAL(indexMerged[i], i + 1);
}

BOOST_AUTO_TEST_CASE(Test_IndexTuple_FourElement_Align) {
  using TupleType = chem::detail::IndexTuple<TwoU1, 4>;
  TupleType indexUnsorted(std::array<int, 4>{3, 4, 8, 2});
  indexUnsorted.align(true);
  BOOST_CHECK_EQUAL(indexUnsorted[0], 8);
  BOOST_CHECK_EQUAL(indexUnsorted[1], 2);
  BOOST_CHECK_EQUAL(indexUnsorted[2], 4);
  BOOST_CHECK_EQUAL(indexUnsorted[3], 3);
}

BOOST_AUTO_TEST_CASE(Test_IndexTuple_SixElement_Align_1) {
  using TupleType = chem::detail::IndexTuple<TwoU1, 6>;
  TupleType indexUnsorted(std::array<int, 6>{3, 1, 4, 5, 6, 2});
  indexUnsorted.align(true);
  BOOST_CHECK_EQUAL(indexUnsorted[0], 6);
  BOOST_CHECK_EQUAL(indexUnsorted[1], 2);
  BOOST_CHECK_EQUAL(indexUnsorted[2], 5);
  BOOST_CHECK_EQUAL(indexUnsorted[3], 4);
  BOOST_CHECK_EQUAL(indexUnsorted[4], 3);
  BOOST_CHECK_EQUAL(indexUnsorted[5], 1);
}

BOOST_AUTO_TEST_CASE(Test_IndexTuple_SixElement_Align_2) {
  using TupleType = chem::detail::IndexTuple<TwoU1, 6>;
  TupleType indexUnsorted(std::array<int, 6>{3, 1, 2, 3, 3, 3});
  indexUnsorted.align(true);
  BOOST_CHECK_EQUAL(indexUnsorted[0], 3);
  BOOST_CHECK_EQUAL(indexUnsorted[1], 3);
  BOOST_CHECK_EQUAL(indexUnsorted[2], 3);
  BOOST_CHECK_EQUAL(indexUnsorted[3], 2);
  BOOST_CHECK_EQUAL(indexUnsorted[4], 3);
  BOOST_CHECK_EQUAL(indexUnsorted[5], 1);
}

#endif  // HAVE_TwoU1