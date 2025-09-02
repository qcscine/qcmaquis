/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN
// check includes first
#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <complex>
#include <iostream>
// #include "dmrg/models/generate_mpo/1D_mpo_maker.hpp"
// #include "dmrg/mp_tensors/mpo_times_mps.hpp"
// #include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/model.h"
// #include "dmrg/sim/matrix_types.h"
// #include "dmrg/sim/matrix_types.h"
// #include "dmrg/mp_tensors/mpo.h"
// #include "dmrg/mp_tensors/mps.h"
// #include "dmrg/models/generate_mpo.hpp"
//

typedef boost::mpl::list<
#ifdef HAVE_TrivialGroup
    TrivialGroup
#ifdef HAVE_U1
    ,
#endif
#ifdef HAVE_NU1
    ,
#endif
#endif
#ifdef HAVE_U1
    U1
#ifdef HAVE_NU1
    ,
#endif
#endif
#ifdef HAVE_NU1
    NU1_template<5>
#endif
    >
    symmetries;

BOOST_AUTO_TEST_CASE_TEMPLATE(Test_RegisterOperatorsWithTags, S, symmetries) {
  // type defs
  using Matrix = alps::numeric::matrix<std::complex<double>>;
  using base = model_impl<Matrix, S>;
  using op_t = typename base::op_t;
  using table_type = typename base::table_type;
  using tag_type = typename base::tag_type;
  using value_type = typename Matrix::value_type;
  // creating matrices
  alps::numeric::matrix<std::complex<double>> identity(4, 4, 0.),
      identityScaled(4, 4, 0.), notIdentity(4, 4, 0.);
  double scalingFactor = 2.;
  identity(0, 0) = 1.;
  identityScaled(0, 0) = scalingFactor;
  for (int n = 1; n < 4; n++) {
    identity(n, n) = 1.;
    identityScaled(n, n) = scalingFactor;
    notIdentity(n, n - 1) = 1.;
  }
  op_t identity_oploc, identityScaled_oploc, notIdentity_oploc;
  // conditionals based on charge types
  // TrivialGroup
  if (symm_traits::SymmetryNameTrait<S>::symmName() == "none") {
    typename S::charge C = S::IdentityCharge;
    identity_oploc.insert_block(identity, C, C);
    identityScaled_oploc.insert_block(identityScaled, C, C);
    notIdentity_oploc.insert_block(notIdentity, C, C);
  }
  // other vibrational symmetries
  if (symm_traits::SymmetryNameTrait<S>::symmName() == "u1" ||
      symm_traits::SymmetryNameTrait<S>::symmName() == "nu1") {
    typename S::charge C0 = typename S::charge(0);
    typename S::charge C1 = typename S::charge(1);
    identity_oploc.insert_block(identity, C1, C1);
    identityScaled_oploc.insert_block(identityScaled, C1, C1);
    notIdentity_oploc.insert_block(notIdentity, C1, C0);
  }
  // pointer to the tag handler
  std::shared_ptr<TagHandler<Matrix, S>> tag_handler;
  tag_handler = std::make_shared<table_type>();
  // if bosonic
  auto TagIdentity =
      tag_handler->register_op(identity_oploc, tag_detail::bosonic);
  auto TagIdentityScaled =
      tag_handler->register_op(identityScaled_oploc, tag_detail::bosonic);
  auto TagNotIdentity =
      tag_handler->register_op(notIdentity_oploc, tag_detail::bosonic);
  // checks whether tags are differnet
  BOOST_CHECK_NE(TagIdentity, TagIdentityScaled);
  BOOST_CHECK_NE(TagIdentity, TagNotIdentity);
  BOOST_CHECK_NE(TagIdentityScaled, TagNotIdentity);
  // checking register
  bool isIdentityPesent = tag_handler->hasRegistered(identity_oploc);
  bool isIdentityScaledPesent =
      tag_handler->hasRegistered(identityScaled_oploc);
  bool isNotIdentityPesent = tag_handler->hasRegistered(notIdentity_oploc);
  // checks whether tags are correctly registered
  BOOST_CHECK_EQUAL(isIdentityPesent, true);
  BOOST_CHECK_EQUAL(isIdentityScaledPesent, true);
  BOOST_CHECK_EQUAL(isNotIdentityPesent, true);
}
