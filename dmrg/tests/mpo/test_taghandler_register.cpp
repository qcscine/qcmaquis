/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#define BOOST_TEST_MAIN
 //check includes first
#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <complex>
#include "dmrg/models/generate_mpo/1D_mpo_maker.hpp"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps_initializers_helper.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/generate_mpo.hpp"

#include <iostream>
#include <boost/test/included/unit_test.hpp>
//


typedef boost::mpl::list<
#ifdef HAVE_TrivialGroup
TrivialGroup
#endif
#ifdef HAVE_U1
, U1
#endif
#ifdef HAVE_NU1
, NU1
#endif
#ifdef HAVE_TwoU1PG
, TwoU1PG
#endif
#ifdef HAVE_TwoU1
, TwoU1
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
#ifdef HAVE_SU2U1
, SU2U1
#endif
> symmetries;

#ifdef DMRG_VIBRONIC
#ifdef DMRG_VIBRATIONAL
BOOST_AUTO_TEST_CASE_TEMPLATE(Test_RegisterOperatorsWithTags, S, symmetries)
{
#ifdef HAVE_TrivialGroup
  alps::numeric::matrix<std::complex<double>> identity(4, 4, 0.), identityScaled(4, 4, 0.), notIdentity(4, 4, 0.);
  double scalingFactor = 2.;
  identity(0, 0) = 1.;
  identityScaled(0, 0) = scalingFactor;
  for (int n = 1; n < 4; n++){
    identity(n,n) = 1.;
    identityScaled(n,n) = scalingFactor;
    notIdentity(n,n-1) = 1.;
  }
  //trivial symmetry
  if(symm_traits::SymmetryNameTrait<S>::symmName() == "none"){
    model_impl<alps::numeric::matrix<std::complex<double>>, TrivialGroup>::op_t identity_oploc, identityScaled_oploc, notIdentity_oploc;
    TrivialGroup::charge C = TrivialGroup::IdentityCharge;
    identity_oploc.insert_block(identity, C, C);
    identityScaled_oploc.insert_block(identityScaled, C, C);
    notIdentity_oploc.insert_block(notIdentity, C, C);
    //pointer to the tag handler
    std::shared_ptr<TagHandler<alps::numeric::matrix<std::complex<double>>, TrivialGroup>> tag_handler; 
    tag_handler = std::make_shared<model_impl<alps::numeric::matrix<std::complex<double>>, TrivialGroup>::table_type>(); 
  } //restricted vibrational symmetries
  if(symm_traits::SymmetryNameTrait<S>::symmName() == "u1" || symm_traits::SymmetryNameTrait<S>::symmName() == "nu1" ){
    model_impl<alps::numeric::matrix<std::complex<double>>, S>::op_t identity_oploc, identityScaled_oploc, notIdentity_oploc;
    identity_oploc.insert_block(identity, 1, 1);
    identityScaled_oploc.insert_block(identityScaled, 1, 1);
    notIdentity_oploc.insert_block(notIdentity, 1, 0);
    //pointer to the tag handler
    std::shared_ptr<TagHandler<alps::numeric::matrix<std::complex<double>>, S>> tag_handler; 
    tag_handler = std::make_shared<model_impl<alps::numeric::matrix<std::complex<double>>, S>::table_type>(); 
  }
  else{ //electronic symmetries
    model_impl<alps::numeric::matrix<std::complex<double>>, S>::op_t identity_oploc, identityScaled_oploc, notIdentity_oploc;
    identity_oploc.insert_block(alps::numeric::matrix<std::complex<double>>(1, 1, 1), 0, 0);
    identityScaled_oploc.insert_block(alps::numeric::matrix<std::complex<double>>(scalingFactor, scalingFactor, scalingFactor), 0, 0);
    notIdentity_oploc.insert_block(alps::numeric::matrix<std::complex<double>>(1, 1, 1), 1, 1);
    //pointer to the tag handler
    std::shared_ptr<TagHandler<alps::numeric::matrix<std::complex<double>>, S>> tag_handler; 
    tag_handler = std::make_shared<model_impl<alps::numeric::matrix<std::complex<double>>, S>::table_type>();
  }
  //bosonic symmetries
  if(symm_traits::SymmetryNameTrait<S>::symmName() == "none" || symm_traits::SymmetryNameTrait<S>::symmName() == "u1" || symm_traits::SymmetryNameTrait<S>::symmName() == "nu1"){
    auto TagIdentity = tag_handler->register_op(identity_oploc, tag_detail::bosonic);
    auto TagIdentityScaled = tag_handler->register_op(identityScaled_oploc, tag_detail::bosonic);
    auto TagNotIdentity = tag_handler->register_op(notIdentity_oploc, tag_detail::bosonic);
  }
  //fermionic symmetries
  if(symm_traits::SymmetryNameTrait<S>::symmName() == "TwoU1PG" || symm_traits::SymmetryNameTrait<S>::symmName() == "TwoU1" || symm_traits::SymmetryNameTrait<S>::symmName() == "SU2U1PG" || symm_traits::SymmetryNameTrait<S>::symmName() == "SU2U1"){
    auto TagIdentity = tag_handler->register_op(identity_oploc, tag_detail::fermionic);
    auto TagIdentityScaled = tag_handler->register_op(identityScaled_oploc, tag_detail::fermionic);
    auto TagNotIdentity = tag_handler->register_op(notIdentity_oploc, tag_detail::fermionic);
  }
  //checks if alreagy registered
  bool isIdentityPesent = tag_handler->hasRegistered(identity_oploc);
  bool isIdentityScaledPesent = tag_handler->hasRegistered(identityScaled_oploc);
  bool isNotIdentityPesent = tag_handler->hasRegistered(notIdentity_oploc);
  BOOST_CHECK_EQUAL(isIdentityPesent, true);
  BOOST_CHECK_EQUAL(isIdentityScaledPesent, true);
  BOOST_CHECK_EQUAL(isNotIdentityPesent, true);
  //checks whether tags are differnet
  BOOST_CHECK_NE(TagIdentity, TagIdentityScaled);
  BOOST_CHECK_NE(TagIdentity, TagNotIdentity);
  BOOST_CHECK_NE(TagIdentityScaled, TagNotIdentity);
}

