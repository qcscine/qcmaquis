/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE MPS_INITIALIZER_VIBRATIONAL

#ifdef DMRG_VIBRATIONAL
#include "Fixtures/NModeFixture.h"
#include "Fixtures/WatsonFixture.h"
#endif
#ifdef DMRG_VIBRONIC
#include "Fixtures/VibronicFixture.h"
#endif

#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps_initializers_helper.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/generate_mpo.hpp"

#include <iostream>
#include <boost/test/included/unit_test.hpp>

#ifdef DMRG_VIBRATIONAL

BOOST_FIXTURE_TEST_CASE(Test_Vibrational_Initializer_Helper_NU1, NModeFixture)
{
#ifdef HAVE_NU1
  using Symmetry = NU1_template<2>;
  using IndexType = Index<Symmetry>;
  std::vector<IndexType> physCharges;
  std::vector<int> inputVec(2), siteTypes;
  inputVec[0] = 1;
  inputVec[1] = 2;
  // Populates the physical indices
  auto lattice = Lattice(parametersFADTwoBody);
  int latticeSize = lattice.size();
  auto nModeModel = Model<matrix, NU1_template<2>>(lattice, parametersFADTwoBody);
  for (int iSite = 0; iSite < latticeSize; iSite++)
    siteTypes.push_back(lattice.get_prop<int>("type", iSite));
  for (int iType = 0; iType < lattice.getMaxType(); iType++)
    physCharges.push_back(nModeModel.phys_dim(iType));
  auto outputVector = HelperClassBasisVectorConverter<Symmetry>::GenerateIndexFromString(parametersFADTwoBody, inputVec, physCharges,
                                                                                         siteTypes, latticeSize);
  for (int iSite = 0; iSite < outputVector.size(); iSite++) {
    if (iSite == 1) {
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite][0])[0], 1);
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite][0])[1], 0);
    }
    else if (iSite == 13) {
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite][0])[0], 0);
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite][0])[1], 1);
    }
    else {
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite][0])[0], 0);
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite][0])[1], 0);
    }
    BOOST_CHECK_EQUAL(boost::get<1>(outputVector[iSite][0]), 0);
  }
#endif
}

#ifdef HAVE_NU1

/** @brief Verifies that the energy obtained initializing the MPS with an ONV is correct */
BOOST_FIXTURE_TEST_CASE(Test_Vibrational_Initializer_OneMode_Energy_NU1, NModeFixture)
{
  using Symmetry = NU1_template<1>;
  parametersFADOneBody.set("init_type", "basis_state_generic");
  parametersFADOneBody.set("init_basis_state", "0");
  // Populates the physical indices
  auto lattice = Lattice(parametersFADOneBody);
  int latticeSize = lattice.size();
  auto nModeModel = Model<matrix, Symmetry>(lattice, parametersFADOneBody);
  auto mps = MPS<matrix, Symmetry>(latticeSize, *(nModeModel.initializer(lattice, parametersFADOneBody)));
  auto mpo = make_mpo(lattice, nModeModel);
  auto energy = expval(mps, mpo)/norm(mps);
  // The energy is taken from the integral provides as input in the fixture class.
  BOOST_CHECK_CLOSE(energy, -2.359242429009664e+03, 1.0E-10);
}

/** @brief Verifies that the energy obtained initializing the MPS with an ONV is correct */
BOOST_FIXTURE_TEST_CASE(Test_Vibrational_Initializer_OneMode_Energy_FromBinary_NU1, NModeFixture)
{
  using Symmetry = NU1_template<1>;
  parametersFADOneBodyBinary.set("init_type", "basis_state_generic");
  parametersFADOneBodyBinary.set("init_basis_state", "10");
  // Populates the physical indices
  auto lattice = Lattice(parametersFADOneBodyBinary);
  int latticeSize = lattice.size();
  auto nModeModel = Model<matrix, Symmetry>(lattice, parametersFADOneBodyBinary);
  auto mps = MPS<matrix, Symmetry>(latticeSize, *(nModeModel.initializer(lattice, parametersFADOneBodyBinary)));
  auto mpo = make_mpo(lattice, nModeModel);
  auto energy = expval(mps, mpo)/norm(mps);
  // The energy is taken from the integral provides as input in the fixture class.
  BOOST_CHECK_CLOSE(energy, 1.408367346423375e+03, 1.0E-10);
}

/** @brief Same as above, but for the two-mode PESs */
BOOST_FIXTURE_TEST_CASE(Test_Vibrational_Initializer_TwoMode_Energy_NU1, NModeFixture)
{
  using Symmetry = NU1_template<2>;
  parametersFADTwoBody.set("init_type", "basis_state_generic");
  parametersFADTwoBody.set("init_basis_state", "2,3");
  // Populates the physical indices
  auto lattice = Lattice(parametersFADTwoBody);
  int latticeSize = lattice.size();
  auto nModeModel = Model<matrix, Symmetry>(lattice, parametersFADTwoBody);
  auto mps = MPS<matrix, Symmetry>(latticeSize, *(nModeModel.initializer(lattice, parametersFADTwoBody)));
  auto mpo = make_mpo(lattice, nModeModel);
  auto energy = expval(mps, mpo)/norm(mps);
  auto refEnergy = 6.996161115711967e+02 + 1.801678060826892e+03 - 2.258583526759012e+01;
  // The energy is taken from the integral provides as input in the fixture class.
  BOOST_CHECK_CLOSE(energy, refEnergy, 1.0E-10);
}

/** @brief Verifies that changing the modals order does not alter the energy */
BOOST_FIXTURE_TEST_CASE(Test_Vibrational_Initializer_TwoMode_NU1_ArbitrarySorting, NModeFixture)
{
  using Symmetry = NU1_template<5>;
  // Conventional sorting
  parametersFADTwoBodyFingerPrint.set("init_type", "basis_state_generic");
  parametersFADTwoBodyFingerPrint.set("init_basis_state", "0,0,0,0,0");
  auto lattice = Lattice(parametersFADTwoBodyFingerPrint);
  auto nModeModel = Model<matrix, Symmetry>(lattice, parametersFADTwoBodyFingerPrint);
  auto mpo = make_mpo(lattice, nModeModel);
  auto mps = MPS<matrix, Symmetry>(lattice.size(), *(nModeModel.initializer(lattice, parametersFADTwoBodyFingerPrint)));
  auto energy1 = expval(mps, mpo)/overlap(mps, mps);
  // Random sorting
  parametersFADTwoBodyFingerPrint.set("modals_order", "11,3,5,10,19,0,14,4,8,13,18,17,2,12,9,1,6,7,15,16");
  auto latticeFiedler = Lattice(parametersFADTwoBodyFingerPrint);
  auto nModeModelFiedler = Model<matrix, Symmetry>(latticeFiedler, parametersFADTwoBodyFingerPrint);
  auto mpoFiedler = make_mpo(latticeFiedler, nModeModelFiedler);
  auto mpsFiedler = MPS<matrix, Symmetry>(latticeFiedler.size(), *(nModeModelFiedler.initializer(latticeFiedler, parametersFADTwoBodyFingerPrint)));
  auto energy2 = expval(mpsFiedler, mpoFiedler)/overlap(mpsFiedler, mpsFiedler);
  // The energy is taken from the integral provides as input in the fixture class.
  BOOST_CHECK_CLOSE(energy1, energy2, 1.0E-10);
}

#endif // HAVE_NU1

#ifdef HAVE_TrivialGroup

/** @brief Tests the coherent initialization of an MPS */
BOOST_FIXTURE_TEST_CASE(Test_Vibrational_Initializer_Coherent, WatsonFixture)
{
  using Symmetry = TrivialGroup;
  auto lattice = Lattice(parametersEthyleneWatson);
  int latticeSize = lattice.size();
  auto watsonModel = Model<matrix, Symmetry>(lattice, parametersEthyleneWatson);
  auto mpo = make_mpo(lattice, watsonModel);
  // Construction of the ground state
  parametersEthyleneWatson.set("init_type", "basis_state_generic");
  parametersEthyleneWatson.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
  auto mpsGS = MPS<matrix, Symmetry>(latticeSize, *(watsonModel.initializer(lattice, parametersEthyleneWatson)));
  // Construction of the excited state
  parametersEthyleneWatson.set("init_type", "basis_state_generic");
  parametersEthyleneWatson.set("init_basis_state", "1,0,0,0,0,0,0,0,0,0,0,0");
  auto mpsES = MPS<matrix, Symmetry>(latticeSize, *(watsonModel.initializer(lattice, parametersEthyleneWatson)));
  // Construction of the coherent superposition
  parametersEthyleneWatson.set("init_type", "coherent");
  parametersEthyleneWatson.set("init_coeffs", "0.5,0.5");
  parametersEthyleneWatson.set("init_bond_dimension", 6);
  parametersEthyleneWatson.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0|1,0,0,0,0,0,0,0,0,0,0,0");
  auto mpsCoherent = MPS<matrix, Symmetry>(latticeSize, *(watsonModel.initializer(lattice, parametersEthyleneWatson)));
  // The energy is taken from the integral provides as input in the fixture class.
  auto energyGS = expval(mpsGS, mpo)/norm(mpsGS);
  auto energyES = expval(mpsES, mpo)/norm(mpsES);
  auto energyCoherent = expval(mpsCoherent, mpo)/norm(mpsCoherent);
  //
  BOOST_CHECK_CLOSE(energyGS+energyES, 2*energyCoherent, 1.0E-10);
}


#endif // HAVE_NONE
#endif // DMRG_VIBRATIONAL


#ifdef DMRG_VIBRONIC

/** @brief Tests the coherent initialization of an MPS in the excitonicextended model */
BOOST_FIXTURE_TEST_CASE(Test_Vibrational_Initializer_Coherent_ExitonicExtended, VibronicFixture)
{
#ifdef HAVE_U1
  using Symmetry = U1;
  auto lattice = Lattice(parametersSimpleCoherent);
  int latticeSize = lattice.size();
  auto eeModel = Model<matrix, Symmetry>(lattice, parametersSimpleCoherent); //excitonicextended (EE) model
  auto mpo = make_mpo(lattice, eeModel);
  // Construction of the first state (excitation on monomer one)
  parametersSimpleCoherent.set("init_type", "basis_state_generic");
  parametersSimpleCoherent.set("init_basis_state", "1,0,0,0");
  auto mpsStateOne = MPS<matrix, Symmetry>(latticeSize, *(eeModel.initializer(lattice, parametersSimpleCoherent)));
  // Construction of the second state (excitation on monomer two)
  parametersSimpleCoherent.set("init_type", "basis_state_generic");
  parametersSimpleCoherent.set("init_basis_state", "0,0,1,0");
  auto mpsStateTwo = MPS<matrix, Symmetry>(latticeSize, *(eeModel.initializer(lattice, parametersSimpleCoherent)));
  // Construction of the coherent superposition
  parametersSimpleCoherent.set("init_type", "coherent");
  parametersSimpleCoherent.set("init_coeffs", "0.5,0.5");
  parametersSimpleCoherent.set("init_bond_dimension", 5);
  parametersSimpleCoherent.set("init_basis_state", "1,0,0,0|0,0,1,0");
  auto mpsCoherent = MPS<matrix, Symmetry>(latticeSize, *(eeModel.initializer(lattice, parametersSimpleCoherent)));
  // The energy is taken from the integral provides as input in the fixture class.
  auto energyOne = expval(mpsStateOne, mpo)/norm(mpsStateOne);
  auto energyTwo = expval(mpsStateTwo, mpo)/norm(mpsStateTwo);
  auto energyCoherent = expval(mpsCoherent, mpo)/norm(mpsCoherent);
  //
  BOOST_CHECK_CLOSE(energyOne+energyTwo, 2*energyCoherent, 1.0E-10);
#endif // HAVE_U1
}


#ifdef HAVE_TrivialGroup

BOOST_FIXTURE_TEST_CASE(Test_RegisterOperatorsWithTags_TrivialGroup, VibronicFixture)
{
  //create instance of tag_handler
  //generating 4x4 matrices
  alps::numeric::matrix<std::complex<double>> identity(4, 4, 0.), identityScaled(4, 4, 0.), notIdentity(4, 4, 0.);
  double scalingFactor = 2.;
  TrivialGroup::charge C = TrivialGroup::IdentityCharge;
  identity(0, 0) = 1.;
  identityScaled(0, 0) = scalingFactor;
  for (int n = 1; n < 4; n++){
    identity(n,n) = 1.;
    identityScaled(n,n) = scalingFactor;
    notIdentity(n,n-1) = 1.;
  }
  model_impl<alps::numeric::matrix<std::complex<double>>, TrivialGroup>::op_t identity_oploc, identityScaled_oploc, notIdentity_oploc;
  identity_oploc.insert_block(identity, C, C);
  identityScaled_oploc.insert_block(identityScaled, C, C);
  notIdentity_oploc.insert_block(notIdentity, C, C);
  //pointer to the tag handler
  std::shared_ptr<TagHandler<alps::numeric::matrix<std::complex<double>>, TrivialGroup>> tag_handler; //from chatGPT
  tag_handler = std::make_shared<model_impl<alps::numeric::matrix<std::complex<double>>, TrivialGroup>::table_type>(); //from chatGPT
  auto TagIdentity = tag_handler->register_op(identity_oploc, tag_detail::bosonic);
  auto TagIdentityScaled = tag_handler->register_op(identityScaled_oploc, tag_detail::bosonic);
  auto TagNotIdentity = tag_handler->register_op(notIdentity_oploc, tag_detail::bosonic);
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

#endif // HAVE_NONE
#endif // DMRG_VIBRONIC
