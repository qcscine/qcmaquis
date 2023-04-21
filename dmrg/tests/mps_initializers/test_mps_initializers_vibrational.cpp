/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
*               2021 Alberto Baiardi <abaiardi@ethz.ch>
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

#define BOOST_TEST_MODULE MPS_INITIALIZER_VIBRATIONAL

#ifdef DMRG_VIBRATIONAL

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "Fixtures/NModeFixture.h"
#include "Fixtures/WatsonFixture.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps_initializers_helper.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/generate_mpo.hpp"

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
  parametersEthyleneWatson.set("init_coeff", "0.5,0.5");
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
