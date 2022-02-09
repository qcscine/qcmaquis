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
  for (int iType = 0; iType <= lattice.maximum_vertex_type(); iType++)
    physCharges.push_back(nModeModel.phys_dim(iType));
  auto outputVector = HelperClassBasisVectorConverter<Symmetry>::GenerateIndexFromString(inputVec, physCharges,
                                                                                         siteTypes, latticeSize);
  for (int iSite = 0; iSite < outputVector.size(); iSite++) {
    if (iSite == 1) {
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite])[0], 1);
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite])[1], 0);
    }
    else if (iSite == 13) {
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite])[0], 0);
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite])[1], 1);
    }
    else {
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite])[0], 0);
      BOOST_CHECK_EQUAL(boost::get<0>(outputVector[iSite])[1], 0);
    }
    BOOST_CHECK_EQUAL(boost::get<1>(outputVector[iSite]), 0);
  }
#endif
}

#ifdef HAVE_NU1

/** @brief Verifies that the energy obtained initializing the MPS with an ONV is correct */
BOOST_FIXTURE_TEST_CASE(Test_Vibrational_Initializer_OneMode_Energy_NU1, NModeFixture)
{
  using Symmetry = NU1_template<1>;
  parametersFADOneBody.set("init_state", "basis_state_generic");
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
  parametersFADOneBodyBinary.set("init_state", "basis_state_generic");
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
  parametersFADTwoBody.set("init_state", "basis_state_generic");
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

#endif // HAVE_NU1

#endif // DMRG_VIBRATIONAL
