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

#endif // DMRG_VIBRATIONAL
