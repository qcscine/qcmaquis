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

#define BOOST_TEST_MODULE OverlapPropagatorVibrational

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/SweepBasedAlgorithms/OverlapPropagator.h"
#include "Fixtures/WatsonFixture.h"
#include "Fixtures/NModeFixture.h"

#ifdef DMRG_VIBRATIONAL

/**
 * @brief Verifies that the OverlapPropagator works for the Watson vibrational Hamiltonian 
 * The check is done by verifying that the overlap returned by [OverlapPropagator] is the
 * same as the one obtained by calling overlap(mps1, mps2)
 */
BOOST_FIXTURE_TEST_CASE(Test_OverlapPropagator_Vibrational_Watson, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
  // Data Generation
  using OverlapPropagatorType = OverlapPropagator<matrix, TrivialGroup, storage::disk>;
  using MPSType = MPS<matrix, TrivialGroup>;
  auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
  auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatsonHarmonic);
  auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
  parametersEthyleneWatsonHarmonic.set("init_state", "default");
  parametersEthyleneWatsonHarmonic.set("seed", 30031989);
  auto mpsDefault = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
  parametersEthyleneWatsonHarmonic.set("init_state", "const");
  auto mpsConst = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
  mpsDefault.normalize_right();
  mpsConst.normalize_right();
  // Constructs the overlap propagator
  std::vector<MPSType> mpsVector;
  mpsVector.push_back(mpsConst);
  auto overlapPropagator = OverlapPropagatorType(mpsDefault, mpsVector, 0);
  auto referenceOverlap = overlap(mpsConst, mpsDefault);
  for (int iSite = 0; iSite < 12; iSite++) {
    auto orthoMPS = overlapPropagator.getOrthogonalVector(0, iSite);
    overlapPropagator.propagateLeftOverlapBoundaries(iSite, iSite+1);
    auto overlap = ietl::dot(orthoMPS, mpsDefault[iSite]);
    BOOST_CHECK_CLOSE(overlap, referenceOverlap, 1.0E-8);
  }
#endif // HAVE_TrivialGroup
}

#ifdef HAVE_NU1
 
BOOST_FIXTURE_TEST_CASE(Test_OverlapPropagator_Vibrational_NU1, NModeFixture)
{
  using NU1SymmGroup = NU1_template<2>;
  using MPSType = MPS<matrix, NU1SymmGroup>;
  using OverlapPropagatorType = OverlapPropagator<matrix, NU1SymmGroup, storage::disk>;
  parametersFADTwoBody.set("init_state", "const");
  auto nModeLattice = Lattice(parametersFADTwoBody);
  auto nModeModel = Model<matrix, NU1SymmGroup>(nModeLattice, parametersFADTwoBody);
  auto nModeMPO = make_mpo(nModeLattice, nModeModel);
  auto mpsDefault = MPS<matrix, NU1SymmGroup>(nModeLattice.size(), *(nModeModel.initializer(nModeLattice, parametersFADTwoBody)));
  auto mpsVector = std::vector<MPSType>({mpsDefault});
  auto overlapPropagator = OverlapPropagatorType(mpsDefault, mpsVector, 0);
  overlapPropagator.propagateLeftOverlapBoundaries(0, 7);
  auto mpsOrtho = overlapPropagator.getOrthogonalVector(0, 7);
  auto overlap = ietl::dot(mpsOrtho, mpsDefault[7]);
  BOOST_CHECK_CLOSE(overlap, 1., 1.0E-7);
}

#endif // HAVE_NU1

#endif // DMRG_VIBRATIONAL