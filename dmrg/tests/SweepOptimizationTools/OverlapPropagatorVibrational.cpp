/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
  parametersEthyleneWatsonHarmonic.set("init_type", "default");
  parametersEthyleneWatsonHarmonic.set("seed", 30031989);
  auto mpsDefault = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
  parametersEthyleneWatsonHarmonic.set("init_type", "const");
  auto mpsConst = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
  mpsDefault.normalize_right();
  mpsConst.normalize_right();
  // Constructs the overlap propagator
  std::vector<MPSType> mpsVector;
  mpsVector.push_back(mpsConst);
  auto overlapPropagator = OverlapPropagatorType(mpsDefault, mpsVector, 0);
  auto referenceOverlap = overlap(mpsConst, mpsDefault);
  for (int iSite = 0; iSite < 12; iSite++) {
    auto orthoMPS = overlapPropagator.template getOrthogonalVector<SweepOptimizationType::SingleSite>(0, iSite, iSite+1);
    overlapPropagator.updateLeftOverlapBoundaries(iSite+1);
    auto overlap = ietl::dot(orthoMPS, mpsDefault[iSite]);
    BOOST_CHECK_CLOSE(overlap, referenceOverlap, 1.0E-8);
  }
#endif // HAVE_TrivialGroup
}

#endif // DMRG_VIBRATIONAL
