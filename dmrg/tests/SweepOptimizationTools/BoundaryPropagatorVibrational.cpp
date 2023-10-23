/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE BoundaryPropagatorVibrational

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/SweepBasedAlgorithms/BoundaryPropagator.h"
#include "Fixtures/WatsonFixture.h"
#include "Fixtures/NModeFixture.h"

#ifdef DMRG_VIBRATIONAL

/**
 * @brief Verifies that the BoundaryPropagator works for the Watson vibrational Hamiltonian 
 * The check is done by verifying that, once feeded to the SiteProblem constructor, the
 * correct energy is retrieved
 */
BOOST_FIXTURE_TEST_CASE(Test_BoundaryPropagator_Vibrational_Watson, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
  using BoundaryPropagatorType = BoundaryPropagator<matrix, TrivialGroup, storage::disk>;
  using SiteProblemType = SiteProblem<matrix, TrivialGroup>;
  parametersEthyleneWatsonHarmonic.set("init_type", "default");
  auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
  auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatsonHarmonic);
  auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
  auto mpsDefault = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
  mpsDefault.normalize_right();
  mpsDefault.canonize(5);
  auto boundaryPropagator = BoundaryPropagatorType(mpsDefault, watsonHarmonicMPO, 5);
  auto siteProblem = SiteProblemType(boundaryPropagator.getLeftBoundary(5), boundaryPropagator.getRightBoundary(6), watsonHarmonicMPO[5]);
  auto energyFromSP = siteProblem.get_energy(mpsDefault[5]);
  auto energyFromExpval = expval(mpsDefault, watsonHarmonicMPO);
  BOOST_CHECK_CLOSE(energyFromSP, energyFromExpval, 1e-7);
#endif // HAVE_TrivialGroup
}

#ifdef HAVE_NU1

BOOST_FIXTURE_TEST_CASE(Test_BoundaryPropagator_Vibrational_NU1, NModeFixture)
{
  using NU1SymmGroup = NU1_template<2>;
  using BoundaryPropagatorType = BoundaryPropagator<matrix, NU1SymmGroup, storage::disk>;
  using SiteProblemType = SiteProblem<matrix, NU1SymmGroup>;
  parametersFADTwoBody.set("init_type", "const");
  auto nModeLattice = Lattice(parametersFADTwoBody);
  auto nModeModel = Model<matrix, NU1SymmGroup>(nModeLattice, parametersFADTwoBody);
  auto nModeMPO = make_mpo(nModeLattice, nModeModel);
  auto mpsDefault = MPS<matrix, NU1SymmGroup>(nModeLattice.size(), *(nModeModel.initializer(nModeLattice, parametersFADTwoBody)));
  mpsDefault.normalize_right();
  mpsDefault.canonize(5);
  auto boundaryPropagator = BoundaryPropagatorType(mpsDefault, nModeMPO, 5);
  mpsDefault.canonize(7);
  boundaryPropagator.updateLeftBoundary(6);
  boundaryPropagator.updateLeftBoundary(7);
  auto siteProblem = SiteProblemType(boundaryPropagator.getLeftBoundary(7), boundaryPropagator.getRightBoundary(8), nModeMPO[7]);
  auto energyFromSP = siteProblem.get_energy(mpsDefault[7]);
  auto energyFromExpval = expval(mpsDefault, nModeMPO);
  BOOST_CHECK_CLOSE(energyFromSP, energyFromExpval, 1e-7);
}

#endif // HAVE_NU1

#endif // DMRG_VIBRATIONAL