/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE SweepBasedEnergyMinimizationVibrational

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/SweepBasedAlgorithms/SweepBasedEnergyMinimization.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/WatsonFixture.h"

#ifdef DMRG_VIBRATIONAL

/**
 * @brief Checks that the constructor for the sweep-based energy minimization works.
 */
BOOST_FIXTURE_TEST_CASE(Test_SweepBasedEnergyMinimizerSS_Vibrational_Watson, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
  // Data Generation
  using SweepBasedMinimizerSS = SweepBasedEnergyMinimization<matrix, TrivialGroup, storage::disk, SweepOptimizationType::SingleSite>;
  using MPSType = MPS<matrix, TrivialGroup>;
  parametersEthyleneWatsonHarmonic.set("nsweeps", 10);
  parametersEthyleneWatsonHarmonic.set("max_bond_dimension", 100);
  parametersEthyleneWatsonHarmonic.set("alpha_initial", 1.0E-8);
  parametersEthyleneWatsonHarmonic.set("alpha_main", 1.0E-15);
  parametersEthyleneWatsonHarmonic.set("alpha_final", 0.);
  auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
  auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatsonHarmonic);
  auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
  parametersEthyleneWatsonHarmonic.set("init_type", "default");
  parametersEthyleneWatsonHarmonic.set("seed", 30031989);
  auto mpsDefault = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
  mpsDefault.normalize_right();
  auto energyMinimizer = SweepBasedMinimizerSS(mpsDefault, watsonHarmonicMPO, parametersEthyleneWatsonHarmonic, watsonModel, lattice, false);
  energyMinimizer.runSweepSimulation();
  double optimalEnergyFromSweeper = energyMinimizer.getSpecificResult<double>("Energy");
  // Now does the same with the interface
  parametersEthyleneWatsonHarmonic.set("optimization", "singlesite");
  maquis::DMRGInterface<double> interfaceWatson(parametersEthyleneWatsonHarmonic);
  interfaceWatson.optimize();
  double optimalEnergyFromInterface = interfaceWatson.energy();
  BOOST_CHECK_CLOSE(optimalEnergyFromInterface, optimalEnergyFromSweeper, 1.0e-7);
#endif // HAVE_TrivialGroup
}

#ifdef HAVE_TrivialGroup

BOOST_FIXTURE_TEST_CASE(Test_SweepBasedEnergyMinimizerTS_Vibrational_Bilinearly, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
  // Data Generation
  using SweepBasedMinimizerTS = SweepBasedEnergyMinimization<matrix, TrivialGroup, storage::disk, SweepOptimizationType::TwoSite>;
  using MPSType = MPS<matrix, TrivialGroup>;
  parametersBilinearly.set("nsweeps", 10);
  parametersBilinearly.set("max_bond_dimension", 20);
  auto lattice = Lattice(parametersBilinearly);
  auto bilinearlyModel = Model<matrix, TrivialGroup>(lattice, parametersBilinearly);
  auto bilinearlyMPO = make_mpo(lattice, bilinearlyModel);
  parametersBilinearly.set("init_type", "const");
  // TSOptimizer calculation
  auto mpsConst = MPS<matrix, TrivialGroup>(lattice.size(), *(bilinearlyModel.initializer(lattice, parametersBilinearly)));
  mpsConst.normalize_right();
  auto energyMinimizer = SweepBasedMinimizerTS(mpsConst, bilinearlyMPO, parametersBilinearly, bilinearlyModel, lattice, false);
  energyMinimizer.runSweepSimulation();
  double optimalEnergyFromSweeper = energyMinimizer.getSpecificResult<double>("Energy");
  // Interface calculation
  parametersBilinearly.set("optimization", "twosite");
  maquis::DMRGInterface<double> interfaceBilinearly(parametersBilinearly);
  interfaceBilinearly.optimize();
  double optimalEnergyFromInterface = interfaceBilinearly.energy();
  // Final check
  BOOST_CHECK_CLOSE(optimalEnergyFromInterface, optimalEnergyFromSweeper, 1.0e-7);
#endif // HAVE_TrivialGroup
}

#endif // HAVE_TrivialGroup

#endif // DMRG_VIBRATIONAL
