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

#define BOOST_TEST_MODULE VibrationalFEAST

#include <boost/test/included/unit_test.hpp>
#include "Fixtures/WatsonFixture.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/MetaSweepSimulations/FEASTSimulator.h"
#include "dmrg/sim/matrix_types.h"
#include "utils/fpcomparison.h"

BOOST_FIXTURE_TEST_CASE(Test_FEAST_MPS_Getter, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
  using FEASTSimulatorType = FEASTSimulator<TrivialGroup>;
  using ModelType = Model<cmatrix, TrivialGroup>;
  parametersH2COWatson.set("max_bond_dimension", 10);
  parametersH2COWatson.set("feast_num_states", 2);
  parametersH2COWatson.set("feast_max_iter", 5);
  parametersH2COWatson.set("feast_emin", 0.);
  parametersH2COWatson.set("feast_emax", 6000.);
  parametersH2COWatson.set("feast_num_points", 8);
  parametersH2COWatson.set("feast_init_type", "basis_state_generic");
  parametersH2COWatson.set("feast_init_onv", "0,0,0,0,0,0|1,0,0,0,0,0");
  auto vibrationalLattice = Lattice(parametersH2COWatson);
  auto vibrationalModel = ModelType(vibrationalLattice, parametersH2COWatson);
  auto feastSimulator = FEASTSimulatorType(parametersH2COWatson, vibrationalModel, vibrationalLattice);
  // Checks consistency between guess MPS
  auto firstMPS = feastSimulator.getCurrentGuess(0);
  auto secondMPS = feastSimulator.getCurrentGuess(1);
  auto overlapBetweenMPS = overlap(firstMPS, secondMPS);
  BOOST_CHECK_SMALL(maquis::real(overlapBetweenMPS), 1.0E-15);
  // Checks consistency for quadrature points
  auto quadPoints = feastSimulator.getQuadraturePoints();
  BOOST_CHECK_EQUAL(quadPoints.size(), 8);
#endif // HAVE_TrivialGroup
}

#ifdef HAVE_TrivialGroup

/*
BOOST_FIXTURE_TEST_CASE(Test_FEAST_H2CO, WatsonFixture)
{
  using FEASTSimulatorType = FEASTSimulator<TrivialGroup>;
  using ModelType = Model<cmatrix, TrivialGroup>;
  //
  parametersH2COWatson.set("max_bond_dimension", 10);
  parametersH2COWatson.set("optimization", "twosite");
  parametersH2COWatson.set("symmetry", "none");
  parametersH2COWatson.set("init_state", "const");
  parametersH2COWatson.set("nsweeps", 3);
  parametersH2COWatson.set("chkpfile", "GS.H2CO.chkp.h5");
  maquis::DMRGInterface<double> interfaceOptimizerGS(parametersH2COWatson);
  interfaceOptimizerGS.optimize();
  auto energyFromOptimizerGS = interfaceOptimizerGS.energy();
  // Excited-state calculation
  parametersH2COWatson.set("n_ortho_states", 1);
  parametersH2COWatson.set("ortho_states", "GS.H2CO.chkp.h5");
  parametersH2COWatson.set("chkpfile", "ES.H2CO.chkp.h5");
  maquis::DMRGInterface<double> interfaceOptimizerES(parametersH2COWatson);
  interfaceOptimizerES.optimize();
  auto energyFromOptimizerES = interfaceOptimizerES.energy();
  // Cleans up stuff
  boost::filesystem::remove_all("GS.H2CO.chkp.h5");
  boost::filesystem::remove_all("ES.H2CO.chkp.h5");
  // FEAST
  auto eMin = energyFromOptimizerGS - (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  auto eMax = energyFromOptimizerGS + (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  parametersH2COWatson.set("feast_num_states", 1);
  parametersH2COWatson.set("feast_max_iter", 1);
  parametersH2COWatson.set("feast_emin", eMin);
  parametersH2COWatson.set("feast_emax", eMax);
  parametersH2COWatson.set("feast_num_points", 2);
  parametersH2COWatson.set("feast_init_type", "default");
  parametersH2COWatson.set("nsweeps", 1);
  auto vibrationalLattice = Lattice(parametersH2COWatson);
  auto vibrationalModel = ModelType(vibrationalLattice, parametersH2COWatson);
  auto vibrationalMPO = make_mpo(vibrationalLattice, vibrationalModel);
  auto feastSimulator = FEASTSimulatorType(parametersH2COWatson, vibrationalModel, vibrationalLattice);
  feastSimulator.runFeastSimulation(vibrationalMPO);
}
*/

BOOST_FIXTURE_TEST_CASE(Test_FEAST_Bilinearly, WatsonFixture)
{
  using FEASTSimulatorType = FEASTSimulator<TrivialGroup>;
  using ModelType = Model<cmatrix, TrivialGroup>;
  //
  parametersBilinearly.set("init_state", "basis_state_generic");
  parametersBilinearly.set("init_basis_state", "0,0,0,0,0,0");
  parametersBilinearly.set("optimization", "singlesite");
  parametersBilinearly.set("alpha_initial", 1.0E-8);
  parametersBilinearly.set("alpha_initial", 1.0E-15);
  parametersBilinearly.set("alpha_initial", 0.);
  parametersBilinearly.set("nsweeps", 20);
  parametersBilinearly.set("ngrowsweeps", 2);
  parametersBilinearly.set("nmainsweeps", 2);
  parametersBilinearly.set("max_bond_dimension", 20);
  parametersBilinearly.set("MODEL", "watson");
  parametersBilinearly.set("chkpfile", "GS.Bilinearly.chkp.h5");
  maquis::DMRGInterface<double> interfaceOptimizerGS(parametersBilinearly);
  interfaceOptimizerGS.optimize();
  auto energyFromOptimizerGS = interfaceOptimizerGS.energy();
  // Excited-state calculation
  parametersBilinearly.set("n_ortho_states", 1);
  parametersBilinearly.set("ortho_states", "GS.Bilinearly.chkp.h5");
  parametersBilinearly.set("chkpfile", "ES.Bilinearly.chkp.h5");
  maquis::DMRGInterface<double> interfaceOptimizerES(parametersBilinearly);
  interfaceOptimizerES.optimize();
  auto energyFromOptimizerES = interfaceOptimizerES.energy();
  // Cleans up stuff
  boost::filesystem::remove_all("GS.Bilinearly.chkp.h5");
  boost::filesystem::remove_all("ES.Bilinearly.chkp.h5");
  // == DMRG[FEAST] ==
  auto eMin = energyFromOptimizerGS - (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  auto eMax = energyFromOptimizerGS + (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  // FEAST-specific parameters
  parametersBilinearly.set("feast_num_states", 1);
  parametersBilinearly.set("feast_max_iter", 1);
  parametersBilinearly.set("feast_emin", eMin);
  parametersBilinearly.set("feast_emax", eMax);
  parametersBilinearly.set("feast_num_points", 8);
  parametersBilinearly.set("feast_init_type", "basis_state_generic_default");
  parametersBilinearly.set("feast_init_onv", "2,0,1,2,1,2");
  // parametersH2COWatson.set("feast_truncation_type", "end");
  // Setup parameters for the linear system solver.
  parametersBilinearly.set("linsystem_precond", "no");
  parametersBilinearly.set("linsystem_krylov_dim", 50);
  parametersBilinearly.set("linsystem_tol", 1.0E-5);
  parametersBilinearly.set("linsystem_init", "last");
  auto vibrationalLattice = Lattice(parametersBilinearly);
  auto vibrationalModel = ModelType(vibrationalLattice, parametersBilinearly);
  auto vibrationalMPO = make_mpo(vibrationalLattice, vibrationalModel);
  auto feastSimulator = FEASTSimulatorType(parametersBilinearly, vibrationalModel, vibrationalLattice);
  feastSimulator.runFeastSimulation(vibrationalMPO);
  auto feastEnergy = feastSimulator.getVibrationalEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromOptimizerGS, 1.0E-6);
}

#endif // HAVE_TrivialGroup
