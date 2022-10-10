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
#include "Fixtures/NModeFixture.h"
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
  auto vibrationalMPO = make_mpo(vibrationalLattice, vibrationalModel);
  auto feastSimulator = FEASTSimulatorType(parametersH2COWatson, vibrationalModel, vibrationalLattice, vibrationalMPO);
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

BOOST_FIXTURE_TEST_CASE(Test_FEAST_H2CO, WatsonFixture)
{
  using FEASTSimulatorType = FEASTSimulator<TrivialGroup>;
  using ModelType = Model<cmatrix, TrivialGroup>;
  //
  parametersH2COWatsonNoCoriolis.set("max_bond_dimension", 50);
  parametersH2COWatsonNoCoriolis.set("init_state", "const");
  parametersH2COWatsonNoCoriolis.set("optimization", "singlesite");
  parametersH2COWatsonNoCoriolis.set("symmetry", "none");
  parametersH2COWatsonNoCoriolis.set("nsweeps", 20);
  parametersH2COWatsonNoCoriolis.set("ngrowsweeps", 3);
  parametersH2COWatsonNoCoriolis.set("nmainsweeps", 3);
  parametersH2COWatsonNoCoriolis.set("chkpfile", "GS.H2CO.chkp.h5");
  parametersH2COWatsonNoCoriolis.set("Nmax", 6);
  parametersH2COWatsonNoCoriolis.set("alpha_initial", 1.0E-8);
  parametersH2COWatsonNoCoriolis.set("alpha_main", 1.0E-16);
  parametersH2COWatsonNoCoriolis.set("alpha_final", 0.);
  maquis::DMRGInterface<double> interfaceOptimizerGS(parametersH2COWatsonNoCoriolis);
  interfaceOptimizerGS.optimize();
  auto energyFromOptimizerGS = interfaceOptimizerGS.energy();
  // Excited-state calculation
  parametersH2COWatsonNoCoriolis.set("n_ortho_states", 1);
  parametersH2COWatsonNoCoriolis.set("ortho_states", "GS.H2CO.chkp.h5");
  parametersH2COWatsonNoCoriolis.set("chkpfile", "ES.H2CO.chkp.h5");
  maquis::DMRGInterface<double> interfaceOptimizerES(parametersH2COWatsonNoCoriolis);
  interfaceOptimizerES.optimize();
  auto energyFromOptimizerES = interfaceOptimizerES.energy();
  // Cleans up stuff
  boost::filesystem::remove_all("GS.H2CO.chkp.h5");
  boost::filesystem::remove_all("ES.H2CO.chkp.h5");
  // FEAST
  auto eMin = energyFromOptimizerGS - (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  auto eMax = energyFromOptimizerGS + (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  parametersH2COWatsonNoCoriolis.set("feast_num_states", 1);
  parametersH2COWatsonNoCoriolis.set("feast_max_iter", 1);
  parametersH2COWatsonNoCoriolis.set("feast_emin", eMin);
  parametersH2COWatsonNoCoriolis.set("feast_emax", eMax);
  parametersH2COWatsonNoCoriolis.set("feast_num_points", 8);
  parametersH2COWatsonNoCoriolis.set("feast_init_type", "basis_state_generic_const");
  parametersH2COWatsonNoCoriolis.set("feast_init_onv", "3,2,1,1,2,0");
  parametersH2COWatsonNoCoriolis.set("feast_overlap_convergence_threshold", 1.0E-5);
  parametersH2COWatsonNoCoriolis.set("feast_energy_convergence_threshold", 1.0E-5);
  //
  parametersH2COWatsonNoCoriolis.set("linsystem_precond", "no");
  parametersH2COWatsonNoCoriolis.set("linsystem_krylov_dim", 50);
  parametersH2COWatsonNoCoriolis.set("linsystem_tol", 1.0E-5);
  parametersH2COWatsonNoCoriolis.set("linsystem_init", "last");
  auto vibrationalLattice = Lattice(parametersH2COWatsonNoCoriolis);
  auto vibrationalModel = ModelType(vibrationalLattice, parametersH2COWatsonNoCoriolis);
  auto vibrationalMPO = make_mpo(vibrationalLattice, vibrationalModel);
  auto feastSimulator = FEASTSimulatorType(parametersH2COWatsonNoCoriolis, vibrationalModel, vibrationalLattice, vibrationalMPO);
  feastSimulator.runFEAST();
  auto feastEnergy = feastSimulator.getEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromOptimizerGS, 1.0E-6);
}

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
  // Setup parameters for the linear system solver.
  parametersBilinearly.set("linsystem_precond", "no");
  parametersBilinearly.set("linsystem_krylov_dim", 50);
  parametersBilinearly.set("linsystem_tol", 1.0E-5);
  parametersBilinearly.set("linsystem_init", "last");
  auto vibrationalLattice = Lattice(parametersBilinearly);
  auto vibrationalModel = ModelType(vibrationalLattice, parametersBilinearly);
  auto vibrationalMPO = make_mpo(vibrationalLattice, vibrationalModel);
  auto feastSimulator = FEASTSimulatorType(parametersBilinearly, vibrationalModel, vibrationalLattice, vibrationalMPO);
  feastSimulator.runFEAST();
  auto feastEnergy = feastSimulator.getEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromOptimizerGS, 1.0E-6);
}

#endif // HAVE_TrivialGroup

#ifdef HAVE_NU1

/** @brief DMRG[FEAST] test for the n-mode Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_FEAST_FAD_Fingerprint, NModeFixture)
{
  parametersFADTwoBodyFingerPrint.set("nsweeps", 20);
  parametersFADTwoBodyFingerPrint.set("init_state", "const");
  parametersFADTwoBodyFingerPrint.set("optimization", "twosite");
  parametersFADTwoBodyFingerPrint.set("twosite_truncation", "heev_truncate");
  parametersFADTwoBodyFingerPrint.set("alpha_initial", 1.0E-8);
  parametersFADTwoBodyFingerPrint.set("alpha_main", 1.0E-15);
  parametersFADTwoBodyFingerPrint.set("alpha_final", 0.);
  parametersFADTwoBodyFingerPrint.set("ngrowsweeps", 2);
  parametersFADTwoBodyFingerPrint.set("nmainsweeps", 2);
  parametersFADTwoBodyFingerPrint.set("max_bond_dimension", 100);
  parametersFADTwoBodyFingerPrint.set("seed", 30061957);
  // Ground-state calculation
  maquis::DMRGInterface<double> interfaceGroundState(parametersFADTwoBodyFingerPrint);
  interfaceGroundState.optimize();
  auto groundStateEnergy = interfaceGroundState.energy();
  // == DMRG[FEAST] ==
  double energyGap = 10.;
  auto eMin = groundStateEnergy - energyGap;
  auto eMax = groundStateEnergy + energyGap;
  // We can affort a smaller number of sweeps
  parametersFADTwoBodyFingerPrint.set("nsweeps", 4);
  // FEAST-specific parameters
  parametersFADTwoBodyFingerPrint.set("feast_num_states", 1);
  parametersFADTwoBodyFingerPrint.set("feast_max_iter", 2);
  parametersFADTwoBodyFingerPrint.set("feast_emin", eMin);
  parametersFADTwoBodyFingerPrint.set("feast_emax", eMax);
  parametersFADTwoBodyFingerPrint.set("feast_num_points", 8);
  parametersFADTwoBodyFingerPrint.set("feast_init_type", "default");
  parametersFADTwoBodyFingerPrint.set("feast_integral_type", "full");
  // Setup parameters for the linear system solver.
  parametersFADTwoBodyFingerPrint.set("linsystem_precond", "no");
  parametersFADTwoBodyFingerPrint.set("linsystem_krylov_dim", 50);
  parametersFADTwoBodyFingerPrint.set("linsystem_tol", 1.0E-12);
  parametersFADTwoBodyFingerPrint.set("linsystem_init", "last");
  maquis::DMRGInterface<std::complex<double>> interfaceFEAST(parametersFADTwoBodyFingerPrint);
  interfaceFEAST.runFEAST();
  auto feastEnergy = maquis::real(interfaceFEAST.energy());
  BOOST_CHECK_SMALL(std::abs((feastEnergy-groundStateEnergy)/groundStateEnergy), 1.0E-10);
}

#endif