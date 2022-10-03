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

#define BOOST_TEST_MODULE ElectronicFEAST

#include <boost/mpl/list.hpp>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/MetaSweepSimulations/FEASTSimulator.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/H2Fixture.h"
#include "Fixtures/LiHFixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/** @brief Test FEAST for electronic calculations (we target the ground and first excited state) */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_FEAST_Electronic_H2, S, symmetries, H2Fixture)
{
  using FEASTSimulatorType = FEASTSimulator<S>;
  using ModelType = Model<cmatrix, S>;
  //
  parametersH2.set("max_bond_dimension", 10);
  parametersH2.set("init_state", "default");
  parametersH2.set("seed", 19893003);
  parametersH2.set("optimization", "twosite");
  parametersH2.set("nsweeps", 5);
  parametersH2.set("chkpfile", "GS.H2.FEAST.chkp.h5");
  parametersH2.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
  maquis::DMRGInterface<double> interfaceOptimizerGS(parametersH2);
  interfaceOptimizerGS.optimize();
  auto energyFromOptimizerGS = interfaceOptimizerGS.energy();
  // Excited-state calculation
  parametersH2.set("n_ortho_states", 1);
  parametersH2.set("ortho_states", "GS.H2.FEAST.chkp.h5");
  parametersH2.set("chkpfile", "ES.H2.FEAST.chkp.h5");
  maquis::DMRGInterface<double> interfaceOptimizerES(parametersH2);
  interfaceOptimizerES.optimize();
  auto energyFromOptimizerES = interfaceOptimizerES.energy();
  // Cleans up stuff
  boost::filesystem::remove_all("GS.H2.FEAST.chkp.h5");
  boost::filesystem::remove_all("ES.H2.FEAST.chkp.h5");
  // FEAST for ground state
  auto eMin = energyFromOptimizerGS - (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  auto eMax = energyFromOptimizerGS + (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  parametersH2.set("feast_num_states", 1);
  parametersH2.set("feast_max_iter", 1);
  parametersH2.set("feast_emin", eMin);
  parametersH2.set("feast_emax", eMax);
  parametersH2.set("feast_num_points", 8);
  parametersH2.set("feast_init_type", "default");
  parametersH2.set("linsystem_precond", "no");
  parametersH2.set("linsystem_krylov_dim", 50);
  parametersH2.set("linsystem_tol", 1.0E-5);
  parametersH2.set("linsystem_init", "last");
  auto H2Lattice = Lattice(parametersH2);
  auto H2Model = ModelType(H2Lattice, parametersH2);
  auto H2MPO = make_mpo(H2Lattice, H2Model);
  auto feastSimulator = FEASTSimulatorType(parametersH2, H2Model, H2Lattice, H2MPO);
  feastSimulator.runFEAST();
  auto feastEnergy = feastSimulator.getEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromOptimizerGS, 1.0E-6);
  // FEAST for excited state
  eMin = energyFromOptimizerES - 0.0001;
  eMax = energyFromOptimizerES + 0.0001;
  parametersH2.set("feast_max_iter", 5);
  parametersH2.set("feast_emin", eMin);
  parametersH2.set("feast_emax", eMax);
  auto feastSimulatorES = FEASTSimulatorType(parametersH2, H2Model, H2Lattice, H2MPO);
  feastSimulatorES.runFEAST();
  feastEnergy = feastSimulatorES.getEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromOptimizerES, 1.0E-6);
}

#ifdef HAVE_SU2U1PG

/**
 * @brief Test FEAST for electronic calculations (we target the first excited state) 
 * Note that, unlike the other test-cases, here we benchmark against DMRG[IPI] and also
 * repeat the FEAST iterations.
 */
BOOST_FIXTURE_TEST_CASE(Test_FEAST_Electronic_LiH, LiHFixture)
{
  using FEASTSimulatorType = FEASTSimulator<SU2U1PG>;
  using ModelType = Model<cmatrix, SU2U1PG>;
  // Generic parameters
  parametersLiH.set("max_bond_dimension", 50);
  parametersLiH.set("init_state", "hf");
  parametersLiH.set("hf_occ", "4,1,1,1");
  parametersLiH.set("optimization", "twosite");
  parametersLiH.set("symmetry", "su2u1pg");
  parametersLiH.set("init_state", "const");
  // IPI-specific parameters
  parametersLiH.set("ipi_sweep_overlap_threshold", 1.0E-5);
  parametersLiH.set("ipi_sweep_energy_threshold", 1.0E-5);
  parametersLiH.set("ipi_sweeps_per_system", 5);
  parametersLiH.set("ipi_iterations", 8);
  // Linear system parameters (note that the same set of parameters is used also for DMRG[FEAST])
  parametersLiH.set("linsystem_precond", "no");
  parametersLiH.set("linsystem_init", "last");
  parametersLiH.set("linsystem_max_it", 1);
  parametersLiH.set("linsystem_tol", 1.0E-10);
  parametersLiH.set("linsystem_krylov_dim", 20);
  // == IPI SIMULATION ==
  // Ground state, reference taken from test2.cpp == -7.90436
  parametersLiH.set("ipi_shift", -8.); 
  maquis::DMRGInterface<double> optimizerIpiGS(parametersLiH);
  optimizerIpiGS.runInversePowerIteration();
  auto energyFromIpiGS = optimizerIpiGS.energy();
  BOOST_CHECK_SMALL(std::abs(energyFromIpiGS - -7.90436), 1.0E-4);
  // First excited-state, reference taken from test2.cpp == -7.77349
  parametersLiH.set("ipi_shift", -7.8);
  maquis::DMRGInterface<double> optimizerIpiES1(parametersLiH);
  optimizerIpiES1.runInversePowerIteration();
  auto energyFromIpiES1 = optimizerIpiES1.energy();
  BOOST_CHECK_SMALL(std::abs(energyFromIpiES1 - -7.77349), 1.0E-4);
  // Second excited state (no reference from test2.cpp, but DMRG data == -7.275314952)
  parametersLiH.set("ipi_shift", -7.3); 
  maquis::DMRGInterface<double> optimizerIpiES2(parametersLiH);
  optimizerIpiES2.runInversePowerIteration();
  auto energyFromIpiES2 = optimizerIpiES2.energy();
  BOOST_CHECK_SMALL(std::abs(energyFromIpiES2 - -7.275314952), 1.0E-4);
  // FEAST for ground state
  parametersLiH.set("nsweeps", 5);
  parametersLiH.set("feast_num_states", 1);
  parametersLiH.set("feast_max_iter", 3);
  parametersLiH.set("feast_num_points", 8);
  parametersLiH.set("feast_init_type", "default");
  parametersLiH.set("feast_overlap_convergence_threshold", 1.0E-5);
  parametersLiH.set("feast_energy_convergence_threshold", 1.0E-6);
  //
  parametersLiH.set("linsystem_precond", "no");
  parametersLiH.set("linsystem_krylov_dim", 50);
  parametersLiH.set("linsystem_tol", 1.0E-5);
  parametersLiH.set("linsystem_init", "last");
  // Ground state 
  auto eMin = energyFromIpiGS - 0.001;
  auto eMax = energyFromIpiGS + 0.001;
  parametersLiH.set("feast_emin", eMin);
  parametersLiH.set("feast_emax", eMax);
  auto LiHLattice = Lattice(parametersLiH);
  auto LiHModel = ModelType(LiHLattice, parametersLiH);
  auto LiHMPO = make_mpo(LiHLattice, LiHModel);
  auto feastSimulator = FEASTSimulatorType(parametersLiH, LiHModel, LiHLattice, LiHMPO);
  feastSimulator.runFEAST();
  auto feastEnergy = feastSimulator.getEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromIpiGS, 1.0E-5);
  // First excited state
  eMin = energyFromIpiES1 - 0.001;
  eMax = energyFromIpiES1 + 0.001;
  parametersLiH.set("feast_emin", eMin);
  parametersLiH.set("feast_emax", eMax);
  auto feastSimulatorES1 = FEASTSimulatorType(parametersLiH, LiHModel, LiHLattice, LiHMPO);
  feastSimulatorES1.runFEAST();
  feastEnergy = feastSimulatorES1.getEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromIpiES1, 1.0E-5);
  // Simultaneous calculation on ground and excited state (so, 2 roots)
  eMin = energyFromIpiGS - 0.001;
  eMax = energyFromIpiES1 + 0.001;
  parametersLiH.set("feast_emin", eMin);
  parametersLiH.set("feast_emax", eMax);
  parametersLiH.set("feast_num_states", 2);
  auto feastSimulatorTwoStates = FEASTSimulatorType(parametersLiH, LiHModel, LiHLattice, LiHMPO);
  feastSimulatorTwoStates.runFEAST();
  auto feastEnergyGS = feastSimulatorTwoStates.getEnergy(0);
  auto feastEnergyES = feastSimulatorTwoStates.getEnergy(1);
  BOOST_CHECK_CLOSE(feastEnergyGS, energyFromIpiGS, 1.0E-5);
  BOOST_CHECK_CLOSE(feastEnergyES, energyFromIpiES1, 1.0E-5);
}

#endif // HAVE_SU2U1PG