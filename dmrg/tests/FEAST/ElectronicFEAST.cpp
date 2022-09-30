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

#include <boost/test/included/unit_test.hpp>
#include "dmrg/MetaSweepSimulations/FEASTSimulator.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/H2Fixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/** @brief Test FEAST for electronic calculations (we target the first excited state) */
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