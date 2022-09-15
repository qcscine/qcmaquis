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
#include "dmrg/sim/matrix_types.h"
#include "dmrg/MetaSweepSimulations/FEASTSimulator.h"
#include "utils/fpcomparison.h"

BOOST_FIXTURE_TEST_CASE(Test_FEAST_Electronic, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
  using FEASTSimulatorType = FEASTSimulator<matrix, TrivialGroup>;
  using ModelType = Model<matrix, TrivialGroup>;
  //
  parametersH2COWatson.set("max_bond_dimension", 100);
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
  // FEAST
  auto eMin = energyFromOptimizerGS - (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  auto eMax = energyFromOptimizerGS + (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  parametersH2COWatson.set("feast_num_states", 1);
  parametersH2COWatson.set("feast_max_iter", 5);
  parametersH2COWatson.set("feast_emin", eMin);
  parametersH2COWatson.set("feast_emax", eMax);
  parametersH2COWatson.set("feast_num_points", 8);
  auto vibrationalLattice = Lattice(parametersH2COWatson);
  auto vibrationalModel = ModelType(vibrationalLattice, parametersH2COWatson);
  auto feastSimulator = FEASTSimulatorType(parametersH2COWatson, vibrationalModel, vibrationalLattice);
  // Cleans up stuff
  boost::filesystem::remove_all("GS.H2CO.chkp.h5");
  boost::filesystem::remove_all("ES.H2CO.chkp.h5");
#endif // HAVE_TrivialGroup
}