/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE SweepBasedLinearSystemElectronic

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/SweepBasedAlgorithms/SweepBasedLinearSystem.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/BenzeneFixture.h"

/**
 * @brief Checks that the linear system solver works for electronic problems.
 */
BOOST_FIXTURE_TEST_CASE(Test_SweepBasedLinearSystemSS_Electronic_Benzene, BenzeneFixture)
{
#ifdef HAVE_TwoU1PG
  using SweepBasedLinearSolverSS = SweepBasedLinearSystem<matrix, TwoU1PG, storage::disk, SweepOptimizationType::SingleSite>;
  parametersBenzene.set("nsweeps", 10);
  parametersBenzene.set("max_bond_dimension", 100);
  parametersBenzene.set("alpha_initial", 1.0E-8);
  parametersBenzene.set("alpha_main", 1.0E-15);
  parametersBenzene.set("alpha_final", 0.);
  auto benzeneLattice = Lattice(parametersBenzene);
  auto benzeneModel = Model<matrix, TwoU1PG>(benzeneLattice, parametersBenzene);
  auto benzeneMPO = make_mpo(benzeneLattice, benzeneModel);
  parametersBenzene.set("init_type", "hf");
  parametersBenzene.set("hf_occ", "4,4,4,1,1,1");
  auto hfBenzeneMPS = MPS<matrix, TwoU1PG>(benzeneLattice.size(), *(benzeneModel.initializer(benzeneLattice, parametersBenzene)));
  hfBenzeneMPS.normalize_right();
  // Calculates the energy via the interface
  parametersBenzene.set("optimization", "twosite");
  parametersBenzene.set("symmetry", "2u1pg");
  maquis::DMRGInterface<double> interfaceBenzene(parametersBenzene);
  interfaceBenzene.optimize();
  double energyFromInterface = interfaceBenzene.energy();
  // Parameters that are specific for the solution of the linear system.
  parametersBenzene.set("linsystem_precond", "no");
  parametersBenzene.set("linsystem_init", "mps");
  parametersBenzene.set("linsystem_max_it", 1);
  parametersBenzene.set("linsystem_tol", 1.0E-10);
  parametersBenzene.set("linsystem_krylov_dim", 100);
  parametersBenzene.set("linsystem_solver", "GMRES");
  parametersBenzene.set("linsystem_exact_error", "yes");
  // Set the shift of DMRG[IPI] as the energy - 1 Hartree
  parametersBenzene.set("nsweeps", 3);
  parametersBenzene.set("ipi_shift", energyFromInterface-0.1);
  std::vector<double> energyFromIPI;
  // Does the IPI iteration "by hand"
  int nIPI = 10;
  for (int iSweep = 0; iSweep < nIPI; iSweep++) {
    auto linearSolver = SweepBasedLinearSolverSS(hfBenzeneMPS, benzeneMPO, parametersBenzene, benzeneModel, benzeneLattice, false);
    linearSolver.runSweepSimulation();
    energyFromIPI.push_back(linearSolver.template getSpecificResult<double>("Energy"));
  }
  BOOST_CHECK_CLOSE(energyFromInterface, energyFromIPI[nIPI-1], 1.0e-7);
#endif // HAVE_TwoU1PG
}

#ifdef HAVE_SU2U1PG

/** @brief Same as above, but 1) for SU2U1 2) via the interface and 3) with the two-site optimizer */
BOOST_FIXTURE_TEST_CASE(Test_SweepBasedLinearSystemTS_Interface_Electronic_Benzene, BenzeneFixture)
{
  // Prepares the input parameters
  parametersBenzene.set("max_bond_dimension", 100);
  parametersBenzene.set("optimization", "twosite");
  parametersBenzene.set("symmetry", "su2u1pg");
  parametersBenzene.set("init_type", "const");
  // Optimization-specific parameters
  parametersBenzene.set("nsweeps", 10);
  parametersBenzene.set("chkpfile", "GS.Benzene.chkp.h5");
  maquis::DMRGInterface<double> interfaceOptimizerGS(parametersBenzene);
  interfaceOptimizerGS.optimize();
  auto energyFromOptimizerGS = interfaceOptimizerGS.energy();
  // Excited-state calculation
  parametersBenzene.set("n_ortho_states", 1);
  parametersBenzene.set("ortho_states", "GS.Benzene.chkp.h5");
  parametersBenzene.set("chkpfile", "ES.Benzene.chkp.h5");
  maquis::DMRGInterface<double> interfaceOptimizerES(parametersBenzene);
  interfaceOptimizerES.optimize();
  auto energyFromOptimizerES = interfaceOptimizerES.energy();
  // IPI-specific parameters
  double shiftGS = energyFromOptimizerGS-(energyFromOptimizerES-energyFromOptimizerGS)/10.;
  parametersBenzene.set("ipi_shift", shiftGS);
  parametersBenzene.set("ipi_sweep_threshold", 1.0E-5);
  parametersBenzene.set("ipi_sweeps_per_system", 2);
  parametersBenzene.set("ipi_iterations", 10);
  parametersBenzene.set("linsystem_precond", "no");
  parametersBenzene.set("linsystem_init", "mps");
  parametersBenzene.set("linsystem_max_it", 1);
  parametersBenzene.set("linsystem_tol", 1.0E-10);
  parametersBenzene.set("linsystem_krylov_dim", 30);
  parametersBenzene.set("linsystem_solver", "GMRES");
  parametersBenzene.set("chkpfile", "GS.IPI.Benzene.chkp.h5");
  maquis::DMRGInterface<double> interfaceGroundStateIPI(parametersBenzene);
  interfaceGroundStateIPI.runInversePowerIteration();
  auto energyGroundStateIPI = interfaceGroundStateIPI.energy();
  BOOST_CHECK_CLOSE(energyFromOptimizerGS, energyGroundStateIPI, 1.0E-7);
  //
  double shiftES = energyFromOptimizerES-(energyFromOptimizerES-energyFromOptimizerGS)/10.;
  parametersBenzene.set("ipi_shift", shiftES);
  parametersBenzene.set("chkpfile", "ES.IPI.Benzene.chkp.h5");
  maquis::DMRGInterface<double> interfaceExcitedStateIPI(parametersBenzene);
  interfaceExcitedStateIPI.runInversePowerIteration();
  auto energyExcitedStateIPI = interfaceExcitedStateIPI.energy();
  BOOST_CHECK_CLOSE(energyFromOptimizerES, energyExcitedStateIPI, 1.0E-7);
  // Cleans up stuff
  boost::filesystem::remove_all("GS.Benzene.chkp.h5");
  boost::filesystem::remove_all("ES.Benzene.chkp.h5");
  boost::filesystem::remove_all("GS.IPI.Benzene.chkp.h5");
  boost::filesystem::remove_all("ES.IPI.Benzene.chkp.h5");
}


#endif // HAVE_SU2U1PG
