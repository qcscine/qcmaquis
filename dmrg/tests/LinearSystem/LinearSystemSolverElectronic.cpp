/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE LinearSolverElectronic

#include "dmrg/LinearSystem/LinSystemTraitsClass.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/SweepBasedAlgorithms/SweepBasedLinearSystem.h"
#include "Fixtures/BenzeneFixture.h"
#include "utils/fpcomparison.h"
#include <boost/test/included/unit_test.hpp>

/**
 * @brief Validates the linear solver for electronic problems.
 * The test-case is, in this case, Benzene.
 */
BOOST_FIXTURE_TEST_CASE(Test_LinearSolver_Electronic, BenzeneFixture) {
#ifdef HAVE_TwoU1PG
  using SimulatorType = SweepBasedLinearSystem<cmatrix, TwoU1PG, storage::disk, SweepOptimizationType::TwoSite>;
  parametersBenzene.set("max_bond_dimension", 100);
  parametersBenzene.set("init_type", "hf");
  parametersBenzene.set("hf_occ", "4,4,4,1,1,1");
  parametersBenzene.set("nsweeps", 10);
  parametersBenzene.set("truncation_initial", 1.0E-30);
  parametersBenzene.set("truncation_final", 1.0E-30);
  parametersBenzene.set("linsystem_init", "last");
  parametersBenzene.set("linsystem_max_it", 1);
  parametersBenzene.set("linsystem_tol", 1.0E-10);
  parametersBenzene.set("linsystem_krylov_dim", 10);
  parametersBenzene.set("linsystem_exact_error", "yes");
  // Defines variables required by the simulator
  auto benzeneLattice = Lattice(parametersBenzene);
  auto benzeneModel = Model<cmatrix, TwoU1PG>(benzeneLattice, parametersBenzene);
  auto benzeneMPO = make_mpo(benzeneLattice, benzeneModel);
  auto initializer = benzeneModel.initializer(benzeneLattice, parametersBenzene);
  auto rhsMps = MPS<cmatrix, TwoU1PG>(benzeneLattice.size(), *initializer);
  auto lhsMps = rhsMps;
  // Generates the simulator and resets the shift to a meaningless complex number.
  auto zShift = std::complex<double>(392., -521.);
  // std::cout << "ALB BEFORE " << overlap(lhsMps, rhsMps) << std::endl;
  auto simulator = SimulatorType(lhsMps, benzeneMPO, parametersBenzene, benzeneModel, benzeneLattice, true);
  // std::cout << "ALB AFTER " << overlap(lhsMps, rhsMps) << std::endl;
  simulator.setShift(zShift);
  simulator.runSweepSimulation();
  // std::cout << "ALB2" << std::endl;
  // std::cout << lhsMps[0] << std::endl;
  // std::cout << rhsMps[0] << std::endl;
  // std::cout << zShift << std::endl;
  auto error = LinSystemTraitClass<cmatrix, TwoU1PG>::calculateError(lhsMps, rhsMps, benzeneMPO, zShift, benzeneModel, benzeneLattice,
                                                                     benzeneModel.total_quantum_numbers(parametersBenzene), 100);
  BOOST_CHECK_SMALL(std::abs(error), 1.0E-10);
}

#endif // HAVE_TrivialGroup