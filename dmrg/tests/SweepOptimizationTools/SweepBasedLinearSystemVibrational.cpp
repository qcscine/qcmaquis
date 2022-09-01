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

#define BOOST_TEST_MODULE SweepBasedLinearSystemVibrational

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/SweepBasedAlgorithms/SweepBasedLinearSystem.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/WatsonFixture.h"

#ifdef DMRG_VIBRATIONAL

/** @brief Checks that the linear system solver via interface works for vibrational problems. */
BOOST_FIXTURE_TEST_CASE(Test_SweepBasedLinearSystemSS_Vibrational_Bilinearlty, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
  parametersBilinearly.set("nsweeps", 5);
  parametersBilinearly.set("max_bond_dimension", 100);
  parametersBilinearly.set("truncation_initial", 1.0E-20);
  parametersBilinearly.set("truncation_main", 1.0E-15);
  auto bilinearlyLattice = Lattice(parametersBilinearly);
  auto bilinearlyModel = Model<matrix, TrivialGroup>(bilinearlyLattice, parametersBilinearly);
  auto bilinearlyMPO = make_mpo(bilinearlyLattice, bilinearlyModel);
  parametersBilinearly.set("init_state", "const");
  auto bilinearlyMPS = MPS<matrix, TrivialGroup>(bilinearlyLattice.size(),
                                                 *(bilinearlyModel.initializer(bilinearlyLattice, parametersBilinearly)));
  bilinearlyMPS.normalize_right();
  // Calculates the energy via the interface
  parametersBilinearly.set("optimization", "singlesite");
  parametersBilinearly.set("symmetry", "none");
  maquis::DMRGInterface<double> interfaceBilinearly(parametersBilinearly);
  interfaceBilinearly.optimize();
  double energyFromInterface = interfaceBilinearly.energy();
  // Parameters that are specific for the solution of the linear system.
  parametersBilinearly.set("linsystem_precond", "yes");
  parametersBilinearly.set("linsystem_init", "mps");
  parametersBilinearly.set("linsystem_max_it", 1);
  parametersBilinearly.set("linsystem_tol", 1.0E-10);
  parametersBilinearly.set("linsystem_krylov_dim", 20);
  parametersBilinearly.set("linsystem_solver", "GMRES");
  // Set the shift of DMRG[IPI] as the energy shifted by -0.1 (note that this Hamiltonian is unitless)
  parametersBilinearly.set("ipi_shift", energyFromInterface-0.1);
  parametersBilinearly.set("ipi_sweep_threshold", 1.0E-5);
  parametersBilinearly.set("ipi_sweeps_per_system", 2);
  parametersBilinearly.set("ipi_iterations", 10);
  maquis::DMRGInterface<double> interfaceBilinearlyIpi(parametersBilinearly);
  interfaceBilinearlyIpi.runInversePowerIteration();
  auto ipiEnergy = interfaceBilinearlyIpi.energy();
  BOOST_CHECK_CLOSE(energyFromInterface, ipiEnergy, 1.0e-7);
#endif // HAVE_TrivialGroup
}

#endif // DMRG_VIBRATIONAL
