/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE LinearSolverVibrational

#include "dmrg/LinearSystem/linsolver.h"
#include "dmrg/mp_tensors/contractions/engine.h"
#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/WatsonFixture.h"

#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include <boost/mpl/list.hpp>
#include <boost/mpl/push_front.hpp>
#include <iostream>

/**
 * @brief Tests the linear solver for a trivial case.
 * 
 * The trivial case is generated as follows:
 * 
 *  - we construct the MPO representation of the harmonic vibrational
 *    Hamiltonian of ethylene.
 *  - we construct the ONV corresponding to the ground state.
 * 
 * This MPS is trivially an eigenvector of the harmonic MPO and, therefore,
 * it should solve the linear system Hx = MPS. The iterative linear system
 * should, therefore, converge in a single iteration.
 */
BOOST_FIXTURE_TEST_CASE(Test_LinearSolver_Trivial, WatsonFixture) {
#ifdef HAVE_TrivialGroup
  // Types declaration
  using BoundaryType = Boundary<typename storage::constrained<matrix>::type, TrivialGroup>;
  using contr = contraction::Engine<matrix, typename storage::constrained<matrix>::type, TrivialGroup>;
  using BlockMatrix = block_matrix<typename storage::constrained<matrix>::type, TrivialGroup>;
  using OrthoContainer = std::vector< BlockMatrix >;
  using SiteProblem = SiteProblem<matrix, TrivialGroup>;
  using LinSolver = LinSolver<matrix, TrivialGroup>;
  parametersEthyleneWatsonHarmonic.set("Nmax", "4");
  // Prepares the lattice, the model, and the corresponding MPO
  auto vibrationalLattice = Lattice(parametersEthyleneWatsonHarmonic);
  auto model = Model<matrix, TrivialGroup>(vibrationalLattice, parametersEthyleneWatsonHarmonic);
  auto mpo = make_mpo(vibrationalLattice, model);
  auto latticeSize = mpo.length();
  // Generates the MPS.
  parametersEthyleneWatsonHarmonic.set("init_type", "basis_state_generic");
  parametersEthyleneWatsonHarmonic.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
  auto mpsHF = MPS<matrix, TrivialGroup>(vibrationalLattice.size(),
                                         *(model.initializer(vibrationalLattice, parametersEthyleneWatsonHarmonic)));
  mpsHF.normalize_right();
  // Prepares the boundaries
  std::vector<BoundaryType> leftBoundary, rightBoundary;
  leftBoundary.resize(latticeSize+1);
  rightBoundary.resize(latticeSize+1);
  leftBoundary[0] = mpsHF.left_boundary();
  for (int iSite = 0; iSite < latticeSize; iSite++)
    leftBoundary[iSite+1] = contr::overlap_mpo_left_step(mpsHF[iSite], mpsHF[iSite], leftBoundary[iSite], mpo[iSite]);
  rightBoundary[latticeSize] = mpsHF.right_boundary();
  for (int iSite = latticeSize-1; iSite >= 0; iSite--)
    rightBoundary[iSite] = contr::overlap_mpo_right_step(mpsHF[iSite], mpsHF[iSite], rightBoundary[iSite+1], mpo[iSite]);
  // Prepares the overlap boundaries
  OrthoContainer orthoLeft(latticeSize+1), orthoRight(latticeSize+1);
  orthoLeft[0] = mpsHF.left_boundary()[0];
  orthoRight[latticeSize] = mpsHF.right_boundary()[0];
  for (int iSite = 0; iSite < latticeSize; iSite++)
    orthoLeft[iSite+1] = contr::overlap_left_step(mpsHF[iSite], mpsHF[iSite], orthoLeft[iSite]);
  for (int iSite = latticeSize-1; iSite >= 0; iSite--)
    orthoRight[iSite] = contr::overlap_right_step(mpsHF[iSite], mpsHF[iSite], orthoRight[iSite+1]);
  // Prepares the SiteProblem object and the corresponding ortho object.
  auto siteProblem = std::make_shared<SiteProblem>(leftBoundary[0], rightBoundary[1], mpo[0]);
  auto rhs = contraction::site_ortho_boundaries(mpsHF[0], mpsHF[0], orthoLeft[0], orthoRight[1]);
  auto precond = std::make_shared<BlockMatrix>(contraction::Engine<matrix, matrix, TrivialGroup>::diagonal_hamiltonian(leftBoundary[0], rightBoundary[1], mpo[0], mpsHF[0]));
  double zShift = 0.;
  auto linearSolver = LinSolver(siteProblem, mpsHF[0], rhs, zShift, parametersEthyleneWatsonHarmonic, precond, true);
  auto result = linearSolver.res();
  // First check: since H*psi = E*psi, and we set rhs=psi, the solution to the linear system should be the inverse of the energy.
  auto norm = 1./ietl::two_norm(std::get<2>(result));
  BOOST_CHECK_CLOSE(norm, referenceHarmonicEnergy, 1.0E-8);
#endif // HAVE_TrivialGroup
}

#ifdef HAVE_TrivialGroup

/**
 * @brief Validates the linear solver for the complex case.
 * The test-case is, in this case, the bilinearly coupled harmonic oscillator.
 */
BOOST_FIXTURE_TEST_CASE(Test_LinearSolver_Trivial_Complex, WatsonFixture) {
  // Types declaration
  using BoundaryType = Boundary<typename storage::constrained<cmatrix>::type, TrivialGroup>;
  using contr = contraction::Engine<cmatrix, typename storage::constrained<cmatrix>::type, TrivialGroup>;
  using BlockMatrix = block_matrix<typename storage::constrained<cmatrix>::type, TrivialGroup>;
  using OrthoContainer = std::vector< BlockMatrix >;
  using SiteProblem = SiteProblem<cmatrix, TrivialGroup>;
  using LinSolver = LinSolver<cmatrix, TrivialGroup>;
  // Prepares the lattice, the model, and the corresponding MPO
  parametersBilinearly.set("max_bond_dimension", 10);
  parametersBilinearly.set("init_type", "default");
  parametersBilinearly.set("seed", 1989);
  auto vibrationalLattice = Lattice(parametersBilinearly);
  auto model = Model<cmatrix, TrivialGroup>(vibrationalLattice, parametersBilinearly);
  auto mpo = make_mpo(vibrationalLattice, model);
  auto latticeSize = mpo.length();
  // Generates the first MPS
  auto initializerPointerFirst = model.initializer(vibrationalLattice, parametersBilinearly);
  auto firstMPS = MPS<cmatrix, TrivialGroup>(latticeSize, *initializerPointerFirst);
  firstMPS.normalize_right();
  // Generates the second MPS (note that the seed is changed in order to ensure that the lhs != rhs)
  parametersBilinearly.set("seed", 1991);
  auto initializerPointerSecond = model.initializer(vibrationalLattice, parametersBilinearly);
  auto secondMPS = MPS<cmatrix, TrivialGroup>(latticeSize, *initializerPointerSecond);
  secondMPS.normalize_right();
  // Prepares the boundaries
  std::vector<BoundaryType> leftBoundary, rightBoundary;
  leftBoundary.resize(latticeSize+1);
  rightBoundary.resize(latticeSize+1);
  leftBoundary[0] = firstMPS.left_boundary();
  for (int iSite = 0; iSite < latticeSize; iSite++)
    leftBoundary[iSite+1] = contr::overlap_mpo_left_step(firstMPS[iSite], firstMPS[iSite], leftBoundary[iSite], mpo[iSite]);
  rightBoundary[latticeSize] = firstMPS.right_boundary();
  for (int iSite = latticeSize-1; iSite >= 0; iSite--)
    rightBoundary[iSite] = contr::overlap_mpo_right_step(firstMPS[iSite], firstMPS[iSite], rightBoundary[iSite+1], mpo[iSite]);
  // Prepares the overlap boundaries
  OrthoContainer orthoLeft(latticeSize+1), orthoRight(latticeSize+1);
  orthoLeft[0] = firstMPS.left_boundary()[0];
  orthoRight[latticeSize] = firstMPS.right_boundary()[0];
  for (int iSite = 0; iSite < latticeSize; iSite++)
    orthoLeft[iSite+1] = contr::overlap_left_step(firstMPS[iSite], secondMPS[iSite], orthoLeft[iSite]);
  for (int iSite = latticeSize-1; iSite >= 0; iSite--)
    orthoRight[iSite] = contr::overlap_right_step(firstMPS[iSite], secondMPS[iSite], orthoRight[iSite+1]);
  auto siteProblem = std::make_shared<SiteProblem>(leftBoundary[0], rightBoundary[1], mpo[0]);
  auto rhs = contraction::site_ortho_boundaries(firstMPS[0], secondMPS[0], orthoLeft[0], orthoRight[1]);
  // Constructs the linear system solver
  parametersBilinearly.set("linsystem_precond", "no");
  parametersBilinearly.set("linsystem_init", "last");
  parametersBilinearly.set("linsystem_max_it", 1);
  parametersBilinearly.set("linsystem_tol", 1.0E-15);
  parametersBilinearly.set("linsystem_krylov_dim", 100);
  parametersBilinearly.set("linsystem_solver", "GMRES");
  std::shared_ptr<BlockMatrix> precond;
  auto zShift = std::complex<double>(392., -521.);
  auto linearSolver = LinSolver(siteProblem, firstMPS[0], rhs, zShift, parametersBilinearly, precond, true);
  auto result = linearSolver.res();
  firstMPS[0] = std::get<2>(result);
  std::complex<double> lhsTerm = expval(firstMPS, mpo)-zShift*overlap(firstMPS, firstMPS);
  std::complex<double> rhsTerm = overlap(firstMPS, secondMPS);
  BOOST_CHECK_CLOSE(std::abs(lhsTerm), std::abs(rhsTerm), 1.0E-10);
}

#endif // HAVE_TrivialGroup