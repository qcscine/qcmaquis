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

#define BOOST_TEST_MODULE LinearSolver

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
    // Prepares the lattice, the model, and the corresponding MPO
    auto vibrationalLattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto model = Model<matrix, TrivialGroup>(vibrationalLattice, parametersEthyleneWatsonHarmonic);
    auto mpo = make_mpo(vibrationalLattice, model);
    auto latticeSize = mpo.length();
    // Generates the MPS.
    parametersEthyleneWatsonHarmonic.set("init_state", "basis_state_generic");
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
    SiteProblem siteProblem(leftBoundary[0], rightBoundary[1], mpo[0]);
    auto rhs = contraction::site_ortho_boundaries(mpsHF[0], mpsHF[0], orthoLeft[0], orthoRight[1]);
    auto precond = contraction::Engine<matrix, matrix, TrivialGroup>::diagonal_hamiltonian(leftBoundary[0], rightBoundary[1], mpo[0], mpsHF[0]);
    double zShift = 0.;
    auto linearSolver = LinSolver(siteProblem, mpsHF[0], rhs, zShift, parametersEthyleneWatsonHarmonic, precond);
    auto result = linearSolver.res();
    // First check: since H*psi = E*psi, and we set rhs=psi, the solution to the linear system should be
    // 
    auto norm = 1./ietl::two_norm(result.second);
    BOOST_CHECK_CLOSE(norm, referenceHarmonicEnergy, 1.0E-8);
#endif // HAVE_TrivialGroup
}
