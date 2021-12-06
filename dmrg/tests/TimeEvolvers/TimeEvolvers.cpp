/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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

#define BOOST_TEST_MODULE TimeEvolvers

#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include <boost/mpl/list.hpp>
#include <boost/mpl/push_front.hpp>
#include <iostream>
#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/mp_tensors/zerositeproblem.h"
#include "dmrg/mp_tensors/contractions/engine.h"
#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/evolve/TimeEvolvers/timeevolver.h"
#include "Fixtures/TimeEvolversFixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/** Checks that the number of overall blocks of a block_matrix is correct */
BOOST_AUTO_TEST_CASE_TEMPLATE(TimeEvolversAddTime, S, symmetries) {
    // The parameters are not really used
    DmrgParameters p;
    p.set("site_types", "0,0,0,0");
    p.set("L", 4);
    p.set("irrep", 0);
    p.set("nsweeps",2);
    p.set("max_bond_dimension", 100);
    p.set("u1_total_charge1", 1);
    p.set("u1_total_charge2", 1);
    p.set("hamiltonian_units", "Hartree");
    p.set("time_units", "as");
    p.set("propagator_accuracy", 1.0E-10);
    p.set("propagator_maxiter", 10);
    p.set("time_step", 1.);
    p.set("imaginary_time", "no");
    TimeEvolver<cmatrix, S, DmrgParameters> timeEvolver(p);
    auto gotTimeStep = timeEvolver.get_time();
    // Note that here we have an additional 0.5 factor because we propagate by t/2.
    // in the l->r sweep, and then by other t/2. in the r->l sweep.
    auto refValue = 1./(2.*24.18884254);
    BOOST_CHECK_CLOSE(gotTimeStep, refValue, 1e-7);
}

/** Checks that the energy of an MPS is conserved after the propagation of a SiteProblem */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestEnergyConservationSiteproblem, S, symmetries, TestTimeEvolverFixture) {
    // Types declaration
    using BoundaryType = Boundary<typename storage::constrained<cmatrix>::type, S>;
    using contr = contraction::Engine<cmatrix, typename storage::constrained<cmatrix>::type, S>;
    // The parameters are not really used
    DmrgParameters p;
    p.set("integrals_binary", maquis::serialize(integrals));
    p.set("site_types", "0,0,0,0");
    p.set("L", 4);
    p.set("irrep", 0);
    p.set("nsweeps",2);
    p.set("max_bond_dimension", 100);
    p.set("hamiltonian_units", "Hartree");
    p.set("time_units", "as");
    p.set("propagator_accuracy", 1.0E-20);
    p.set("propagator_maxiter", 10);
    p.set("time_step", 1.);
    p.set("imaginary_time", "no");
    p.set("COMPLEX", 1);
    // For SU2U1
    p.set("nelec", 2);
    p.set("spin", 0);
    // For 2U1
    p.set("u1_total_charge1", 1);
    p.set("u1_total_charge2", 1);
    // Prepares the lattice and the model.
    // In both cases, we construct everything as is done in [sim.hpp]
    auto lat = Lattice(p);
    auto model = Model<cmatrix, S>(lat, p);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<cmatrix, S>(lat.size(), *(model.initializer(lat, p)));
    auto latticeSize = mpo.length();
    // Prepares the boundaries
    std::vector<BoundaryType> left, right;
    left.resize(latticeSize+1);
    right.resize(latticeSize+1);
    left[0] = mps.left_boundary();
    for (int iSite = 0; iSite < latticeSize; iSite++)
        left[iSite+1] = contr::overlap_mpo_left_step(mps[iSite], mps[iSite], left[iSite], mpo[iSite]);
    right[latticeSize] = mps.right_boundary();
    for (int iSite = latticeSize-1; iSite >= 0; iSite--)
        right[iSite] = contr::overlap_mpo_right_step(mps[iSite], mps[iSite], right[iSite+1], mpo[iSite]);
    // Calculates the energy from the SiteProblem
    SiteProblem<cmatrix, S> sp0(left[0], right[1], mpo[0]);
    auto initialEnergy = maquis::real(expval(mps, mpo));
    // Prepares the TimeEvolver and forwards propagates it
    auto timeEvolver = TimeEvolver<cmatrix, S, DmrgParameters>(p);
    timeEvolver.evolve(sp0, mps[0], true);
    // Calculates the final energy from the apply method
    auto sigmaVector = sp0.apply(mps[0]);
    auto intermediateEnergy = maquis::real(ietl::dot(mps[0], sigmaVector)) + mpo.getCoreEnergy();
    BOOST_CHECK_CLOSE(initialEnergy, intermediateEnergy, 1e-10);
    // Calculates the final energy as expval
    auto finalEnergy = maquis::real(expval(mps, mpo));
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1e-10);
    // And now with the back-propagation
    timeEvolver.evolve(sp0, mps[0], false);
    finalEnergy = maquis::real(expval(mps, mpo));
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1e-10);
}

/** Checks that the energy of an MPS is conserved after the propagation of a ZeroSiteProblem */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestEnergyConservationZeroSiteproblem, S, symmetries, TestTimeEvolverFixture) {
    // Types declaration
    using BoundaryType = Boundary<typename storage::constrained<cmatrix>::type, S>;
    using contr = contraction::Engine<cmatrix, typename storage::constrained<cmatrix>::type, S>;
    // The parameters are not really used
    DmrgParameters p;
    p.set("integrals_binary", maquis::serialize(integrals));
    p.set("site_types", "0,0,0,0");
    p.set("L", 4);
    p.set("irrep", 0);
    p.set("nsweeps",2);
    p.set("max_bond_dimension", 100);
    p.set("hamiltonian_units", "Hartree");
    p.set("time_units", "as");
    p.set("propagator_accuracy", 1.0E-20);
    p.set("propagator_maxiter", 10);
    p.set("time_step", 1.);
    p.set("imaginary_time", "no");
    p.set("COMPLEX", 1);
    // For SU2U1
    p.set("nelec", 2);
    p.set("spin", 0);
    // For 2U1
    p.set("u1_total_charge1", 1);
    p.set("u1_total_charge2", 1);
    // Prepares the lattice and the model.
    // In both cases, we construct everything as is done in [sim.hpp]
    auto lat = Lattice(p);
    auto model = Model<cmatrix, S>(lat, p);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<cmatrix, S>(lat.size(), *(model.initializer(lat, p)));
    auto latticeSize = mpo.length();
    auto initialEnergy = maquis::real(expval(mps, mpo));
    // Prepares the boundaries
    std::vector<BoundaryType> left, right;
    left.resize(latticeSize+1);
    right.resize(latticeSize+1);
    left[0] = mps.left_boundary();
    for (int iSite = 0; iSite < latticeSize; iSite++)
        left[iSite+1] = contr::overlap_mpo_left_step(mps[iSite], mps[iSite], left[iSite], mpo[iSite]);
    right[latticeSize] = mps.right_boundary();
    for (int iSite = latticeSize-1; iSite >= 0; iSite--)
        right[iSite] = contr::overlap_mpo_right_step(mps[iSite], mps[iSite], right[iSite+1], mpo[iSite]);
    // Calculates the energy from the ZeroSiteProblem
    block_matrix<cmatrix, S> UMat, VMat, QMat, RMat;
    block_matrix<typename alps::numeric::associated_real_diagonal_matrix<cmatrix>::type, S> SMat;
    mps[0].make_left_paired();
    qr(mps[0].data(), QMat, RMat);
    mps[0].replace_left_paired(QMat);
    left[1] = contr::overlap_mpo_left_step(mps[0], mps[0], left[0], mpo[0]);
    ZeroSiteProblem<cmatrix, S> zsp(mpo[0], mpo[1], left[1], right[1]);
    // Prepares the TimeEvolver and forwards propagates it
    auto timeEvolver = TimeEvolver<cmatrix, S, DmrgParameters>(p);
    timeEvolver.evolve(zsp, RMat, true);
    auto sigmaVectorBM = zsp.apply(RMat);
    auto finalEnergy = maquis::real(sigmaVectorBM.scalar_overlap(RMat)) + mpo.getCoreEnergy();
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1e-10);
}

/** Checks that the energy decreases if the imaginary-time propagation is performed */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestEnergyDecreaseSiteproblem, S, symmetries, TestTimeEvolverFixture) {
    // Types declaration
    using BoundaryType = Boundary<typename storage::constrained<cmatrix>::type, S>;
    using contr = contraction::Engine<cmatrix, typename storage::constrained<cmatrix>::type, S>;
    // The parameters are not really used
    DmrgParameters p;
    p.set("integrals_binary", maquis::serialize(integrals));
    p.set("site_types", "0,0,0,0");
    p.set("L", 4);
    p.set("irrep", 0);
    p.set("nsweeps",2);
    p.set("max_bond_dimension", 100);
    p.set("hamiltonian_units", "Hartree");
    p.set("time_units", "as");
    p.set("propagator_accuracy", 1.0E-20);
    p.set("propagator_maxiter", 10);
    p.set("time_step", 1.);
    p.set("imaginary_time", "yes");
    p.set("COMPLEX", 1);
    // For SU2U1
    p.set("nelec", 2);
    p.set("spin", 0);
    // For 2U1
    p.set("u1_total_charge1", 1);
    p.set("u1_total_charge2", 1);
    // Prepares the lattice and the model.
    // In both cases, we construct everything as is done in [sim.hpp]
    auto lat = Lattice(p);
    auto model = Model<cmatrix, S>(lat, p);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<cmatrix, S>(lat.size(), *(model.initializer(lat, p)));
    auto latticeSize = mpo.length();
    mps.normalize_right();
    // Prepares the boundaries
    std::vector<BoundaryType> left, right;
    left.resize(latticeSize+1);
    right.resize(latticeSize+1);
    left[0] = mps.left_boundary();
    for (int iSite = 0; iSite < latticeSize; iSite++)
        left[iSite+1] = contr::overlap_mpo_left_step(mps[iSite], mps[iSite], left[iSite], mpo[iSite]);
    right[latticeSize] = mps.right_boundary();
    for (int iSite = latticeSize-1; iSite >= 0; iSite--)
        right[iSite] = contr::overlap_mpo_right_step(mps[iSite], mps[iSite], right[iSite+1], mpo[iSite]);
    // Calculates the energy from the SiteProblem
    SiteProblem<cmatrix, S> sp0(left[0], right[1], mpo[0]);
    auto initialEnergy = maquis::real(expval(mps, mpo));
    // Prepares the TimeEvolver and forwards propagates it
    auto timeEvolver = TimeEvolver<cmatrix, S, DmrgParameters>(p);
    BOOST_TEST(timeEvolver.isImag());
    timeEvolver.evolve(sp0, mps[0], false);
    mps[0] /= ietl::two_norm(mps[0]);
    // Calculates the final energy as expval
    auto finalEnergy = maquis::real(expval(mps, mpo));
    BOOST_TEST(finalEnergy < initialEnergy);
}