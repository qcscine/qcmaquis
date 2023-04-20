/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

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

#include "test_siteproblem_fixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

BOOST_AUTO_TEST_CASE_TEMPLATE( Test_SiteProblem, S, symmetries)
{
    // Types definition
    using BoundaryType = Boundary<typename storage::constrained<matrix>::type, S>;
    using contr = contraction::Engine<matrix, typename storage::constrained<matrix>::type, S>;
    // The integral file is taken from 
    DmrgParameters p;
    const auto& integrals = TestSiteproblemFixture::integrals;
    p.set("integrals_binary", maquis::serialize(integrals));
    p.set("site_types", "0,0,0,0");
    p.set("L", 4);
    p.set("irrep", 0);
    p.set("nsweeps",2);
    p.set("max_bond_dimension",100);
    // For SU2U1
    p.set("nelec", 2);
    p.set("spin", 0);
    // For 2U1
    p.set("u1_total_charge1", 1);
    p.set("u1_total_charge2", 1);
    // Prepares the lattice and the model.
    // In both cases, we construct everything as is done in [sim.hpp]
    auto lat = Lattice(p);
    auto model = Model<matrix, S>(lat, p);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<matrix, S>(lat.size(), *(model.initializer(lat, p)));
    BOOST_CHECK_EQUAL(mps.size(), mpo.size());
    // Moves the normalization to the second site
    mps.normalize_right();
    auto overlap1 = norm(mps);
    auto overlap2 = mps[0].scalar_overlap(mps[0]);
    BOOST_CHECK_CLOSE(overlap1, overlap2, 1e-10);
    // Calculates the energy
    auto energy = expval(mps, mpo);
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
    SiteProblem<matrix, S> sp0(left[0], right[1], mpo[0]);
    auto sigmaVector = sp0.apply(mps[0]);
    auto energy2 = ietl::dot(sigmaVector, mps[0]) + mpo.getCoreEnergy();
    BOOST_CHECK_CLOSE(energy, energy2, 1e-10);
    energy2 = sp0.get_energy(mps[0]) + mpo.getCoreEnergy();
    BOOST_CHECK_CLOSE(energy, energy2, 1e-10);
    // Next site
    SiteProblem<matrix, S> sp1(left[1], right[2], mpo[1]);
    sigmaVector = sp1.apply(mps[1]);
    energy2 = ietl::dot(sigmaVector, mps[1]) + mpo.getCoreEnergy();
    BOOST_CHECK_CLOSE(energy, energy2, 1e-10);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( Test_ZeroSiteProblem, S, symmetries)
{
    // Types definition
    using BoundaryType = Boundary<typename storage::constrained<matrix>::type, S>;
    using contr = contraction::Engine<matrix, typename storage::constrained<matrix>::type, S>;
    // The integral file is taken from 
    DmrgParameters p;
    const auto& integrals = TestSiteproblemFixture::integrals;
    p.set("integrals_binary", maquis::serialize(integrals));
    p.set("site_types", "0,0,0,0");
    p.set("L", 4);
    p.set("irrep", 0);
    p.set("nsweeps",2);
    p.set("max_bond_dimension", 100);
    // For SU2U1
    p.set("nelec", 2);
    p.set("spin", 0);
    // For 2U1
    p.set("u1_total_charge1", 1);
    p.set("u1_total_charge2", 1);
    // Prepares the lattice and the model.
    // In both cases, we construct everything as is done in [sim.hpp]
    auto lat = Lattice(p);
    auto model = Model<matrix, S>(lat, p);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<matrix, S>(lat.size(), *(model.initializer(lat, p)));
    BOOST_CHECK_EQUAL(mps.size(), mpo.size());
    // Moves the normalization to the second site
    mps.normalize_right();
    auto overlap1 = norm(mps);
    auto overlap2 = mps[0].scalar_overlap(mps[0]);
    BOOST_CHECK_CLOSE(overlap1, overlap2, 1e-10);
    // Calculates the energy
    auto energy = expval(mps, mpo);
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
    // Calculates the energy from the ZeroSiteProblem
    block_matrix<matrix, S> UMat, VMat, QMat, RMat;
    block_matrix<typename alps::numeric::associated_real_diagonal_matrix<matrix>::type, S> SMat;
    mps[0].make_left_paired();
    qr(mps[0].data(), QMat, RMat);
    mps[0].replace_left_paired(QMat);
    left[1] = contr::overlap_mpo_left_step(mps[0], mps[0], left[0], mpo[0]);
    gemm(SMat, VMat, UMat);
    ZeroSiteProblem<matrix, S> zsp(mpo[0], mpo[1], left[1], right[1]);
    auto sigmaVectorBM = zsp.apply(RMat);
    auto energy3 = sigmaVectorBM.scalar_overlap(RMat) + mpo.getCoreEnergy();
    BOOST_CHECK_CLOSE(energy, energy3, 1e-10);
}
