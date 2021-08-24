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
#include <iostream>
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/optimize/ietl_lanczos_solver.h"
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
#include "dmrg/evolve/siteshifter.h"
#include "dmrg/evolve/perturber.h"
#include "Fixtures/TimeEvolversFixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/** Checks the constructor of a site shifter object */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestConstructorSiteShifter, S, symmetries, TestTimeEvolverFixture) {
    // Types declaration
    using BoundaryType = Boundary<typename storage::constrained<matrix>::type, S>;
    using contr = contraction::Engine<matrix, typename storage::constrained<matrix>::type, S>;
    using TimeEvolverType = TimeEvolver<matrix, S, DmrgParameters>;
    using PerturberType = Perturber<matrix, S>;
    using SiteShifterType = SiteShifter<matrix, S, TimeEvolverType, PerturberType>;
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
    // In both cases, we construct everything as is done in [sim.hpp]
    auto lat = Lattice(p);
    auto model = Model<matrix, S>(lat, p);
    auto mpo = make_mpo(lat, model);
    auto mps = MPS<matrix, S>(lat.size(), *(model.initializer(lat, p)));
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
    // Constructs the SiteShfiter object
    auto pointerToPerturber = std::make_shared<PerturberType>(left, right, mpo, p);
    auto timeEvolver = std::make_shared<TimeEvolverType>(p);
    auto siteShifter = SiteShifter<matrix, S, TimeEvolverType, PerturberType>(mps, timeEvolver, pointerToPerturber, 0, true);
    BOOST_TEST(siteShifter.getSite() == 0);
}