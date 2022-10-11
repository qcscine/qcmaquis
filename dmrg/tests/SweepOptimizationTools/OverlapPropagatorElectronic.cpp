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

#define BOOST_TEST_MODULE OverlapPropagatorElectronic

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/SweepBasedAlgorithms/OverlapPropagator.h"
#include "Fixtures/BenzeneFixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/**
 * @brief Checks that the OverlapPropagator object works propery.
 * 
 * The check is done by verifying that, if the object is constructed from two orthogonal MPSs,
 * the overlap between the reference MPS for a given site and the MPSTensor of the other MPS
 * is zero.
 */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestOverlapPropagatorElectronic, S, symmetries, BenzeneFixture)
{
  using OverlapPropagatorType = OverlapPropagator<matrix, S, storage::disk>;
  using MPSType = MPS<matrix, S>;
  parametersBenzene.set("init_state", "hf");
  auto latticeBenzene = Lattice(parametersBenzene);
  auto modelBenzene = Model<matrix, S>(latticeBenzene, parametersBenzene);
  auto mpoBenzene = make_mpo(latticeBenzene, modelBenzene);
  parametersBenzene.set("hf_occ", "4,4,4,1,1,1");
  auto mpsGS = MPS<matrix, S>(latticeBenzene.size(), *(modelBenzene.initializer(latticeBenzene, parametersBenzene)));
  parametersBenzene.set("hf_occ", "4,4,1,4,1,1");
  auto mpsExc1 = MPS<matrix, S>(latticeBenzene.size(), *(modelBenzene.initializer(latticeBenzene, parametersBenzene)));
  parametersBenzene.set("hf_occ", "4,1,4,1,4,1");
  auto mpsExc2 = MPS<matrix, S>(latticeBenzene.size(), *(modelBenzene.initializer(latticeBenzene, parametersBenzene)));
  mpsGS.normalize_right();
  mpsExc1.normalize_right();
  mpsExc2.normalize_right();
  std::vector<MPSType> excitedStateMPSs;
  excitedStateMPSs.push_back(mpsExc1);
  excitedStateMPSs.push_back(mpsExc2);
  auto overlapPropagator = OverlapPropagatorType(mpsGS, excitedStateMPSs, 0);
  // Checks orthogonality
  auto ortho1 = overlapPropagator.template getOrthogonalVector<SweepOptimizationType::SingleSite>(0, 0, 1);
  auto overlap1 = ietl::dot(ortho1, mpsGS[0]);
  BOOST_CHECK_SMALL(overlap1, 1.0E-10);
  auto ortho2 = overlapPropagator.template getOrthogonalVector<SweepOptimizationType::SingleSite>(1, 0, 1);
  auto overlap2 = ietl::dot(ortho1, mpsGS[0]);
  BOOST_CHECK_SMALL(overlap2, 1.0E-10);
}
