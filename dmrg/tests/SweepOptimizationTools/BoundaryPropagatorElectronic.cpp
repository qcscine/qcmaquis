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

#define BOOST_TEST_MODULE BoundaryPropagatorElectronic

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/SweepBasedAlgorithms/BoundaryPropagator.h"
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
 * @brief Checks that the BoundaryPropagator object works propery.
 *
 * The check is done by verifying that the first element of the right boundary
 * and the last element of the left boundary contain the MPS energy.
 *
 */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestConstructorBoundaryPropagator, S, symmetries, BenzeneFixture)
{
  using BoundaryPropagatorType = BoundaryPropagator<matrix, S, storage::disk>;
  parametersBenzene.set("init_state", "const");
  auto latticeBenzene = Lattice(parametersBenzene);
  auto modelBenzene = Model<matrix, S>(latticeBenzene, parametersBenzene);
  auto mpoBenzene = make_mpo(latticeBenzene, modelBenzene);
  auto mpsBenzeneConst = MPS<matrix, S>(latticeBenzene.size(), *(modelBenzene.initializer(latticeBenzene, parametersBenzene)));
  mpsBenzeneConst.normalize_right();
  auto energy = expval(mpsBenzeneConst, mpoBenzene);
  auto boundaryPropagator = BoundaryPropagatorType(mpsBenzeneConst, mpoBenzene);
  boundaryPropagator.updateRightBoundary(0);
  // Simple checks
  auto lastRightBoundary = boundaryPropagator.getRightBoundary(0);
  BOOST_CHECK_EQUAL(lastRightBoundary.aux_dim(), 1);
  BOOST_CHECK_EQUAL(lastRightBoundary[0].n_blocks(), 1);
  // Energy check
  auto energyFromRightBoundary = lastRightBoundary[0].trace() + mpoBenzene.getCoreEnergy();
  BOOST_CHECK_CLOSE(energyFromRightBoundary, energy, 1.0E-8);
  // Does the same for the left part
  for (int iSite = 1; iSite <= latticeBenzene.size(); iSite++)
    boundaryPropagator.updateLeftBoundary(iSite);
  auto lastLeftBoundary = boundaryPropagator.getLeftBoundary(latticeBenzene.size());
  BOOST_CHECK_EQUAL(lastLeftBoundary.aux_dim(), 1);
  BOOST_CHECK_EQUAL(lastLeftBoundary[0].n_blocks(), 1);
  auto energyFromLeftBoundary = lastLeftBoundary[0].trace() + mpoBenzene.getCoreEnergy();
  BOOST_CHECK_CLOSE(energyFromLeftBoundary, energy, 1.0E-8);
}
