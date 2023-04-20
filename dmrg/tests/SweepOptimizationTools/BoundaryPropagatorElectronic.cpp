/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE BoundaryPropagatorElectronic

#include <iostream>
#include <boost/mpl/list.hpp>
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
  parametersBenzene.set("init_type", "const");
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
