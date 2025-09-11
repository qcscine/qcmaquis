/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE TestEnergyConvergence

#include <filesystem>
#include <boost/mpl/list.hpp>
#include <boost/test/included/unit_test.hpp>
#include "Fixtures/N2Fixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
    TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
    ,
    SU2U1PG
#endif
    >
    symmetries;

/** @brief Test that the 'conv_thresh' works*/
BOOST_FIXTURE_TEST_CASE_TEMPLATE(TestConvThresh, S, symmetries, N2Fixture) {
  parameters.set("max_bond_dimension", 2000);
  parameters.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
  parameters.set("conv_thresh", 1e-9);
  parameters.set("nsweeps", 1);

  // run only one sweep
  maquis::DMRGInterface<double> oneSweepInterface(parameters);
  oneSweepInterface.optimize();
  auto oneSweepEnergy = oneSweepInterface.energy();

  // converge energy
  parameters.set("nsweeps", 100);
  maquis::DMRGInterface<double> convInterface(parameters);
  convInterface.optimize();
  auto convEnergy = convInterface.energy();

  // ensure one sweep does not already converge the MPS
  BOOST_CHECK_GT(std::abs(convEnergy - oneSweepEnergy), 1e-7);
  BOOST_TEST(convEnergy == refEnergy, boost::test_tools::tolerance(1e-8));
}
