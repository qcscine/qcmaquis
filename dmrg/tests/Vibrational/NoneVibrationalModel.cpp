// Test file for Nmodecopact

#define BOOST_TEST_MODULE MODEL_VIBRATIONAL_PAIRED

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/vibrational/none/NModeModelPairedOperators.hpp"
#include "Fixtures/NModeFixture.h"
#include "maquis_dmrg.h"
#include "dmrg/sim/matrix_types.h"

/** Test for the integral parser with the one-body Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_OneBody_Paired, NModeFixture) {
#ifdef HAVE_TrivialGroup
  auto integrals = Vibrational::detail::NModeIntegralParser<double>(
      parametersFADOneBodyPaired, lattice_factory(parametersFADOneBodyPaired)
  );
  // Checks sizes
  BOOST_CHECK_EQUAL(integrals.first.size(), 78);
  BOOST_CHECK_EQUAL(integrals.second.size(), 78);
  BOOST_CHECK_EQUAL(integrals.first[0].size(), 12);
  // Check numeric values
  BOOST_CHECK_CLOSE(-2.359242429009664e+03, integrals.second[0], 1e-12);
  BOOST_CHECK_EQUAL(1, integrals.first[0][0]);
  BOOST_CHECK_EQUAL(0, integrals.first[0][1]);
  BOOST_CHECK_EQUAL(1, integrals.first[0][2]);
  BOOST_CHECK_EQUAL(0, integrals.first[0][3]);
  BOOST_CHECK_EQUAL(-1, integrals.first[0][4]);
  BOOST_CHECK_EQUAL(-1, integrals.first[0][5]);
  BOOST_CHECK_EQUAL(-1, integrals.first[0][6]);
  BOOST_CHECK_EQUAL(-1, integrals.first[0][7]);
  BOOST_CHECK_EQUAL(-1, integrals.first[0][8]);
  BOOST_CHECK_EQUAL(-1, integrals.first[0][9]);
  BOOST_CHECK_EQUAL(-1, integrals.first[0][10]);
  BOOST_CHECK_EQUAL(-1, integrals.first[0][11]);
#endif  // TrivialGroup
}

#ifdef HAVE_TrivialGroup
/** Test for the integral parser with the two-body Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_TwoBody_Paired, NModeFixture) {
  auto integrals = Vibrational::detail::NModeIntegralParser<double>(
      parametersFADTwoBodyPaired, lattice_factory(parametersFADTwoBodyPaired)
  );
  // Checks sizes
  BOOST_CHECK_EQUAL(integrals.first.size(), 4845);
  BOOST_CHECK_EQUAL(integrals.second.size(), 4845);
  // Check numeric values
  int idx = 4844;
  BOOST_CHECK_CLOSE(2.904591355839877e+04, integrals.second[idx], 1e-12);
  BOOST_CHECK_EQUAL(1, integrals.first[idx][0]);
  BOOST_CHECK_EQUAL(10, integrals.first[idx][1]);
  BOOST_CHECK_EQUAL(1, integrals.first[idx][2]);
  BOOST_CHECK_EQUAL(10, integrals.first[idx][3]);
  BOOST_CHECK_EQUAL(2, integrals.first[idx][4]);
  BOOST_CHECK_EQUAL(10, integrals.first[idx][5]);
  BOOST_CHECK_EQUAL(2, integrals.first[idx][6]);
  BOOST_CHECK_EQUAL(10, integrals.first[idx][7]);
  BOOST_CHECK_EQUAL(-1, integrals.first[idx][8]);
  BOOST_CHECK_EQUAL(-1, integrals.first[idx][9]);
  BOOST_CHECK_EQUAL(-1, integrals.first[idx][10]);
  BOOST_CHECK_EQUAL(-1, integrals.first[idx][11]);
}

/** Test for the threshold functionality of the integral parser */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_Threshold_Paired, NModeFixture) {
  auto integrals = Vibrational::detail::NModeIntegralParser<double>(
      parametersFADOneBodyPaired, lattice_factory(parametersFADOneBodyPaired)
  );
  BOOST_CHECK_EQUAL(integrals.first.size(), 78);
  BOOST_CHECK_EQUAL(integrals.second.size(), 78);
  parametersFADOneBodyPaired.set("integral_cutoff", 1.0E-6);
  integrals = Vibrational::detail::NModeIntegralParser<double>(
      parametersFADOneBodyPaired, lattice_factory(parametersFADOneBodyPaired)
  );
  BOOST_CHECK_EQUAL(integrals.first.size(), 39);
  BOOST_CHECK_EQUAL(integrals.second.size(), 39);
}

/** Tests the [create_terms] method */
BOOST_FIXTURE_TEST_CASE(Test_Model_Create_Terms_Paired, NModeFixture) {
  auto lattice = Lattice(parametersFADOneBodyPaired);
  auto nModeModel =
      NModeModelPaired<tmatrix<double>>(lattice, parametersFADOneBodyPaired, false);
  auto sizeBefore = nModeModel.hamiltonian_terms().size();
  BOOST_CHECK_EQUAL(sizeBefore, 0);
  nModeModel.create_terms();
  auto sizeAfter = nModeModel.hamiltonian_terms().size();
  BOOST_CHECK_EQUAL(sizeAfter, 78);
}

/** Checks consistency for the physical dimensions for a 1-mode system */
BOOST_FIXTURE_TEST_CASE(Test_Model_PhysDim_OneMode_Paired, NModeFixture) {
  auto lattice = Lattice(parametersFADOneBodyPaired);
  auto nModeModel =
      NModeModelPaired<tmatrix<double>>(lattice, parametersFADOneBodyPaired, false);
  const auto& physicalDimensions0 = nModeModel.phys_dim(0);
  BOOST_CHECK_EQUAL(physicalDimensions0.sum_of_sizes(), 39);
}


#endif  // HAVE_TrivialGroup
