/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

 #define BOOST_TEST_MODULE MODEL_EXCITONIC_NMODE

 #include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/vibrational/u1/ExcitonicModel.hpp"
#include "Fixtures/VibronicFixture.h"
#include "maquis_dmrg.h"
#include "dmrg/sim/matrix_types.h"

/** Test for the integral parser for the trivial vibronic Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_ExcitonicNmode, VibronicFixture) {
#ifdef HAVE_U1
  auto lattice = Lattice(parametersNmodeThiopheneOneBody);
  auto integrals = Vibrational::detail::parseIntegralExcitonic<double>(
      parametersNmodeThiopheneOneBody, lattice
  );
  // Checks sizes
  BOOST_CHECK_EQUAL(integrals.first.size(), 84);
  BOOST_CHECK_EQUAL(integrals.second.size(), 84);
  BOOST_CHECK_EQUAL(integrals.first[0].size(), 2);
#endif  // HAVE_U1
}

/** Test for the thiophene dimer based on the excitonic hamiltonian*/
BOOST_FIXTURE_TEST_CASE(Test_Thiophene_Energy_ExcitonicNmode, VibronicFixture) {
#ifdef HAVE_U1
  auto lattice = Lattice(parametersNmodeThiopheneOneBody);
  auto integrals = Vibrational::detail::parseIntegralExcitonic<double>(
      parametersNmodeThiopheneOneBody, lattice
  );
  // Creates the interface for a single excitation
  maquis::DMRGInterface<std::complex<double>> interface_single(
      parametersNmodeThiopheneOneBody
  );
  interface_single.optimize();
  double interface_single_optimizedEnergy = interface_single.energy().real();
  BOOST_CHECK_CLOSE(interface_single_optimizedEnergy, 0.58483639508900598, 1.0E-10);
#endif  // HAVE_U1
}