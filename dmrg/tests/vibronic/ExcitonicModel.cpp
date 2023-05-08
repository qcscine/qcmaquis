/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE MODEL_EXCITONIC_U1

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/vibrational/u1/ExcitonicModel.hpp"
#include "Fixtures/VibronicFixture.h"
#include "maquis_dmrg.h"
#include "dmrg/sim/matrix_types.h"

/** Test for the integral parser for the trivial vibronic Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_Excitonic, VibronicFixture)
{
#ifdef HAVE_U1
    auto lattice = Lattice(parametersExcitonicAggregate);
    auto integrals = Vibrational::detail::parseIntegralExcitonic<double>(parametersExcitonicAggregate, lattice);
    // Checks sizes
    BOOST_CHECK_EQUAL(integrals.first.size(), 30);
    BOOST_CHECK_EQUAL(integrals.second.size(), 30);
    BOOST_CHECK_EQUAL(integrals.first[0].size(), 2);
#endif // HAVE_U1
}
