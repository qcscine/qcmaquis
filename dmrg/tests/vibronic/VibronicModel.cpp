/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE MODEL_VIBRONIC_U1

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/vibrational/u1/VibronicModel.hpp"
#include "Fixtures/VibronicFixture.h"
#include "maquis_dmrg.h"
#include "dmrg/sim/matrix_types.h"

/** Test for the integral parser for the trivial vibronic Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_Vibronic, VibronicFixture)
{
#ifdef HAVE_U1
    auto lattice = Lattice(parametersFakeVibronic);
    auto integrals = Vibrational::detail::parseIntegralVibronic<double>(parametersFakeVibronic, lattice);
    // Checks sizes
    BOOST_CHECK_EQUAL(integrals.first.size(), 6);
    BOOST_CHECK_EQUAL(integrals.second.size(), 6);
    BOOST_CHECK_EQUAL(integrals.first[0].size(), 4);
#endif // HAVE_U1
}

#ifdef HAVE_U1

/** Tests the [create_terms] method for the fake vibronic model */
BOOST_FIXTURE_TEST_CASE(Test_VibronicModel_Create_Terms, VibronicFixture)
{
    auto lattice = Lattice(parametersFakeVibronic);
    auto fakeVibronicModel = VibronicModel<tmatrix<double>>(lattice, parametersFakeVibronic);
    auto sizeBefore = fakeVibronicModel.hamiltonian_terms().size();
    BOOST_CHECK_EQUAL(sizeBefore, 0);
    fakeVibronicModel.create_terms();
    auto sizeAfter = fakeVibronicModel.hamiltonian_terms().size();
    BOOST_CHECK_EQUAL(sizeAfter, 6);
}

/** Checks consistency for the physical dimensions for the "trivial" vibronic model */
BOOST_FIXTURE_TEST_CASE(Test_Model_PhysDim_FakeVibronic, VibronicFixture)
{
    auto lattice = Lattice(parametersFakeVibronic);
    auto fakeVibronicModel = VibronicModel<tmatrix<double>>(lattice, parametersFakeVibronic);
    BOOST_CHECK_EQUAL(fakeVibronicModel.phys_dim(0).sum_of_sizes(), 6);
    BOOST_CHECK_EQUAL(fakeVibronicModel.phys_dim(1).sum_of_sizes(), 2);
}

/** Checks consistency for the overall QN for the fake vibronic model */
BOOST_FIXTURE_TEST_CASE(Test_Model_TotalQN_FakeVibronic, VibronicFixture)
{
    auto lattice = Lattice(parametersFakeVibronic);
    auto fakeVibronicModel = VibronicModel<tmatrix<double>>(lattice, parametersFakeVibronic);
    auto totalQN = fakeVibronicModel.total_quantum_numbers(parametersFakeVibronic);
    BOOST_CHECK_EQUAL(totalQN, 1);
}

/** Simple check on tags */
BOOST_FIXTURE_TEST_CASE(Test_Model_Tag_SimpleCheck_FakeVibronic, VibronicFixture)
{
    auto lattice = Lattice(parametersFakeVibronic);
    auto fakeVibronicModel = VibronicModel<tmatrix<double>>(lattice, parametersFakeVibronic);
    auto identityTag = fakeVibronicModel.identity_matrix_tag(0);
    auto fillingTag = fakeVibronicModel.filling_matrix_tag(0);
    BOOST_CHECK(identityTag == fillingTag);
    identityTag = fakeVibronicModel.identity_matrix_tag(1);
    fillingTag = fakeVibronicModel.filling_matrix_tag(1);
    BOOST_CHECK(identityTag == fillingTag);
}

/** Check on symbolic operator getter */
BOOST_FIXTURE_TEST_CASE(Test_Model_Symbolic_Operator_FakeVibronic, VibronicFixture)
{
    auto lattice = Lattice(parametersFakeVibronic);
    auto fakeVibronicModel = VibronicModel<tmatrix<double>>(lattice, parametersFakeVibronic);
    BOOST_CHECK(fakeVibronicModel.filling_matrix_tag(0) == fakeVibronicModel.get_operator_tag("fill", 0));
    BOOST_CHECK(fakeVibronicModel.filling_matrix_tag(1) == fakeVibronicModel.get_operator_tag("fill", 1));
}

#endif // HAVE_U1
