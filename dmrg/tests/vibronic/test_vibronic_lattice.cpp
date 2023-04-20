/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <boost/filesystem/operations.hpp>
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/models/lattice/VibronicLattice.hpp"
#include "Fixtures/VibronicFixture.h"

#ifdef HAVE_U1

/** @brief Checks that the constructor for the vibronic lattice work properly */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_Vibronic_Pyrazine, VibronicFixture)
{
    auto lattice = VibronicLattice(parametersVibronic);
    auto size = lattice.size();
    BOOST_CHECK_EQUAL(size, 26);
}

/** @brief Checks that the constructor for the excitonic lattice works properly */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_Vibronic_Excitonic, VibronicFixture)
{
    auto lattice = VibronicLattice(parametersExcitonicAggregate);
    auto size = lattice.size();
    BOOST_CHECK_EQUAL(size, 66);
}

BOOST_FIXTURE_TEST_CASE(Test_Lattice_GetProp_Vibronic_Pyrazine, VibronicFixture)
{
    auto lattice = VibronicLattice(parametersVibronic);
    auto posOfThirdMode = lattice.get_prop<int>("vibindex", 0, 3);
    // Note that the first two positions are taken by the electronic DOF, so the
    // fourth mode (iMode=3) should be in the sixth position (5 starting from 0)
    BOOST_CHECK_EQUAL(posOfThirdMode, 5);
    // The second excited state should instead be in the first position.
    auto posOfSecondExcitedState = lattice.get_prop<int>("eleindex", 0, 1);
    BOOST_CHECK_EQUAL(posOfSecondExcitedState, 1);
}

BOOST_FIXTURE_TEST_CASE(Test_Lattice_GetProp_Excitonic, VibronicFixture)
{
    parametersExcitonicAggregate.set("vibronic_sorting", "intertwined");
    auto lattice = VibronicLattice(parametersExcitonicAggregate);
    auto posOfSecondModeOfSecondMonomer = lattice.get_prop<int>("vibindex", 1, 2);
    // The index 0 is the first excite state of the first monomer, then 1-10 are the next modes,
    // then 11 is the excited states of the second monomer.
    BOOST_CHECK_EQUAL(posOfSecondModeOfSecondMonomer, 14);
    auto posOfExcitedStates = lattice.get_prop<int>("eleindex", 5, 0);
    BOOST_CHECK_EQUAL(posOfExcitedStates, 55);
    // Now changes the sorting
    parametersExcitonicAggregate.set("vibronic_sorting", "firstele");
    auto newLattice = VibronicLattice(parametersExcitonicAggregate);
    posOfSecondModeOfSecondMonomer = newLattice.get_prop<int>("vibindex", 1, 2);
    BOOST_CHECK_EQUAL(posOfSecondModeOfSecondMonomer, 18);
    posOfExcitedStates = newLattice.get_prop<int>("eleindex", 5, 0);
    BOOST_CHECK_EQUAL(posOfExcitedStates, 5);
}

#endif // HAVE_U1
