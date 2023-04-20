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
#include "dmrg/models/lattice/WatsonLattice.hpp"
#include "Fixtures/WatsonFixture.h"

#ifdef HAVE_TrivialGroup

/** @brief Checks that the constructor for the Watson-based Hamiltonian works properly */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_Watson_Ethylene, WatsonFixture)
{
    auto lattice = WatsonLattice(parametersEthyleneWatson);
    auto size = lattice.size();
    BOOST_CHECK_EQUAL(size, 12);
}

/** @brief Checks that the constructor for the Watson-based Hamiltonian works properly */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_MaxVertexType_Watson_Ethylene, WatsonFixture)
{
    auto lattice = WatsonLattice(parametersEthyleneWatson);
    auto size = lattice.size();
    auto numSiteTypes = lattice.getMaxType();
    BOOST_CHECK_EQUAL(size, numSiteTypes);
}

/** @brief Checks that the property getter works properly */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_PropertyGetter_Watson_Ethylene, WatsonFixture)
{
    auto lattice = WatsonLattice(parametersEthyleneWatson);
    for (int iSite = 0; iSite < lattice.size(); iSite++) {
        auto siteType = lattice.get_prop<int>("type", iSite);
        BOOST_CHECK_EQUAL(siteType, iSite);
    }
}

#endif // HAVE_TrivialGroup
