/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics,ETH Zurich
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
 *
 * This software is part of the ALPS Applications,published under the ALPS
 * Application License; you can use,redistribute it and/or modify it under
 * the terms of the license,either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not,the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS",WITHOUT WARRANTY OF ANY KIND,EXPRESS OR
 * IMPLIED,INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE,TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY,WHETHER IN CONTRACT,TORT OR OTHERWISE,
 * ARISING FROM,OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <boost/filesystem/operations.hpp>
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/models/lattice/NModeLattice.hpp"
#include "dmrg/models/lattice/WatsonLattice.hpp"
#include "Fixtures/NModeFixture.h"
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

#ifdef HAVE_NU1

/** @brief Checks the size of the lattice */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_2ModeSystem, NModeFixture)
{
    auto lattice = NModeLattice(parametersTwoMode);
    auto size = lattice.size();
    BOOST_CHECK_EQUAL(size, 24);
}

/** @brief Checks the size of the lattice */
BOOST_FIXTURE_TEST_CASE(Test_Site_Types_2ModeSystem, NModeFixture)
{
    auto lattice = NModeLattice(parametersTwoMode);
    auto typeOfSites = lattice.getMaxType();
    BOOST_CHECK_EQUAL(typeOfSites, 2);
}

/** @brief Checks the size of the lattice for the 4-mode input */
BOOST_FIXTURE_TEST_CASE(Test_Site_Types_4ModeSystem, NModeFixture)
{
    auto lattice = NModeLattice(parametersFourMode);
    auto typeOfSites = lattice.getMaxType();
    BOOST_CHECK_EQUAL(typeOfSites, 4);
}

/** @brief Checks the partition of the lattice for the 4-mode input */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Partition_4ModeSystem, NModeFixture)
{
    auto lattice = NModeLattice(parametersFourMode);
    int posOfFirstType = lattice.get_prop<int>("sublatticePos", 0);
    BOOST_CHECK_EQUAL(posOfFirstType, 0);
    posOfFirstType = lattice.get_prop<int>("sublatticePos", 2);
    BOOST_CHECK_EQUAL(posOfFirstType, 7);
}

/** @brief Checks that the lattice is contruscted correctly from parameters data */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_From_Parameters, NModeFixture)
{
    auto lattice = lattice_factory(parametersFADOneBody);
    BOOST_CHECK_EQUAL(lattice->size(), 39);
}

/** @brief Validates the implementation of a generic sorting for an n-mode Lattice */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_NMode_Arbitrary_Sorting, NModeFixture) {
    DmrgParameters parametersArbitrarySorting;
    std::string order = "22,27,24,49,48,52,50,51,56,55,86,84,85,57,54,53,47,28,23,46,26,25,20,19,17,18,16,15,14,13,12,10,11,9,8,7,6,5,4,3,2,1,0,79,77,76,74,72,61,60,62,63,70,64,69,";
    order += "65,67,68,66,59,71,58,73,29,30,33,34,32,31,35,36,37,75,38,81,39,40,41,78,82,42,80,43,83,44,45,21";
    parametersArbitrarySorting.set("modals_order", order);
    parametersArbitrarySorting.set("symmetry", "nu1");
    parametersArbitrarySorting.set("model_library", "coded");
    parametersArbitrarySorting.set("lattice_library", "coded");
    parametersArbitrarySorting.set("LATTICE", "nmode lattice");
    parametersArbitrarySorting.set("L", 87);
    parametersArbitrarySorting.set("nmode_num_modes", 3);
    parametersArbitrarySorting.set("nmode_num_basis", "29,29,29");
    parametersArbitrarySorting.set("MODEL", "nmode");
    auto lattice = NModeLattice(parametersArbitrarySorting);
    // Checks that all elements appear only once
    std::set<int> visitedElements;
    for (int iMode = 0; iMode < 3; iMode++) {
        for (int iModal = 0; iModal < 29; iModal++) {
            auto tmpPos = lattice.get_prop<int>("absolutePositionInLattice", iMode, iModal);
            BOOST_CHECK(visitedElements.find(tmpPos) == visitedElements.end());
            visitedElements.insert(tmpPos);
        }
    }
    BOOST_CHECK_EQUAL(visitedElements.size(), 87);
}

#endif // HAVE_NU1
