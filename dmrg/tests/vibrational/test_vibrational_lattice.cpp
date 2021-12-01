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
    auto size = lattice.maximum_vertex_type();
    BOOST_CHECK_EQUAL(size, 1);
}

/** @brief Checks that the property getter works properly */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_PropertyGetter_Watson_Ethylene, WatsonFixture)
{
    auto lattice = WatsonLattice(parametersEthyleneWatson);
    for (int iSite = 0; iSite < lattice.size(); iSite++) {
        auto siteType = lattice.get_prop<int>("type", iSite);
        BOOST_CHECK_EQUAL(siteType, 0);
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
    auto typeOfSites = lattice.maximum_vertex_type();
    BOOST_CHECK_EQUAL(typeOfSites, 1);
}

/** @brief Checks the size of the lattice for the 4-mode input */
BOOST_FIXTURE_TEST_CASE(Test_Site_Types_4ModeSystem, NModeFixture)
{
    auto lattice = NModeLattice(parametersFourMode);
    auto typeOfSites = lattice.maximum_vertex_type();
    BOOST_CHECK_EQUAL(typeOfSites, 3);
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

#endif // HAVE_NU1