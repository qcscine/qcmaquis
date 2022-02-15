/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Alberto Baiardi <robin.feldmann@phys.chem.ethz.ch>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#define BOOST_TEST_MODULE MODEL_VIBRATIONAL_NONE

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/vibrational/none/model.hpp"
#include "Fixtures/WatsonFixture.h"
#include "maquis_dmrg.h"
#include "dmrg/sim/matrix_types.h"

// Note that here below at least one test is defined also if TrivialGroup is not available
// (otherwise boost will complain)

/** Checks consistency for the physical dimensions for the ethylene Watson Hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Model_PhysDim_Ethylene, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
    auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto nModeModel = WatsonHamiltonian<matrix>(lattice, parametersEthyleneWatsonHarmonic, false);
    int siteType = 0;
    const auto& physicalDimensions0 = nModeModel.phys_dim(siteType);
    BOOST_CHECK_EQUAL(physicalDimensions0.sum_of_sizes(), 8);
#endif // HAVE_TrivialGroup
}

#ifdef HAVE_TrivialGroup

/** Simple check on tags */
BOOST_FIXTURE_TEST_CASE(Test_Model_Tag_SimpleCheck_Ethylene, WatsonFixture)
{
    auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto watsonModel = WatsonHamiltonian<matrix>(lattice, parametersEthyleneWatsonHarmonic, false);
    auto identityTag = watsonModel.identity_matrix_tag(0);
    auto fillingTag = watsonModel.filling_matrix_tag(0);
    // The nMode Hamiltonian is bosonic, so the tag should be the same
    BOOST_CHECK(identityTag == fillingTag);
}

/** Check on symbolic operator getter */
BOOST_FIXTURE_TEST_CASE(Test_Model_Symbolic_Operator_Ethylene, WatsonFixture)
{
    auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto watsonModel = WatsonHamiltonian<matrix>(lattice, parametersEthyleneWatsonHarmonic, false);
    BOOST_CHECK(watsonModel.filling_matrix_tag(0) == watsonModel.get_operator_tag("fill", 0));
}

/** Check consistency in the dimension of the integral container */
BOOST_FIXTURE_TEST_CASE(Test_Model_Watson_Ethylene_IntegralContainer, WatsonFixture) 
{
    auto lattice = lattice_factory(parametersEthyleneWatsonHarmonic);
    auto integrals = Vibrational::detail::WatsonIntegralParser<double>(parametersEthyleneWatsonHarmonic, lattice);
    BOOST_CHECK_EQUAL(integrals.size(), 24);
    for (const auto iElements: integrals) {
        for (int iSite = 2; iSite < 6; iSite++)
            BOOST_CHECK_EQUAL(iElements.first[iSite], 0);
    }
}

/** Check consistency in the dimension of the integral container */
BOOST_FIXTURE_TEST_CASE(Test_Model_Watson_Ethylene_TermsSize, WatsonFixture) 
{
    auto lattice = lattice_factory(parametersEthyleneWatsonHarmonic);
    auto watsonModel = WatsonHamiltonian<matrix>(lattice, parametersEthyleneWatsonHarmonic, false);
    watsonModel.create_terms();
    BOOST_CHECK_EQUAL(watsonModel.hamiltonian_terms().size(), 24);
}

/** Check consistency in the dimension of the integral container */
BOOST_FIXTURE_TEST_CASE(Test_Model_Watson_Ethylene_PhysDim, WatsonFixture) 
{
    auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto watsonModel = WatsonHamiltonian<matrix>(lattice, parametersEthyleneWatsonHarmonic, false);
    const auto& physicalDimensions0 = watsonModel.phys_dim(0);
    BOOST_CHECK_EQUAL(physicalDimensions0.sum_of_sizes(), parametersEthyleneWatsonHarmonic["Nmax"]);
}

#endif // HAVE_TrivialGroup
