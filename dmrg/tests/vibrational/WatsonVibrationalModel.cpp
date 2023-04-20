/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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

/** Checks consistency for the physical dimensions for the ethylene Watson Hamiltonian with NMax vector initialization */
BOOST_FIXTURE_TEST_CASE(Test_Model_PhysDim_Ethylene_NMaxVec, WatsonFixture)
{
    parametersEthyleneWatsonHarmonic.set("Nmax", "8");
    auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto nModeModel = WatsonHamiltonian<matrix>(lattice, parametersEthyleneWatsonHarmonic, false);
    int siteType = 5;
    const auto& physicalDimensions5 = nModeModel.phys_dim(siteType);
    BOOST_CHECK_EQUAL(physicalDimensions5.sum_of_sizes(), 8);
}


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
    auto integrals = Vibrational::detail::WatsonIntegralParser<double>(parametersEthyleneWatsonHarmonic, lattice, WatsonCoordinateType::CartesianNormalModes,
                                                                       6, 6, 6);
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
    int nMax = parametersEthyleneWatsonHarmonic["Nmax"];
    BOOST_CHECK_EQUAL(physicalDimensions0.sum_of_sizes(), nMax);
}

#endif // HAVE_TrivialGroup
