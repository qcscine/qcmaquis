/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#ifdef DMRG_VIBRATIONAL

#include <boost/test/included/unit_test.hpp>
#include "Fixtures/WatsonFixture.h"
#include "maquis_dmrg.h"

/**
 * @brief Watson-based calculation on the harmonic ethylene PES
 */
BOOST_FIXTURE_TEST_CASE(Test_vDMRG_Calculation_Ethylene_Harmonic, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
    using InterfaceType = maquis::DMRGInterface<double>;
    // Adds the final input parameters
    parametersEthyleneWatsonHarmonic.set("init_type", "const");
    parametersEthyleneWatsonHarmonic.set("nsweeps", 20);
    parametersEthyleneWatsonHarmonic.set("max_bond_dimension", 20);
    parametersEthyleneWatsonHarmonic.set("MODEL", "watson");
    // Creates the interface
    InterfaceType interface(parametersEthyleneWatsonHarmonic);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), referenceHarmonicEnergy, 1.0E-5);
    // Checks that the overlap of the final wave function with the hf determinant is = 1.
    auto targetOverlap = interface.getCICoefficient("0,0,0,0,0,0,0,0,0,0,0,0");
    BOOST_CHECK_CLOSE(std::abs(targetOverlap), 1.0, 1.0E-12);
    // Checks that overlap calculations cannot be performed for invalid reference determinants
    BOOST_CHECK_THROW(interface.getCICoefficient("0,0,0"), std::runtime_error);
#endif
}

#ifdef HAVE_TrivialGroup

/**
 * @brief Watson-based calculation on the harmonic ethylene PES, but with input given only
 * on two rows.
 */
BOOST_FIXTURE_TEST_CASE(Test_vDMRG_Calculation_Ethylene_Harmonic_TwoRows, WatsonFixture)
{
    using InterfaceType = maquis::DMRGInterface<double>;
    // Adds the final input parameters
    parametersEthyleneWatsonHarmonicTwoRows.set("init_type", "const");
    parametersEthyleneWatsonHarmonicTwoRows.set("nsweeps", 20);
    parametersEthyleneWatsonHarmonicTwoRows.set("max_bond_dimension", 20);
    parametersEthyleneWatsonHarmonicTwoRows.set("MODEL", "watson");
    // Creates the interface
    InterfaceType interface(parametersEthyleneWatsonHarmonicTwoRows);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), referenceHarmonicEnergy, 1.0E-5);
    // Checks that the overlap of the final wave function with the hf determinant is = 1.
    auto targetOverlap = interface.getCICoefficient("0,0,0,0,0,0,0,0,0,0,0,0");
    BOOST_CHECK_CLOSE(std::abs(targetOverlap), 1.0, 1.0E-12);
}

/**
 * @brief Watson-based calculation on the ethylene PES
 * The PES has been taken from the vHBCI reference work, which is
 * J. Chem. Phys., 154, 074104 (2021).
 * Note that the calculation uses the single-site optimizer.
 */
BOOST_FIXTURE_TEST_CASE(Test_vDMRG_Calculation_Ethylene_Sextic_SingleSite, WatsonFixture)
{
    using InterfaceType = maquis::DMRGInterface<double>;
    // Adds the final input parameters
    parametersEthyleneWatson.set("init_type", "basis_state_generic");
    parametersEthyleneWatson.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
    parametersEthyleneWatson.set("nsweeps", 20);
    parametersEthyleneWatson.set("max_bond_dimension", 50);
    parametersEthyleneWatson.set("MODEL", "watson");
    parametersEthyleneWatson.set("optimization", "singlesite");
    parametersEthyleneWatson.set("ngrowsweeps", 2);
    parametersEthyleneWatson.set("nmainsweeps", 2);
    parametersEthyleneWatson.set("alpha_initial", 1.0E-8);
    parametersEthyleneWatson.set("alpha_main", 1.0E-15);
    parametersEthyleneWatson.set("alpha_final", 0.);
    parametersEthyleneWatson.set("Nmax", 7);
    parametersEthyleneWatson.set("truncation_initial", 1.0E-10);
    parametersEthyleneWatson.set("truncation_final", 1.0E-8);
    // Creates the interface
    InterfaceType interface(parametersEthyleneWatson);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), 11011.61, 1.0E-2);
}

/**
 * @brief Watson-based calculation on the ethylene PES
 * The PES has been taken from the vHBCI reference work, which is
 * J. Chem. Phys., 154, 074104 (2021).
 * Note that the calculation uses the single-site optimizer.
 */
BOOST_FIXTURE_TEST_CASE(Test_vDMRG_Calculation_Ethylene_Sextic_TwoSite, WatsonFixture)
{
    using InterfaceType = maquis::DMRGInterface<double>;
    // Adds the final input parameters
    parametersEthyleneWatson.set("init_type", "basis_state_generic");
    parametersEthyleneWatson.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
    parametersEthyleneWatson.set("nsweeps", 10);
    parametersEthyleneWatson.set("max_bond_dimension", 20);
    parametersEthyleneWatson.set("MODEL", "watson");
    parametersEthyleneWatson.set("optimization", "twosite");
    parametersEthyleneWatson.set("ngrowsweeps", 2);
    parametersEthyleneWatson.set("nmainsweeps", 2);
    parametersEthyleneWatson.set("truncation_initial", 1.0E-16);
    parametersEthyleneWatson.set("truncation_final", 1.0E-10);
    // Creates the interface
    InterfaceType interface(parametersEthyleneWatson);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), 11011.61, 1.0E-2);
}

#endif // HAVE_TrivialGroup

#endif // DMRG_VIBRATIONAL
