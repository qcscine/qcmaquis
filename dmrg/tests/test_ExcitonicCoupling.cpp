/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#ifdef DMRG_VIBRONIC

#include <boost/test/included/unit_test.hpp>
#include "Fixtures/VibronicFixture.h"
#include "maquis_dmrg.h"

/** @brief Watson-based calculation on an excitonic Hamiltonian. */
BOOST_FIXTURE_TEST_CASE(Test_Vibronic_Excitonic_SingleSite, VibronicFixture)
{
#ifdef HAVE_U1
    using InterfaceType = maquis::DMRGInterface<double>;
    // Adds the final input parameters
    parametersExcitonicAggregate.set("init_type", "default");
    //parametersExcitonicAggregate.set("init_basis_state", "1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
    parametersExcitonicAggregate.set("optimization", "singlesite");
    parametersExcitonicAggregate.set("nsweeps", 10);
    parametersExcitonicAggregate.set("max_bond_dimension", 20);
    parametersExcitonicAggregate.set("integral_file", "integral_file_Excitonic_Harmonic");
    parametersExcitonicAggregate.set("ngrowsweeps", 2);
    parametersExcitonicAggregate.set("nmainsweeps", 2);
    parametersExcitonicAggregate.set("alpha_initial", 1.0E-8);
    parametersExcitonicAggregate.set("alpha_main", 1.0E-15);
    parametersExcitonicAggregate.set("alpha_final", 0.);
    // Note that here we are imposing that J=0
    parametersExcitonicAggregate.set("J_coupling", 0.);
    // Creates the interface
    InterfaceType interface(parametersExcitonicAggregate);
    interface.optimize();
    // The reference value can be calculated based on the Harmonic approximation.
    auto refEnergy = 4812.5;
    BOOST_CHECK_CLOSE(interface.energy(), refEnergy, 1.0E-5);
#endif
}

#ifdef HAVE_U1

/** @brief Watson-based calculation on an excitonic Hamiltonian. */
BOOST_FIXTURE_TEST_CASE(Test_Vibronic_Excitonic_TwoSite, VibronicFixture)
{
    using InterfaceType = maquis::DMRGInterface<double>;
    // Adds the final input parameters
    parametersExcitonicAggregate.set("init_type", "default");
    //parametersExcitonicAggregate.set("init_basis_state", "1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
    parametersExcitonicAggregate.set("optimization", "twosite");
    parametersExcitonicAggregate.set("nsweeps", 10);
    parametersExcitonicAggregate.set("max_bond_dimension", 20);
    parametersExcitonicAggregate.set("integral_file", "integral_file_Excitonic_Harmonic");
    parametersExcitonicAggregate.set("J_coupling", 0.);
    parametersExcitonicAggregate.set("J_excitation", 1000.);
    // Creates the interface
    InterfaceType interface(parametersExcitonicAggregate);
    interface.optimize();
    // The reference value can be calculated based on the Harmonic approximation.
    // Note that here we must add the vertical energy.
    auto refEnergy = 4812.5 + 1000;
    BOOST_CHECK_CLOSE(interface.energy(), refEnergy, 1.0E-5);
}

#endif // HAVE_U1

#endif // DMRG_VIBRONIC
