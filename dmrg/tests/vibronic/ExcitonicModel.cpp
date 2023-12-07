/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE MODEL_EXCITONIC_NU1

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

/** Test for the thiophene dimer based on the excitonic hamiltonian*/
BOOST_FIXTURE_TEST_CASE(Test_Thiophene_Dimer_Excitonic, VibronicFixture)
{
    #ifdef HAVE_U1
    auto lattice = Lattice(parametersVibronicThiopheneDimer);
    auto integrals = Vibrational::detail::parseIntegralExcitonic<double>(parametersVibronicThiopheneDimer, lattice);
    // Checks sizes
    BOOST_CHECK_EQUAL(integrals.first.size(), 3);
    BOOST_CHECK_EQUAL(integrals.second.size(), 3);
    BOOST_CHECK_EQUAL(integrals.first[0].size(), 2);
    // Creates the interface for a single excitation
    maquis::DMRGInterface<std::complex<double>> interface_single(parametersVibronicThiopheneDimer); 
    interface_single.optimize();
    double interface_single_optimizedEnergy = interface_single.energy().real();
    BOOST_CHECK_CLOSE(interface_single_optimizedEnergy, 0.0052660000, 1.0E-10);
    // Checks energy before and after time evolution. Value should be twice the ZPE of the system (given by the second line in the integral file)
    #ifdef DMRG_TD
    maquis::DMRGInterface<std::complex<double>> interface_single_TD(parametersVibronicThiopheneDimer);
    double interface_single_TD_initialEnergy = interface_single_TD.energy().real();
    interface_single_TD.evolve();
    double interface_single_TD_finalEnergy = interface_single_TD.energy().real();
    BOOST_CHECK_CLOSE(interface_single_TD_initialEnergy, 0.0052660000, 1.0E-10);
    BOOST_CHECK_CLOSE(interface_single_TD_finalEnergy, 0.0052660000, 1.0E-10);
    #endif //DMRG_TD
    //check for consistency when more than one exciton is created
    parametersVibronicThiopheneDimer.set("vibronic_num_excitons", 2); 
    parametersVibronicThiopheneDimer.set("init_basis_state", "1,0,1,0"); 
    // Creates the interface for a double excitation
    maquis::DMRGInterface<std::complex<double>> interface_double(parametersVibronicThiopheneDimer); 
    interface_double.optimize();
    double interface_double_optimizedEnergy = interface_double.energy().real();
    BOOST_CHECK_CLOSE(interface_double_optimizedEnergy, 0.0052660000, 1.0E-10);
    // Checks energy before and after time evolution. Value should be twice the ZPE of the system, since the ground and excited state PES is identical.
    #ifdef DMRG_TD
    maquis::DMRGInterface<std::complex<double>> interface_double_TD(parametersVibronicThiopheneDimer);
    double interface_double_TD_initialEnergy = interface_double_TD.energy().real();
    interface_double_TD.evolve();
    double interface_double_TD_finalEnergy = interface_double_TD.energy().real();
    BOOST_CHECK_CLOSE(interface_double_TD_initialEnergy, 0.0052660000, 1.0E-10);
    BOOST_CHECK_CLOSE(interface_double_TD_finalEnergy, 0.0052660000, 1.0E-10);
    #endif //DMRG_TD
    #endif //HAVE_U1
}

BOOST_FIXTURE_TEST_CASE(Test_Energy_Various_Nmax, VibronicFixture)
{
#ifdef HAVE_U1  
    auto lattice = Lattice(parametersTestNmax);
    auto integrals = Vibrational::detail::parseIntegralExcitonicExtended<double>(parametersTestNmax, lattice);
    //Checks sizes
    BOOST_CHECK_EQUAL(integrals.first.size(), 8);
    BOOST_CHECK_EQUAL(integrals.second.size(), 8);
    BOOST_CHECK_EQUAL(integrals.first[0].size(), 4);
    maquis::DMRGInterface<std::complex<double>> interface(parametersTestNmax); 
    interface.optimize();
    double interface_optimizedEnergy = interface.energy().real();
    BOOST_CHECK_CLOSE(interface_optimizedEnergy, 7.0, 1.0E-10);
    // Checks energy before and after time evolution. 
    #ifdef DMRG_TD
    maquis::DMRGInterface<std::complex<double>> interface_TD(parametersTestNmax);
    double interface_TD_initialEnergy = interface_TD.energy().real();
    interface_TD.evolve();
    double interface_TD_finalEnergy = interface_TD.energy().real();
    BOOST_CHECK_CLOSE(interface_TD_initialEnergy, 7.0, 1.0E-10);
    BOOST_CHECK_CLOSE(interface_TD_finalEnergy, 7.0, 1.0E-10);
    #endif //DMRG_TD
    #endif //HAVE_U1

}




/** Test for the integral parser of the extended excitonic hamiltonian */
BOOST_FIXTURE_TEST_CASE(Test_Integral_Parser_ExtendedExcitonic, VibronicFixture)
{
#ifdef HAVE_U1
    auto lattice = Lattice(parametersExcitonicExtendedAggregate); //add this to vibronic fixture
    auto integrals = Vibrational::detail::parseIntegralExcitonicExtended<double>(parametersExcitonicExtendedAggregate, lattice);
    //Checks sizes
    BOOST_CHECK_EQUAL(integrals.first.size(), 6); //checks number of rows
    BOOST_CHECK_EQUAL(integrals.second.size(), 6); //checks number of rows
    BOOST_CHECK_EQUAL(integrals.first[0].size(), 4); //check if two operators tags, an electronic state and site specific tag have been read in
#endif //HAVE_U1
}



