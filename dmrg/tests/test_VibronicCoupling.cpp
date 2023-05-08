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

/** @brief Watson-based calculation on a 4-mode vibronic Hamiltonian. */
BOOST_FIXTURE_TEST_CASE(Test_Vibronic_SingleSite, VibronicFixture)
{
#ifdef HAVE_U1
    using InterfaceType = maquis::DMRGInterface<double>;
    // Adds the final input parameters
    parametersVibronicPyrazineRedDim.set("init_type", "default");
    parametersVibronicPyrazineRedDim.set("optimization", "singlesite");
    parametersVibronicPyrazineRedDim.set("nsweeps", 10);
    parametersVibronicPyrazineRedDim.set("max_bond_dimension", 20);
    parametersVibronicPyrazineRedDim.set("ngrowsweeps", 2);
    parametersVibronicPyrazineRedDim.set("nmainsweeps", 2);
    parametersVibronicPyrazineRedDim.set("alpha_initial", 1.0E-8);
    parametersVibronicPyrazineRedDim.set("alpha_main", 1.0E-15);
    parametersVibronicPyrazineRedDim.set("alpha_final", 0.);
    // Creates the interface
    InterfaceType interface(parametersVibronicPyrazineRedDim);
    interface.optimize();
    // The reference value can be calculated based on the Harmonic approximation.
    auto refEnergy = -2036.1425;
    BOOST_CHECK_CLOSE(interface.energy(), refEnergy, 1.0E-5);
#endif
}

#endif // DMRG_VIBRONIC
