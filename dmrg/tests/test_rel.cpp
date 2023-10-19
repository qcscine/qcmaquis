/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include "utils/io.hpp" // has to be first include because of impi
#include <iostream>
#include "Fixtures/RelativisticFixture.h"

/** @brief N2 energy with 4c-DMRG */
BOOST_FIXTURE_TEST_CASE( Test_Relativistic, RelativisticFixture )
{
    // only enable the test if we compile the support for the U1DG symmetry
#ifdef HAVE_U1DG
    maquis::DMRGInterface<std::complex<double>> interface(parametersComplex);
    interface.optimize();
    // Checks energy
    BOOST_CHECK_CLOSE(std::real(interface.energy()), -1.0780470133e+02 , 1e-7);
    const typename maquis::DMRGInterface<std::complex<double>>::meas_with_results_type& meas = interface.onerdm();
    // Calculate the trace of the 1-RDM and check if it adds up to the number of electrons
    std::complex<double> value = 0.0;
    for (int i = 0; i < meas.first.size(); i++)
        if (meas.first[i][0] == meas.first[i][1]) // sum up diagonal elements for the trace
            value += meas.second[i];
    BOOST_CHECK_CLOSE(value.real(), 3.0 , 1e-7);
#endif
}
