/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2019 by Leon Freitag <lefreita@ethz.ch>
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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
