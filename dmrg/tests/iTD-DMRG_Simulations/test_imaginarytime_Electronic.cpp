/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
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
#include "utils/io.hpp"
#include <iostream>
#include "maquis_dmrg.h"
#include "Fixtures/BenzeneFixture.h"

/**
 * @brief Tests that iTD-DMRG and TI-DMRG give the same energy.
 * 
 * The data are obtained for CAS(6, 6) and based on the cc-pVDZ basis set.
 * The data are stored in the BenzeneFixture class.
 */
BOOST_FIXTURE_TEST_CASE( TestImaginaryTime, BenzeneFixture )
{
    std::vector<std::string> symmetries;
    #ifdef HAVE_SU2U1PG
    symmetries.push_back("su2u1pg");
    #endif
    #ifdef HAVE_SU2U1
    symmetries.push_back("su2u1");
    #endif
    #ifdef HAVE_TwoU1PG
    symmetries.push_back("2u1pg");
    #endif
    #ifdef HAVE_TwoU1
    symmetries.push_back("2u1");
    #endif

    for (auto&& s: symmetries) {
        maquis::cout << "Running imaginary-time evolution test for symmetry " << s << std::endl;
        parametersBenzeneImaginaryTime.set("symmetry", s);
        parametersBenzene.set("symmetry", s);
        maquis::DMRGInterface<double> realInterface(parametersBenzene);
        maquis::DMRGInterface<std::complex<double>> complexInterface(parametersBenzeneImaginaryTime);
        realInterface.optimize();
        complexInterface.evolve();
        // Test energy conservation
        auto TIEnergy = std::real(realInterface.energy());
        auto iTDEnergy = std::real(complexInterface.energy());
        BOOST_CHECK_CLOSE(TIEnergy, iTDEnergy, 1.0E-10);
    }
}