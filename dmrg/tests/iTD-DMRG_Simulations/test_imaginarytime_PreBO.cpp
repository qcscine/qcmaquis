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
#include "Fixtures/PreBOTimeEvolversFixture.h"

/**
 * @brief Tests that the energy is conserved along a "true" PreBO TD-DMRG propagation.
 */
BOOST_FIXTURE_TEST_CASE( TestImaginaryTimePreBO, PreBOTestTimeEvolverFixture )
{
#ifdef DMRG_PREBO
    // Generic settings
    parametersPreBOComplex.set("optimization", "twosite");
    parametersPreBOReal.set("optimization", "twosite");
    maquis::DMRGInterface<double> realInterface(parametersPreBOReal);
    maquis::DMRGInterface<std::complex<double>> complexInterface(parametersPreBOComplex);
    maquis::cout << "Running conventional DMRG optimization test for PreBO model" << std::endl;
    realInterface.optimize();
    maquis::cout << "Running imaginary-time evolution for PreBO model " << std::endl;
    complexInterface.evolve();
    // Test energy conservation
    auto TIEnergy = std::real(realInterface.energy());
    auto iTDEnergy = std::real(complexInterface.energy());
    BOOST_CHECK_CLOSE(TIEnergy, iTDEnergy, 1.0E-10);
#endif // DMRG_PREBO
}