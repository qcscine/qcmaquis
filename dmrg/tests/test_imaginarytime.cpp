/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021- by Alberto Baiardi <abaiardi@ethz.ch>
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
#include "Fixtures/TimeEvolversFixture.h"
#include "Fixtures/PreBOTimeEvolversFixture.h"
#include "Fixtures/VibronicFixture.h"

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

#ifdef HAVE_U1DG

/**
 * @brief Tests that the energy is conserved along a relativistic TD-DMRG propagation.
 * The data are obtained for N2+ and the 3-21G basis set.
 */
BOOST_FIXTURE_TEST_CASE( TestImaginaryTimeRelativistic, TestTimeEvolverFixture )
{
    // Two-site evolutions
    parametersRelativistic.set("optimization", "twosite");
    parametersRelativistic.set("time_step", 10.);
    parametersRelativistic.set("nsweeps", 40);
    // TD
    maquis::DMRGInterface<std::complex<double>> interfaceTD(parametersRelativistic);
    interfaceTD.evolve();
    auto energyTD = std::real(interfaceTD.energy());
    // TI
    maquis::DMRGInterface<std::complex<double>> interfaceTI(parametersRelativistic);
    interfaceTI.optimize();
    auto energyTI = std::real(interfaceTI.energy());
    // The threshold is here a bit looser because the iTD-DMRG convergence is rather slow
    BOOST_CHECK_CLOSE(energyTD, energyTI, 1.0E-8);
}

#endif // HAVE_U1DG

#ifdef DMRG_PREBO

/**
 * @brief Tests that the energy is conserved along a "true" PreBO TD-DMRG propagation.
 */
BOOST_FIXTURE_TEST_CASE( TestImaginaryTimePreBO, PreBOTestTimeEvolverFixture )
{
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
}

#endif // DMRG_PREBO

#ifdef DMRG_VIBRONIC

/**
 * @brief Tests that the energy obtained with iTD-DMRG and DMRG is coherent for a vibronic Hamiltonian.
 */
BOOST_FIXTURE_TEST_CASE( TestImaginaryTimeVibronic, VibronicFixture )
{
#ifdef HAVE_U1
    parametersVibronicPyrazineRedDimFull.set("init_state", "basis_state_generic");
    parametersVibronicPyrazineRedDimFull.set("init_basis_state", "1,0,0,0,0,0");
    parametersVibronicPyrazineRedDimFull.set("nsweeps", 20);
    parametersVibronicPyrazineRedDimFull.set("max_bond_dimension", 20);
    parametersVibronicPyrazineRedDimFull.set("time_step", 0.1);
    parametersVibronicPyrazineRedDimFull.set("time_units", "as");
    parametersVibronicPyrazineRedDimFull.set("propagator_maxiter", 40);
    parametersVibronicPyrazineRedDimFull.set("TD_backpropagation", "no");
    parametersVibronicPyrazineRedDimFull.set("imaginary_time", "yes");
    maquis::DMRGInterface<double> realInterface(parametersVibronicPyrazineRedDimFull);
    maquis::DMRGInterface<std::complex<double>> complexInterface(parametersVibronicPyrazineRedDimFull);
    maquis::cout << "Running conventional DMRG optimization test for Vibronic model" << std::endl;
    realInterface.optimize();
    maquis::cout << "Running imaginary-time evolution for Vibronic model " << std::endl;
    complexInterface.evolve();
    // Test energy conservation
    auto TIEnergy = std::real(realInterface.energy());
    auto iTDEnergy = std::real(complexInterface.energy());
    BOOST_CHECK_CLOSE(TIEnergy, iTDEnergy, 1.0E-10);
#endif
}

#endif // DMRG_VIBRONIC