/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021-2022 by Alberto Baiardi <abaiardi@ethz.ch>
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
 * @brief Tests that the energy is conserved along a TD-DMRG propagation.
 * 
 * The data are obtained for CAS(6, 6) and based on the cc-pVDZ basis set.
 * The data are stored in the BenzeneFixture class.
 */
BOOST_FIXTURE_TEST_CASE( TestRealTime, BenzeneFixture )
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
        parametersBenzeneRealTime.set("symmetry", s);
        // Single-site evolution
        maquis::cout << "Running SS real-time evolution test for symmetry " << s << std::endl;
        parametersBenzeneRealTime.set("optimization", "singlesite");
        maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersBenzeneRealTime);
        auto initialEnergy = std::real(interfaceSS.energy());
        interfaceSS.evolve();
        auto finalEnergy = std::real(interfaceSS.energy());
        BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-8);
        // Two-site evolution
        maquis::cout << "Running TS real-time evolution test for symmetry " << s << std::endl;
        parametersBenzeneRealTime.set("optimization", "twosite");
        maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersBenzeneRealTime);
        initialEnergy = std::real(interfaceTS.energy());
        interfaceTS.evolve();
        finalEnergy = std::real(interfaceTS.energy());
        BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
    }
}

#ifdef HAVE_U1DG

/**
 * @brief Tests that the energy is conserved along a relativistic TD-DMRG propagation.
 * The data are obtained for N2+ and the 3-21G basis set.
 */
BOOST_FIXTURE_TEST_CASE( TestRealTimeRelativistic, TestTimeEvolverFixture )
{
    // Generic settings
    parametersRelativistic.set("imaginary_time", "no");
    parametersRelativistic.set("TD_backpropagation", "yes");
    parametersRelativistic.set("time_units", "as");
    parametersRelativistic.set("time_step", 1.0E+00);
    // Single-site evolution
    parametersRelativistic.set("optimization", "singlesite");
    maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersRelativistic);
    auto initialEnergy = std::real(interfaceSS.energy());
    interfaceSS.evolve();
    auto finalEnergy = std::real(interfaceSS.energy());
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
    // Two-site evolutions
    parametersRelativistic.set("optimization", "twosite");
    maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersRelativistic);
    initialEnergy = std::real(interfaceTS.energy());
    interfaceTS.evolve();
    finalEnergy = std::real(interfaceTS.energy());
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
}

#endif

#ifdef DMRG_PREBO

/**
 * @brief Tests that the energy is conserved along a PreBO TD-DMRG propagation.
 * 
 * Note that this example uses the BO hamiltonian for H2 in a minimal basis set - i.e.,
 * only electrons are included.
 */
BOOST_FIXTURE_TEST_CASE( TestRealTimePreBOWithBOHamiltonian, PreBOTestTimeEvolverFixture )
{
    // Generic settings
    parametersBOComplex.set("imaginary_time", "no");
    parametersBOComplex.set("TD_backpropagation", "yes");
    parametersBOComplex.set("time_units", "as");
    parametersBOComplex.set("time_step", 1.0E+00);
    // Single-site evolution
    maquis::cout << "Running SS real-time evolution test for PreBO model " << std::endl;
    parametersBOComplex.set("optimization", "singlesite");
    maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersBOComplex);
    auto initialEnergy = std::real(interfaceSS.energy());
    interfaceSS.evolve();
    auto finalEnergy = std::real(interfaceSS.energy());
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
    // Two-site evolutions
    maquis::cout << "Running TS real-time evolution test for PreBO model" << std::endl;
    parametersBOComplex.set("optimization", "twosite");
    maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersBOComplex);
    initialEnergy = std::real(interfaceTS.energy());
    interfaceTS.evolve();
    finalEnergy = std::real(interfaceTS.energy());
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
}

/**
 * @brief Tests that the energy is conserved along a "true" PreBO TD-DMRG propagation.
 */
BOOST_FIXTURE_TEST_CASE( TestRealTimePreBO, PreBOTestTimeEvolverFixture )
{
    // Generic settings
    parametersPreBOComplex.set("imaginary_time", "no");
    parametersPreBOComplex.set("TD_backpropagation", "yes");
    parametersPreBOComplex.set("time_units", "as");
    parametersPreBOComplex.set("time_step", 1.0E+00);
    // Single-site evolution
    maquis::cout << "Running SS real-time evolution test for PreBO model " << std::endl;
    parametersPreBOComplex.set("optimization", "singlesite");
    maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersPreBOComplex);
    auto initialEnergy = std::real(interfaceSS.energy());
    interfaceSS.evolve();
    auto finalEnergy = std::real(interfaceSS.energy());
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
    // Two-site evolutions
    maquis::cout << "Running TS real-time evolution test for PreBO model" << std::endl;
    parametersPreBOComplex.set("optimization", "twosite");
    maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersPreBOComplex);
    initialEnergy = std::real(interfaceTS.energy());
    interfaceTS.evolve();
    finalEnergy = std::real(interfaceTS.energy());
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
}

#endif // DMRG_PREBO

#ifdef DMRG_VIBRONIC

/** @brief Tests that the energy is conserved along an Excitonic TD-DMRG propagation. */
BOOST_FIXTURE_TEST_CASE( TestRealTimeExcitonic, VibronicFixture )
{
#ifdef HAVE_U1
    // Generic settings
    parametersExcitonicAggregateTwoSites.set("imaginary_time", "no");
    parametersExcitonicAggregateTwoSites.set("TD_backpropagation", "yes");
    parametersExcitonicAggregateTwoSites.set("propagator_maxiter", 30);
    parametersExcitonicAggregateTwoSites.set("propagator_accuracy", 1.0E-10);
    parametersExcitonicAggregateTwoSites.set("time_units", "as");
    parametersExcitonicAggregateTwoSites.set("time_step", 0.01);
    parametersExcitonicAggregateTwoSites.set("nsweeps", 2);
    parametersExcitonicAggregateTwoSites.set("max_bond_dimension", 50);
    parametersExcitonicAggregateTwoSites.set("init_state", "basis_state_generic");
    parametersExcitonicAggregateTwoSites.set("init_basis_state", "1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
    // Single-site evolution
    maquis::cout << "Running SS real-time evolution test for Excitonic model " << std::endl;
    parametersExcitonicAggregateTwoSites.set("optimization", "singlesite");
    maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersExcitonicAggregateTwoSites);
    auto initialEnergy = std::real(interfaceSS.energy());
    interfaceSS.evolve();
    auto finalEnergy = std::real(interfaceSS.energy());
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
    // Two-site evolutions
    maquis::cout << "Running TS real-time evolution test for Excitonic model" << std::endl;
    parametersExcitonicAggregateTwoSites.set("optimization", "twosite");
    maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersExcitonicAggregateTwoSites);
    initialEnergy = std::real(interfaceTS.energy());
    interfaceTS.evolve();
    finalEnergy = std::real(interfaceTS.energy());
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
#endif // HAVE_U1
}

#endif // DMRG_VIBRONIC
