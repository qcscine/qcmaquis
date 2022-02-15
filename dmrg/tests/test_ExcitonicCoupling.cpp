/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifdef DMRG_VIBRONIC

#include <boost/test/included/unit_test.hpp>
#include "Fixtures/VibronicFixture.h"
#include "maquis_dmrg.h"

/** @brief Watson-based calculation on an excitonic Hamiltonian. */
BOOST_FIXTURE_TEST_CASE(Test_Vibronic_Excitonic_SingleSite, VibronicFixture)
{
#ifdef HAVE_U1
    using InterfaceType = maquis::DMRGInterface<double, Hamiltonian::Excitonic>;
    // Adds the final input parameters
    parametersExcitonicAggregate.set("init_state", "default");
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
    using InterfaceType = maquis::DMRGInterface<double, Hamiltonian::Excitonic>;
    // Adds the final input parameters
    parametersExcitonicAggregate.set("init_state", "default");
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
