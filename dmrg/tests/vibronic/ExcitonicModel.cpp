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
    // Checks energy before and after time evolution. Value should be twice the ZPE of the system (given by the second line in the integral file)
    double interface_single_initialEnergy = interface_single.energy().real();
    interface_single.evolve();
    double interface_single_finalEnergy = interface_single.energy().real();
    BOOST_CHECK_CLOSE(interface_single_initialEnergy, 0.0052660000, 1.0E-10);
    BOOST_CHECK_CLOSE(interface_single_finalEnergy, 0.0052660000, 1.0E-10);
    //check for consistency when more than one exciton is created
    parametersVibronicThiopheneDimer.set("vibronic_num_excitons", 2); 
    parametersVibronicThiopheneDimer.set("init_basis_state", "1,0,1,0"); 
    // Creates the interface for a double excitation
    maquis::DMRGInterface<std::complex<double>> interface_double(parametersVibronicThiopheneDimer); 
    // Checks energy before and after time evolution. Value should be twice the ZPE of the system, since the ground and excited state PES is identical.
    double interface_double_initialEnergy = interface_double.energy().real();
    interface_double.evolve();
    double interface_double_finalEnergy = interface_double.energy().real();
    BOOST_CHECK_CLOSE(interface_double_initialEnergy, 0.0052660000, 1.0E-10);
    BOOST_CHECK_CLOSE(interface_double_finalEnergy, 0.0052660000, 1.0E-10);
    #endif //HAVE_U1
}

