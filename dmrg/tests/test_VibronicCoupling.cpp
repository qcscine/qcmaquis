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

/** @brief Watson-based calculation on a 4-mode vibronic Hamiltonian. */
BOOST_FIXTURE_TEST_CASE(Test_Vibronic_SingleSite, VibronicFixture)
{
#ifdef HAVE_U1
    using InterfaceType = maquis::DMRGInterface<double, Hamiltonian::Vibronic>;
    // Adds the final input parameters
    parametersVibronicPyrazineRedDim.set("init_state", "default");
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
