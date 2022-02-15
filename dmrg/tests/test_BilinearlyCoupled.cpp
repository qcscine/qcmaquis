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

#ifdef DMRG_VIBRATIONAL

#include <boost/test/included/unit_test.hpp>
#include "Fixtures/WatsonFixture.h"
#include "maquis_dmrg.h"

/**
 * @brief Watson-based calculation on a bilinearly-coupled Harmonic Hamiltonian
 * Note that the reference energy is taken from the work
 * "Calculating vibrational spectra with sum of product basis functions without
 * storing full-dimensional vectors or matrices"
 * by the group of Tucker Carrington
 */
BOOST_FIXTURE_TEST_CASE(Test_vDMRG_Calculation_BilinearlyCoupled, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
    using InterfaceType = maquis::DMRGInterface<double, Hamiltonian::VibrationalCanonical>;
    // Adds the final input parameters
    parametersBilinearly.set("init_state", "basis_state_generic");
    parametersBilinearly.set("init_basis_state", "0,0,0,0,0,0");
    parametersBilinearly.set("optimization", "singlesite");
    parametersBilinearly.set("alpha_initial", 1.0E-8);
    parametersBilinearly.set("alpha_initial", 1.0E-15);
    parametersBilinearly.set("alpha_initial", 0.);
    parametersBilinearly.set("nsweeps", 20);
    parametersBilinearly.set("ngrowsweeps", 2);
    parametersBilinearly.set("nmainsweeps", 2);
    parametersBilinearly.set("max_bond_dimension", 20);
    parametersBilinearly.set("MODEL", "watson");
    // Creates the interface
    InterfaceType interface(parametersBilinearly);
    interface.optimize();
    // The reference value can be calculated based on the Harmonic approximation
    BOOST_CHECK_CLOSE(interface.energy(), 3.8164041, 1.0E-5);
#endif
}

#endif // DMRG_VIBRATIONAL
