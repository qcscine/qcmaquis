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
#include "Fixtures/NModeFixture.h"
#include "maquis_dmrg.h"

/**
 * @brief N-mode calculation on two-mode PES of FAD.
 * The basis set has been generated with a DVR primitive basis.
 */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_2ModeSystem, NModeFixture)
{
#ifdef HAVE_NU1
    // Adds the final input parameters
    parametersFADTwoBody.set("init_state", "const");
    parametersFADTwoBody.set("nsweeps", 20);
    parametersFADTwoBody.set("max_bond_dimension",100);
    parametersFADTwoBody.set("MODEL", "nmode");
    // Creates the interface
    maquis::DMRGInterface<double, Hamiltonian::VibrationalNMode> interface(parametersFADTwoBody);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), -1499.5871477508479, 1.0E-5);
#endif // HAS_NU1
}

#ifdef HAVE_NU1

/** @brief Same as above, but includes measurements */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_2ModeSystem_MeasOneParticle, NModeFixture)
{
    // Adds the final input parameters
    parametersFADTwoBody.set("init_state", "default");
    parametersFADTwoBody.set("seed", 16071991);
    parametersFADTwoBody.set("nsweeps", 100);
    parametersFADTwoBody.set("max_bond_dimension", 20);
    parametersFADTwoBody.set("MODEL", "nmode");
    parametersFADTwoBody.set("MEASURE[One Modal RDM]", "1");
    // Creates the interface
    maquis::DMRGInterface<double, Hamiltonian::VibrationalNMode> interface(parametersFADTwoBody);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), -1499.5871477508479, 1.0E-5);
    // Measurements
    interface.measure();
    // std::cout << interface.getMeasurement("onemodalRDM_0_00").second.size() << std::endl;
    for (int iSite = 0; iSite < 22; iSite++) {
        BOOST_CHECK_EQUAL(interface.getMeasurement("onemodalRDM_00").first[iSite][0],
                          interface.getMeasurement("onemodalRDM_11").first[iSite][0]);
        auto meas1 = interface.getMeasurement("onemodalRDM_00").second[iSite];
        auto meas2 = interface.getMeasurement("onemodalRDM_11").second[iSite];
        BOOST_CHECK_CLOSE(meas1+meas2, 1., 1.0E-16);
    }
}

/** 
 * @brief Measure the two-modal RDM and checks sanity of the modal entropies
 * The sanity check is done by verifying that the subadditivity property is verified
 * (i.e., that S_{pq} < S_p + S_q, see Chem. Phys. 323, 519 (2006))
 */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_2ModeSystem_Subadditivity, NModeFixture)
{
    // Adds the final input parameters
    parametersFADTwoBody.set("init_state", "default");
    parametersFADTwoBody.set("seed", 16071991);
    parametersFADTwoBody.set("nsweeps", 100);
    parametersFADTwoBody.set("max_bond_dimension", 20);
    parametersFADTwoBody.set("MODEL", "nmode");
    parametersFADTwoBody.set("MEASURE[One Modal RDM]", "1");
    parametersFADTwoBody.set("MEASURE[Two Modal RDM]", "1");
    // Creates the interface
    maquis::DMRGInterface<double, Hamiltonian::VibrationalNMode> interface(parametersFADTwoBody);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), -1499.5871477508479, 1.0E-5);
    // Measurements
    interface.measure();
    // TODO COMPLETE!!
}

#endif // HAS_NU1

#endif // DMRG_VIBRATIONAL
