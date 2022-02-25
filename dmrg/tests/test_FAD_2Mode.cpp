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
#include "dmrg/sim/matrix_types.h"
// Needed for the diagonalization
#include "dmrg/utils/utils.hpp"
#include "utils/timings.h"
#include "utils/traits.hpp"
#include "utils/bindings.hpp"

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
    parametersFADTwoBody.set("nsweeps", 10);
    parametersFADTwoBody.set("max_bond_dimension", 20);
    parametersFADTwoBody.set("MODEL", "nmode");
    parametersFADTwoBody.set("MEASURE[One Modal RDM]", "1");
    parametersFADTwoBody.set("MEASURE[Two Modal RDM]", "1");
    // Creates the interface
    maquis::DMRGInterface<double, Hamiltonian::VibrationalNMode> interface(parametersFADTwoBody);
    interface.optimize();
    interface.measure();
    // == TWO-MODAL ENTROPY ==
    tmatrix<double> twoMatrixEntanglement(22, 22, 0.);
    BOOST_CHECK_EQUAL(interface.getMeasurement("twomodeRDM_00").first.size(),
                      interface.getMeasurement("twomodeRDM_11").first.size());
    BOOST_CHECK_EQUAL(interface.getMeasurement("twomodeRDM_11").first.size(),
                      interface.getMeasurement("twomodeRDM_12").first.size());
    BOOST_CHECK_EQUAL(interface.getMeasurement("twomodeRDM_21").first.size(),
                      interface.getMeasurement("twomodeRDM_12").first.size());
    BOOST_CHECK_EQUAL(interface.getMeasurement("twomodeRDM_21").first.size(),
                      interface.getMeasurement("twomodeRDM_22").first.size());
    BOOST_CHECK_EQUAL(interface.getMeasurement("twomodeRDM_33").first.size(),
                      interface.getMeasurement("twomodeRDM_22").first.size());
    auto overallSize = interface.getMeasurement("twomodeRDM_00").first.size();
    for (int iElement = 0; iElement < overallSize; iElement++) {
        // Gets the indices
        int iRow = interface.getMeasurement("twomodeRDM_00").first[iElement][0];
        int iCol = interface.getMeasurement("twomodeRDM_00").first[iElement][1];
        // Extracts the matrix elements
        auto meas00 = interface.getMeasurement("twomodeRDM_00").second[iElement];
        auto meas11 = interface.getMeasurement("twomodeRDM_11").second[iElement];
        auto meas12 = interface.getMeasurement("twomodeRDM_12").second[iElement];
        auto meas21 = interface.getMeasurement("twomodeRDM_21").second[iElement];
        auto meas22 = interface.getMeasurement("twomodeRDM_22").second[iElement];
        auto meas33 = interface.getMeasurement("twomodeRDM_33").second[iElement];
        // Constructs the two-modal entanglement matrix
        tmatrix<double> entanglementMatrix(4, 4, 0.), evecs(4, 4, 0.);
        std::vector<double> evals(4, 0.);
        entanglementMatrix(0, 0) = meas00;
        entanglementMatrix(1, 1) = meas11;
        entanglementMatrix(1, 2) = meas12;
        entanglementMatrix(2, 1) = meas21;
        entanglementMatrix(2, 2) = meas22;
        entanglementMatrix(3, 3) = meas33;
        heev(entanglementMatrix, evecs, evals);
        // Calculates the two-orbital entropy
        double entropy = 0.;
        for (int i = 0; i < 4; i++)
            if (std::abs(evals[i]) > 1.0E-16)
                entropy -= evals[i]*std::log(evals[i]);
        twoMatrixEntanglement(iRow, iCol) = entropy;
    }
    // == ONE-MODAL ENTROPY ==
    std::vector<double> vectorEntanglement(22, 0.);
    for (int iElement = 0; iElement < 22; iElement++) {
        auto meas0 = interface.getMeasurement("onemodalRDM_00").second[iElement];
        auto meas1 = interface.getMeasurement("onemodalRDM_11").second[iElement];
        if (std::abs(meas0) > 1.0E-16)
            vectorEntanglement[iElement] -= meas0*std::log(meas0);
        if (std::abs(meas1) > 1.0E-16)
            vectorEntanglement[iElement] -= meas1*std::log(meas1);
    }
    // Checks subadditivity
    for (int iRow = 0; iRow < 22; iRow++) {
        for (int iCol = 0; iCol < 22; iCol++) {
            BOOST_TEST(twoMatrixEntanglement(iRow, iCol) < vectorEntanglement[iRow] + vectorEntanglement[iCol]);
        }
    }
}

#endif // HAS_NU1

#endif // DMRG_VIBRATIONAL
