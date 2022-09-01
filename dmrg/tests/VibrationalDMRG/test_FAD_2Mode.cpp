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

#ifdef DMRG_VIBRATIONAL

#include <random>
#include <boost/filesystem.hpp>
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
BOOST_FIXTURE_TEST_CASE(Test_vDMRG_FAD_2ModeHamiltonian, NModeFixture)
{
#ifdef HAVE_NU1
    // Adds the final input parameters
    parametersFADTwoBody.set("init_state", "const");
    parametersFADTwoBody.set("nsweeps", 20);
    parametersFADTwoBody.set("max_bond_dimension", 100);
    parametersFADTwoBody.set("MODEL", "nmode");
    parametersFADTwoBody.set("truncation_initial", 1.0E-20);
    parametersFADTwoBody.set("truncation_final", 1.0E-16);
    // SingleSite
    parametersFADTwoBody.set("optimization", "singlesite");
    parametersFADTwoBody.set("alpha_initial", 1.0E-8);
    parametersFADTwoBody.set("alpha_initial", 1.0E-15);
    parametersFADTwoBody.set("alpha_initial", 0.);
    maquis::DMRGInterface<double> interfaceSS(parametersFADTwoBody);
    interfaceSS.optimize();
    BOOST_CHECK_CLOSE(interfaceSS.energy(), -1499.5871477508479, 1.0E-5);
    // TwoSite
    parametersFADTwoBody.set("optimization", "twosite");
    maquis::DMRGInterface<double> interfaceTS(parametersFADTwoBody);
    interfaceTS.optimize();
    BOOST_CHECK_CLOSE(interfaceTS.energy(), -1499.5871477508479, 1.0E-5);
#endif // HAS_NU1
}

#ifdef HAVE_NU1

/**
 * @brief N-mode calculation on one-mode PES of FAD.
 * The basis set has been generated with a DVR primitive basis.
 * Note that we also check the excitation energy.
 */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_1ModeHamiltonian_ExcitedState, NModeFixture)
{
    // Adds the final input parameters
    parametersFADOneBodyBinary.set("init_state", "const");
    parametersFADOneBodyBinary.set("nsweeps", 20);
    parametersFADOneBodyBinary.set("max_bond_dimension",100);
    parametersFADOneBodyBinary.set("MODEL", "nmode");
    parametersFADOneBodyBinary.set("chkpfile", "GS.checkpoint.h5");
    parametersFADOneBodyBinary.set("resfule", "GS.results.h5");
    // Creates the interface and checks the resulting energy
    maquis::DMRGInterface<double> interface(parametersFADOneBodyBinary);
    interface.optimize();
    auto gsEnergy = interface.energy();
    BOOST_CHECK_CLOSE(gsEnergy, -2.359242429009664e+03, 1.0E-5);
    // Checks that the overlap of the final wave function with the hf determinant is = 1.
    auto targetOverlap = interface.getCICoefficient("0");
    BOOST_CHECK_CLOSE(std::abs(targetOverlap), 1.0, 1.0E-16);
    // Checks that overlap calculations cannot be performed for invalid reference determinants
    BOOST_CHECK_THROW(interface.getCICoefficient("0,0,0"), std::runtime_error);
    // Excited-state calculations
    parametersFADOneBodyBinary.set("chkpfile", "ES.checkpoint.h5");
    parametersFADOneBodyBinary.set("resfule", "ES.results.h5");
    parametersFADOneBodyBinary.set("n_ortho_states", 1);
    parametersFADOneBodyBinary.set("ortho_states", "GS.checkpoint.h5");
    // Creates a new interface object and reruns the optimization
    maquis::DMRGInterface<double> interfaceES(parametersFADOneBodyBinary);
    interfaceES.optimize();
    auto esEnergy = interfaceES.energy();
    BOOST_CHECK_CLOSE(esEnergy-gsEnergy, 0.4450425713052937, 1.0E-5);
    boost::filesystem::remove_all("GS.results.h5");
    boost::filesystem::remove_all("GS.checkpoint.h5");
    boost::filesystem::remove_all("ES.results.h5");
    boost::filesystem::remove_all("ES.checkpoint.h5");
}

/** @brief Check that the ground-state energy is the same irrespectively on the modal order */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_FingerprintHamiltonian_Sorting, NModeFixture)
{
    // Adds the final input parameters
    parametersFADTwoBodyFingerPrint.set("init_state", "const");
    parametersFADTwoBodyFingerPrint.set("nsweeps", 20);
    parametersFADTwoBodyFingerPrint.set("max_bond_dimension", 100);
    parametersFADTwoBodyFingerPrint.set("twosite_truncation", "heev_truncate");
    parametersFADTwoBodyFingerPrint.set("alpha_initial", 1.0E-8);
    parametersFADTwoBodyFingerPrint.set("alpha_main", 1.0E-15);
    parametersFADTwoBodyFingerPrint.set("alpha_final", 0.);
    parametersFADTwoBodyFingerPrint.set("ngrowsweeps", 2);
    parametersFADTwoBodyFingerPrint.set("nmainsweeps", 2);
    // Creates the interface and checks the resulting energy
    maquis::DMRGInterface<double> interfaceConventionalSorting(parametersFADTwoBodyFingerPrint);
    interfaceConventionalSorting.optimize();
    auto conventionalSortingEnergy = interfaceConventionalSorting.energy();
    // Generates randomly a reshuffling
    int latticeSize = parametersFADTwoBodyFingerPrint["L"];
    std::vector<int> newOrder(latticeSize);
    for (int iSite = 0; iSite < latticeSize; iSite++)
        newOrder[iSite] = iSite;
    std::shuffle(newOrder.begin(), newOrder.end(), std::default_random_engine());
    std::string inputOrder = "";
    for (int iElement = 0; iElement < newOrder.size(); iElement++) {
        inputOrder += std::to_string(newOrder[iElement]);
        if (iElement != newOrder.size()-1)
            inputOrder += ",";
    }
    parametersFADTwoBodyFingerPrint.set("modals_order", inputOrder);
    maquis::DMRGInterface<double> interfaceRandomSorting(parametersFADTwoBodyFingerPrint);
    interfaceRandomSorting.optimize();
    auto randomSortingEnergy = interfaceRandomSorting.energy();
    BOOST_CHECK_CLOSE(conventionalSortingEnergy, randomSortingEnergy, 1.0E-7);
}

/** @brief Same as above, but includes measurements */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_2ModeSystem_MeasOneParticle, NModeFixture)
{
    // Adds the final input parameters
    parametersFADTwoBody.set("init_state", "const");
    parametersFADTwoBody.set("seed", 16071991);
    parametersFADTwoBody.set("nsweeps", 100);
    parametersFADTwoBody.set("max_bond_dimension", 20);
    parametersFADTwoBody.set("MODEL", "nmode");
    parametersFADTwoBody.set("MEASURE[One Modal RDM]", "1");
    // Creates the interface
    maquis::DMRGInterface<double> interface(parametersFADTwoBody);
    interface.optimize();
    BOOST_CHECK_CLOSE(interface.energy(), -1499.5871477508479, 1.0E-5);
    // Measurements
    interface.measure();
    std::vector<double> firstMeas(22, 0);
    std::vector<double> secondMeas(22, 0);
    for (int iSite = 0; iSite < 22; iSite++) {
        auto idx1 = interface.getMeasurement("onemodalRDM_00").first[iSite][0];
        auto idx2 = interface.getMeasurement("onemodalRDM_11").first[iSite][0];
        firstMeas[idx1] = interface.getMeasurement("onemodalRDM_00").second[iSite];
        secondMeas[idx2] = interface.getMeasurement("onemodalRDM_11").second[iSite];
    }
    // Final check
    for (int iSite = 0; iSite < 22; iSite++)
        BOOST_CHECK_CLOSE(firstMeas[iSite]+secondMeas[iSite], 1., 1.0E-16);
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
    maquis::DMRGInterface<double> interface(parametersFADTwoBody);
    interface.optimize();
    interface.measure();
    // == TWO-MODAL ENTROPY ==
    tmatrix<double> twoMatrixEntanglement(22, 22, 0.), twoMode00(22, 22, 0.), twoMode11(22, 22, 0.),
        twoMode22(22, 22, 0.), twoMode33(22, 22, 0.), twoMode12(22, 22, 0.), twoMode21(22, 22, 0.);
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
    // Extracts the two-modal entropy and calculates the entanglement
    auto overallSize = interface.getMeasurement("twomodeRDM_00").first.size();
    for (int iElement = 0; iElement < overallSize; iElement++) {
        // 0-0
        int iRow = interface.getMeasurement("twomodeRDM_00").first[iElement][0];
        int iCol = interface.getMeasurement("twomodeRDM_00").first[iElement][1];
        twoMode00(iRow, iCol) =  interface.getMeasurement("twomodeRDM_00").second[iElement];
        // 1-1
        iRow = interface.getMeasurement("twomodeRDM_11").first[iElement][0];
        iCol = interface.getMeasurement("twomodeRDM_11").first[iElement][1];
        twoMode11(iRow, iCol) =  interface.getMeasurement("twomodeRDM_11").second[iElement];
        // 1-2
        iRow = interface.getMeasurement("twomodeRDM_12").first[iElement][0];
        iCol = interface.getMeasurement("twomodeRDM_12").first[iElement][1];
        twoMode12(iRow, iCol) =  interface.getMeasurement("twomodeRDM_12").second[iElement];
        // 2-1
        iRow = interface.getMeasurement("twomodeRDM_21").first[iElement][0];
        iCol = interface.getMeasurement("twomodeRDM_21").first[iElement][1];
        twoMode21(iRow, iCol) =  interface.getMeasurement("twomodeRDM_21").second[iElement];
        // 2-2
        iRow = interface.getMeasurement("twomodeRDM_22").first[iElement][0];
        iCol = interface.getMeasurement("twomodeRDM_22").first[iElement][1];
        twoMode22(iRow, iCol) =  interface.getMeasurement("twomodeRDM_22").second[iElement];
        // 3-3
        iRow = interface.getMeasurement("twomodeRDM_33").first[iElement][0];
        iCol = interface.getMeasurement("twomodeRDM_33").first[iElement][1];
        twoMode33(iRow, iCol) =  interface.getMeasurement("twomodeRDM_33").second[iElement];
    }
    // Calculates the entanglement
    for (int iRow = 0; iRow < 22; iRow++) {
        for (int iCol = 0; iCol < 22; iCol++) {
            // Constructs the two-modal entanglement matrix
            tmatrix<double> entanglementMatrix(4, 4, 0.), evecs(4, 4, 0.);
            std::vector<double> evals(4, 0.);
            entanglementMatrix(0, 0) = twoMode00(iRow, iCol);
            entanglementMatrix(1, 1) = twoMode11(iRow, iCol);
            entanglementMatrix(1, 2) = twoMode12(iRow, iCol);
            entanglementMatrix(2, 1) = twoMode21(iRow, iCol);
            entanglementMatrix(2, 2) = twoMode22(iRow, iCol);
            entanglementMatrix(3, 3) = twoMode33(iRow, iCol);
            heev(entanglementMatrix, evecs, evals);
            // Calculates the two-orbital entropy
            double entropy = 0.;
            for (int i = 0; i < 4; i++)
                if (std::abs(evals[i]) > 1.0E-16)
                    entropy -= evals[i]*std::log(evals[i]);
            twoMatrixEntanglement(iRow, iCol) = entropy;
        }
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

/**
 * @brief Measure the one-mode RDM and checks sanity of the trace.
 * Note that the trace of the one-mode RDM should be 1 for each sub-block of modals
 * belonging to the same mode. The overall trace should then be equal to the overall number
 * of modes.
 */
BOOST_FIXTURE_TEST_CASE(Test_Lattice_Size_2ModeSystem_ModeRDM, NModeFixture)
{
    // Adds the final input parameters
    parametersFADTwoBody.set("init_state", "default");
    parametersFADTwoBody.set("seed", 30031989);
    parametersFADTwoBody.set("nsweeps", 10);
    parametersFADTwoBody.set("max_bond_dimension", 20);
    parametersFADTwoBody.set("MODEL", "nmode");
    parametersFADTwoBody.set("MEASURE[One Mode RDM]", "1");
    // Creates the interface
    maquis::DMRGInterface<double> interface(parametersFADTwoBody);
    interface.optimize();
    interface.measure();
    // == TWO-MODAL ENTROPY ==
    tmatrix<double> oneModeRDM(22, 22, 0.);
    auto overallSize = interface.getMeasurement("onemodeRDM").first.size();
    for (int iElement = 0; iElement < overallSize; iElement++) {
        // 0-0
        int iRow = interface.getMeasurement("onemodeRDM").first[iElement][0];
        int iCol = interface.getMeasurement("onemodeRDM").first[iElement][1];
        oneModeRDM(iRow, iCol) =  interface.getMeasurement("onemodeRDM").second[iElement];
    }
    // Calculates the overall trace
    double trace = 0.;
    for (int idx = 0; idx < 22; idx++)
        trace += oneModeRDM(idx, idx);
    BOOST_CHECK_CLOSE(trace, 2., 1.0E-10);
    // Calculates the mode-resolved trace
    double traceMode1 = 0., traceMode2 = 0.;
    for (int idx = 0; idx < 11; idx++) {
        traceMode1 += oneModeRDM(idx, idx);
        traceMode2 += oneModeRDM(idx+11, idx+11);
    }
    BOOST_CHECK_CLOSE(traceMode1, 1., 1.0E-10);
    BOOST_CHECK_CLOSE(traceMode2, 1., 1.0E-10);
    // Checks that the off-diagonal elements of the one-mode RDM are zero.
    for (int iRow = 0; iRow < 11; iRow++) {
        for (int iCol = 11; iCol < 22; iCol++) {
            BOOST_CHECK_CLOSE(oneModeRDM(iRow, iCol), 0., 1.0E-10);
            BOOST_CHECK_CLOSE(oneModeRDM(iCol, iRow), 0., 1.0E-10);
        }
    }
}

#endif // HAS_NU1

#endif // DMRG_VIBRATIONAL
