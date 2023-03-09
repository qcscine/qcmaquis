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


#define BOOST_TEST_MODULE TranscorrelatedNO

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "Fixtures/TranscorrelatedFixture.h"


/**
 * @brief Checks consistency between transcorrelated with three body and normal-ordered transcorrelated
 *         implementations.
 *
 * Energies are calculated for both the transcorrelated Hamiltonian with full three body interactions and the
 * normal-ordered Hamiltonian. These should give equal results.
 */
//BOOST_FIXTURE_TEST_CASE(TestTCMolecular_HE_NormalOrdered, TranscorrelatedFixture)
//{
//    parametersHeTranscorrelated.set("truncation_initial", 1e-50);
//    parametersHeTranscorrelated.set("truncation_final", 1e-50);
//    parametersHeTranscorrelated.set("init_state", "default");
//    parametersHeTranscorrelated.set("site_types", "0,0,0,0,0");
//    parametersHeTranscorrelated.set("max_bond_dimension", 100);
//    parametersHeTranscorrelated.set("propagator_accuracy", 1.0E-10);
//    parametersHeTranscorrelated.set("propagator_maxiter", 10);
//    parametersHeTranscorrelated.set("nsweeps", 0);
//    parametersHeTranscorrelated.set("hamiltonian_units", "Hartree");
//    parametersHeTranscorrelated.set("time_units", "fs");
//    parametersHeTranscorrelated.set("TD_backpropagation", "no");
//    parametersHeTranscorrelated.set("symmetry", "2u1");
//    parametersHeTranscorrelated.set("time_step", 0.1);
//    parametersHeTranscorrelated.set("optimization", "singlesite");
//    parametersHeTranscorrelated.set("transcorrelated_nsweeps_TI", 0);
//    parametersHeTranscorrelated.set("transcorrelated_nsweeps_TC", 20);
//    parametersHeTranscorrelated.set("integral_file", CONVENTIONAL_HE_FCIDUMP_PATH);
//    parametersHeTranscorrelated.set("transcorrelated_integral_file", TRANSCORRELATED_HE_FCIDUMP_PATH);
//
//    //auto energyDMRG = maquis::real(interface.energy());
//
//    parametersHeTranscorrelated.set("transcorrelated_3body_normal_ordered", "yes");
//    maquis::DMRGInterface<double> NO_interface(parametersHeTranscorrelated);
//    NO_interface.runTranscorrelated();
//
//    auto energyDMRGNO = maquis::real(NO_interface.energy());
//
//    parametersHeTranscorrelated.set("transcorrelated_3body_normal_ordered", "no");
//    maquis::DMRGInterface<double> interface(parametersHeTranscorrelated);
//    interface.runTranscorrelated();
//
//    auto energyDMRG = maquis::real(interface.energy());
//
//    BOOST_CHECK_SMALL(std::abs(energyDMRG-energyDMRGNO), 1.0E-10);
//
//    maquis::cout << "Energies: " << energyDMRG << " " << energyDMRGNO  << std::endl;
//}

/**
 * @brief Transcorrelated Normal Ordered DMRG calculation on Be.
 * The reference energy has been generated in this case with the Owl CC
 * code by Max Moerchen.
 */
BOOST_FIXTURE_TEST_CASE(TestTCMolecular_Be_VersusCCNO, TranscorrelatedFixture)
{
    parametersBeTranscorrelatedTwoBody.set("propagator_accuracy", 1.0E-10);
    parametersBeTranscorrelatedTwoBody.set("propagator_maxiter", 10);
    parametersBeTranscorrelatedTwoBody.set("hamiltonian_units", "Hartree");
    parametersBeTranscorrelatedTwoBody.set("time_units", "fs");
    parametersBeTranscorrelatedTwoBody.set("TD_backpropagation", "no");
    parametersBeTranscorrelatedTwoBody.set("symmetry", "2u1");
    parametersBeTranscorrelatedTwoBody.set("time_step", 0.1);
    parametersBeTranscorrelatedTwoBody.set("chkpfile", "Be.tcDMRGNO.checkpoint.h5");
    parametersBeTranscorrelatedTwoBody.set("transcorrelated_nsweeps_TI", 0);
    parametersBeTranscorrelatedTwoBody.set("transcorrelated_nsweeps_TC", 10);
    parametersBeTranscorrelatedTwoBody.set("integral_file", "IntegralFile_Be_Conventional");
    parametersBeTranscorrelatedTwoBody.set("transcorrelated_integral_file", TRANSCORRELATED_BE_FCIDUMP_PATH);
    parametersBeTranscorrelatedTwoBody.set("transcorrelated_3body_normal_ordered", "yes");
    parametersBeTranscorrelatedTwoBody.set("transcorrelated_3body", "yes");
    maquis::DMRGInterface<double> interface(parametersBeTranscorrelatedTwoBody);
    interface.runTranscorrelated();
    auto energyDMRG = maquis::real(interface.energy());

    parametersBeTranscorrelatedTwoBody.set("transcorrelated_3body", "no");
    maquis::DMRGInterface<double> approxInterface(parametersBeTranscorrelatedTwoBody);
    approxInterface.runTranscorrelated();
    auto energyDRMGno3b = maquis::real(approxInterface.energy());

    parametersBeTranscorrelatedTwoBody.set("transcorrelated_3body_normal_ordered", "no");
    parametersBeTranscorrelatedTwoBody.set("transcorrelated_3body", "yes");
    maquis::DMRGInterface<double> tcInterface(parametersBeTranscorrelatedTwoBody);
    tcInterface.runTranscorrelated();
    auto energytcDMRG = maquis::real(tcInterface.energy());

    auto refEnergy = -14.6505807967243;
    BOOST_CHECK_SMALL(std::abs(energyDMRG-energytcDMRG), 1.0E-9);

    std::cout << "Difference with and withour 3B: " << std::abs(energyDMRG - energyDRMGno3b) << std::endl;
}