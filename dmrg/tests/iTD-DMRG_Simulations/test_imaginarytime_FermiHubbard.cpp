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

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <boost/filesystem.hpp>

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/evolve/TimeEvolutionSweep.h"
#include "utils/fpcomparison.h"
#include "utils/io.hpp"
#include "dmrg/utils/time_stopper.h"
#include "Fixtures/TranscorrelatedFixture.h"
#include "maquis_dmrg.h"

/**
 * @brief iTD-DMRG calculation on the real-space Fermi-Hubbard model.
 * 
 * Here we check that, for the 2x2 real-space Fermi-Hubbard Hamiltonian, iTD-DMRG
 * and TI-DMRG return the same energy.
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimevsTIFermiHubbard_RealSpace2x2, TranscorrelatedFixture)
{
#ifdef HAVE_TwoU1
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("nsweeps", 10);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("max_bond_dimension", 50);
    // TI
    maquis::DMRGInterface<double> interfaceTI(parameters2x2_RealSpace_U4_2Alpha1Beta);
    interfaceTI.optimize();
    auto energyTI = interfaceTI.energy();
    // TD
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("time_step", 10.);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("propagator_maxiter", 10);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("imaginary_time", "yes");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("TD_backpropagation", "no");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("simulation_type", "TD");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("COMPLEX", 1);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("time_units", "fs");
    maquis::DMRGInterface<std::complex<double>> interfaceTD(parameters2x2_RealSpace_U4_2Alpha1Beta);
    interfaceTD.evolve();
    auto energyTD = std::real(interfaceTD.energy());
    BOOST_CHECK_CLOSE(energyTD, energyTI, 1.0E-8);
#endif // HAVE_TwoU1
}

#ifdef HAVE_TwoU1

/** 
 * @brief Test on tcDMRG with no correlation parameter.
 * 
 * Checks that tcDMRG with J=0 returns the same energy as iTD-DMRG.
 * We take as a reference the 2x2 real-space Fermi-Hubbard Hamiltonian,
 * and intentionally set the bond dimension m to 10 in order to ensure
 * that the equality holds also for a low bond dimension.
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_RealSpace2x2_Transcorrelated_J0, TranscorrelatedFixture)
{
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("nsweeps", 10);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("max_bond_dimension", 10);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("propagator_maxiter", 10);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("imaginary_time", "yes");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("TD_backpropagation", "no");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("simulation_type", "TD");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("COMPLEX", 1);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("time_units", "fs");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("time_step", 10.);
    // iTD-DMRG
    maquis::DMRGInterface<std::complex<double>> interfaceTD(parameters2x2_RealSpace_U4_2Alpha1Beta);
    interfaceTD.evolve();
    auto energyTD = std::real(interfaceTD.energy());
    // tcDMRG
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("transcorrelated_hamiltonian", "yes");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("J_Transcorrelated", 0.);
    maquis::DMRGInterface<std::complex<double>> interfaceTC(parameters2x2_RealSpace_U4_2Alpha1Beta);
    interfaceTC.evolve();
    auto energyTC = std::real(interfaceTC.energy());
    BOOST_CHECK_CLOSE(energyTD, energyTC, 1.0E-8);
}

/** 
 * @brief Checks consistency of calculations with +/- the same correlation parameter.
 * 
 * For a sufficiently large bond dimension, tcDMRG should converge to the same energy
 * independently on the correlation parameter J. We check that here by running independent 
 * calculations with J and -J, and check that the final energy is the same.
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_RealSpace2x2_Transcorrelated_OppositeJ, TranscorrelatedFixture)
{
    std::vector<double> vectorOfJValues = {0.1, 0.5, 1.0};
    for (const auto& iJ: vectorOfJValues) {
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("nsweeps", 25);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("max_bond_dimension", 50);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("init_state", "const");
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("optimization", "twosite");
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("propagator_maxiter", 10);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("imaginary_time", "yes");
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("TD_backpropagation", "no");
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("simulation_type", "TD");
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("COMPLEX", 1);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("time_units", "as");
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("time_step", 100.);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("transcorrelated_hamiltonian", "yes");
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("J_Transcorrelated", iJ);
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("hamiltonian_units", "Hartree");
        maquis::DMRGInterface<std::complex<double>> interfaceTCPlus(parameters2x2_RealSpace_U4_2Alpha1Beta);
        interfaceTCPlus.evolve();
        auto energyTCPlus = std::real(interfaceTCPlus.energy());
        //
        parameters2x2_RealSpace_U4_2Alpha1Beta.set("J_Transcorrelated", -iJ);
        maquis::DMRGInterface<std::complex<double>> interfaceTCMinus(parameters2x2_RealSpace_U4_2Alpha1Beta);
        interfaceTCMinus.evolve();
        auto energyTCMinus = std::real(interfaceTCMinus.energy());
        BOOST_CHECK_CLOSE(energyTCPlus, energyTCMinus, 1.0E-8);
    }
}

/**
 * @brief Compares the results of a conventional and transcorrelated DMRG calculation.
 * 
 * We use as a reference the two-dimensional real-space Fermi-Hubbard Hamiltonian.
 * We take a simple lattice (2x2, with 3 electrons), for which TI-DMRG and tcDMRG
 * are expected to converge to the same limit with a relatively low bond dimension
 * of m=100. We then do the following checks:
 * 
 *  1) We first calculate the energy with TI-DMRG, tcDMRG[SS], and tcDMRG[TS] and 
 *     verify that the results that we obtain is the same. Note that the reference
 *     energy is taken from PRB, 88, 125132 (2015).
 *  2) We then explicitly apply the correlator onto the MPS, and verify that the
 *     resulting MPS gives the same energy as TI-DMRG, if the energy is evaluated
 *     with the *un-transcorrelated* MPO.
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_CheckLeftVsRight, TranscorrelatedFixture)
{
    // Type definition
    using SingleSiteTimeEvolution = SingleSiteTimeEvolution<cmatrix, TwoU1, storage::disk>;
    using TwoSiteTimeEvolution = TwoSiteTimeEvolution<cmatrix, TwoU1, storage::disk>;
    // == TI-DMRG CALCULATION ==
    auto parametersTI = parameters2x2_RealSpace_U4_2Alpha1Beta;
    parametersTI.set("nsweeps", 20);
    parametersTI.set("max_bond_dimension", 100);
    maquis::DMRGInterface<double> interfaceTI(parametersTI);
    interfaceTI.optimize();
    auto energyTI = interfaceTI.energy();
    BOOST_CHECK_SMALL(energyTI - -1.60463*4, 1.0E-3);
    // == TC-DMRG CALCULATION ==
    // Note that we do the check for various J values
    std::vector<double> vectorOfTranscorrelationParameters = {-1., -0.5, -0.1, 0.1, 0.5, 1.};
    for (const auto& jValue: vectorOfTranscorrelationParameters) {
        // Setup of the parameter object
        auto parametersTCDMRG = parameters2x2_RealSpace_U4_2Alpha1Beta;
        parametersTCDMRG.set("max_bond_dimension", 100);
        parametersTCDMRG.set("init_state", "const");
        parametersTCDMRG.set("alpha_initial", 1.0E-8);
        parametersTCDMRG.set("alpha_main", 1.0E-15);
        parametersTCDMRG.set("alpha_final", 1.0E-30);
        parametersTCDMRG.set("propagator_maxiter", 10);
        parametersTCDMRG.set("propagator_accuracy", 1.0E-10);
        parametersTCDMRG.set("imaginary_time", "yes");
        parametersTCDMRG.set("TD_backpropagation", "no");
        parametersTCDMRG.set("simulation_type", "TD");
        parametersTCDMRG.set("COMPLEX", 1);
        parametersTCDMRG.set("time_units", "as");
        parametersTCDMRG.set("time_step", 10.);
        parametersTCDMRG.set("transcorrelated_hamiltonian", "yes");
        parametersTCDMRG.set("J_Transcorrelated", jValue);
        parametersTCDMRG.set("hamiltonian_units", "Hartree");
        // Setup of the simulation parameters.
        // Note that we don't use the interface in order to be able to access the MPS.
        int nSweeps = 20;
        time_stopper stop_callback(10000.);
        auto lat = Lattice(parametersTCDMRG);
        auto model = Model<cmatrix, TwoU1>(lat, parametersTCDMRG);
        auto mpo = make_mpo(lat, model);
        auto mpsLeftSS = MPS<cmatrix, TwoU1>(lat.size(), *(model.initializer(lat, parametersTCDMRG)));
        auto mpsLeftTS = mpsLeftSS;
        // Single-site time evolution
        auto ssEvolverLeft = SingleSiteTimeEvolution(mpsLeftSS, mpo, parametersTCDMRG, stop_callback);
        for (int iSweep = 0; iSweep < nSweeps; iSweep++)
            ssEvolverLeft.evolve_sweep(iSweep);
        auto energySS = std::real(expval(mpsLeftSS, mpo));
        BOOST_CHECK_CLOSE(energySS, energyTI, 1.0E-8);
        // Two-site time evolution
        auto tsEvolverLeft = TwoSiteTimeEvolution(mpsLeftTS, mpo, parametersTCDMRG, stop_callback);
        for (int iSweep = 0; iSweep < nSweeps; iSweep++)
            tsEvolverLeft.evolve_sweep(iSweep);
        auto energyTS = std::real(expval(mpsLeftTS, mpo));
        BOOST_CHECK_CLOSE(energyTS, energyTI, 1.0E-8);
        // == CHECKS TRANSCORRELATION PROPERTY ==
        // We first construct proper
        auto latOriginal = Lattice(parameters2x2_RealSpace_U4_2Alpha1Beta);
        auto modelOriginal = Model<cmatrix, TwoU1>(latOriginal, parameters2x2_RealSpace_U4_2Alpha1Beta);
        auto mpoOriginal = make_mpo(latOriginal, modelOriginal);
        // Single-site
        auto mpsTimesMPOSS = mpsLeftSS;
        for (int iSite = 0; iSite < mpsTimesMPOSS.size(); iSite++)
            mpsTimesMPOSS[iSite].scaleByExponentialProductOfCharges(-jValue);
        // Two-site
        auto mpsTimesMPOTS = mpsLeftTS;
        for (int iSite = 0; iSite < mpsTimesMPOTS.size(); iSite++)
            mpsTimesMPOTS[iSite].scaleByExponentialProductOfCharges(-jValue);
        auto energySSFromExponential = std::real(expval(mpsTimesMPOSS, mpoOriginal)/overlap(mpsTimesMPOSS, mpsTimesMPOSS));
        auto energyTSFromExponential = std::real(expval(mpsTimesMPOTS, mpoOriginal)/overlap(mpsTimesMPOTS, mpsTimesMPOTS));
        BOOST_CHECK_CLOSE(energySSFromExponential, energyTI, 1.0E-8);
        BOOST_CHECK_CLOSE(energyTSFromExponential, energyTI, 1.0E-8);
    }
}

/**
 * @brief iTD-DMRG calculation on the momentum-space Fermi-Hubbard model.
 * 
 * Here we check that, for the 2x2 momentum-space Fermi-Hubbard Hamiltonian, iTD-DMRG
 * and TI-DMRG return the same energy.
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimevsTIFermiHubbard_Momentum2x2, TranscorrelatedFixture)
{
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("nsweeps", 10);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("max_bond_dimension", 300);
    // TI
    maquis::DMRGInterface<double> interfaceTI(parameters2x2_MomentumSpace_U4_2Alpha1Beta);
    interfaceTI.optimize();
    auto energyTI = interfaceTI.energy();
    // TD
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("time_step", 10.);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("propagator_maxiter", 10);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("imaginary_time", "yes");
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("TD_backpropagation", "no");
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("simulation_type", "TD");
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("COMPLEX", 1);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("time_units", "fs");
    maquis::DMRGInterface<std::complex<double>> interfaceTD(parameters2x2_MomentumSpace_U4_2Alpha1Beta);
    interfaceTD.evolve();
    auto energyTD = std::real(interfaceTD.energy());
    BOOST_CHECK_CLOSE(energyTD, energyTI, 1.0E-8);
}

/** 
 * @brief Checks consistency of calculations with +/- the same correlation parameter.
 * 
 * For a sufficiently large bond dimension, tcDMRG should converge to the same energy
 * independently on the correlation parameter J. We check that here by running independent 
 * calculations with different J values, and check that the final energy is the same.
 * Here we use as a reference the *momentum*-space Fermi-Hubbard Hamiltonian
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_RealSpace2x2_Transcorrelated_DifferentJ, TranscorrelatedFixture)
{
    std::vector<std::pair<double, double>> jValues = {std::make_pair(0.3, -0.2),
                                                      std::make_pair(0.5, -0.1),
                                                      std::make_pair(0.1, -1.2)};
    for (const auto& iJPair: jValues) {
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("nsweeps", 25);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("max_bond_dimension", 50);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("init_state", "const");
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("optimization", "twosite");
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("propagator_maxiter", 10);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("imaginary_time", "yes");
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("TD_backpropagation", "no");
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("simulation_type", "TD");
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("COMPLEX", 1);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("time_units", "as");
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("time_step", 100.);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("transcorrelated_hamiltonian", "yes");
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("J_Transcorrelated", iJPair.first);
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("hamiltonian_units", "Hartree");
        maquis::DMRGInterface<std::complex<double>> interfaceTCPlus(parameters2x2_MomentumSpace_U4_2Alpha1Beta);
        interfaceTCPlus.evolve();
        auto energyTCPlus = std::real(interfaceTCPlus.energy());
        //
        parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("J_Transcorrelated", iJPair.second);
        maquis::DMRGInterface<std::complex<double>> interfaceTCMinus(parameters2x2_MomentumSpace_U4_2Alpha1Beta);
        interfaceTCMinus.evolve();
        auto energyTCMinus = std::real(interfaceTCMinus.energy());
        BOOST_CHECK_CLOSE(energyTCPlus, energyTCMinus, 1.0E-8);
    }
}

#endif // HAVE_TwoU1