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

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/utils/time_stopper.h"
#include "dmrg/evolve/TimeEvolutionSweep.h"
#include "utils/fpcomparison.h"
#include "utils/io.hpp"
#include "Fixtures/TranscorrelatedFixture.h"
#include "maquis_dmrg.h"

/**
 * @brief iTD-DMRG calculation on the real-space Fermi-Hubbard model.
 *
 * Here we check that, for the 2x2 real-space Fermi-Hubbard Hamiltonian, iTD-DMRG
 * and TI-DMRG return the same energy (nothing here is transcorrelated).
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimevsTIFermiHubbard_RealSpace2x2, TranscorrelatedFixture)
{
#if defined(HAVE_TwoU1) and defined(DMRG_TD)
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("nsweeps", 10);
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("max_bond_dimension", 50);
  // TI-DMRG
  maquis::DMRGInterface<double> interfaceTI(parameters2x2_RealSpace_U4_2Alpha1Beta);
  interfaceTI.optimize();
  auto energyTI = interfaceTI.energy();
  // iTD-DMRG
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

/**
 * @brief Test on tcDMRG with no correlation parameter.
 *
 * Checks that tcDMRG with J=0 returns the same energy as iTD-DMRG.
 * We take as a reference the 2x2 real-space Fermi-Hubbard Hamiltonian,
 * and intentionally set the bond dimension m to 10 (i.e., to a low value)
 * in order to ensure that the equality holds also for a low bond dimension.
 *
 * Note that the parameters of iTD-DMRG and tcDMRG are the same, so the
 * simulation is effectively the same.
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_RealSpace2x2_Transcorrelated_J0, TranscorrelatedFixture)
{
#if defined(HAVE_TwoU1) and defined(DMRG_TC)
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("nsweeps", 10);
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("max_bond_dimension", 10);
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("propagator_maxiter", 10);
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("imaginary_time", "yes");
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("TD_backpropagation", "no");
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("time_units", "fs");
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("time_step", 10.);
  // iTD-DMRG
  maquis::DMRGInterface<std::complex<double>> interfaceTD(parameters2x2_RealSpace_U4_2Alpha1Beta);
  interfaceTD.evolve();
  auto energyTD = std::real(interfaceTD.energy());
  // tcDMRG
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("transcorrelated_nsweeps_TI", 0);
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("transcorrelated_nsweeps_TC", 10);
  parameters2x2_RealSpace_U4_2Alpha1Beta.set("J_Transcorrelated", 0.);
  maquis::DMRGInterface<double> interfaceTC(parameters2x2_RealSpace_U4_2Alpha1Beta);
  interfaceTC.runTranscorrelated();
  auto energyTC = std::real(interfaceTC.energy());
  BOOST_CHECK_CLOSE(energyTD, energyTC, 1.0E-8);
#endif // HAVE_TwoU1 and DMRG_TC
}

#if defined(HAVE_TwoU1) and defined(DMRG_TC)

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
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("transcorrelated_nsweeps_TI", 10);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("transcorrelated_nsweeps_TC", 25);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("max_bond_dimension", 50);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("init_state", "const");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("optimization", "twosite");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("propagator_maxiter", 10);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("TD_backpropagation", "no");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("time_units", "as");
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("time_step", 100.);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("J_Transcorrelated", iJ);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("hamiltonian_units", "Hartree");
    maquis::DMRGInterface<std::complex<double>> interfaceTCPlus(parameters2x2_RealSpace_U4_2Alpha1Beta);
    interfaceTCPlus.runTranscorrelated();
    auto energyTCPlus = std::real(interfaceTCPlus.energy());
    //
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("J_Transcorrelated", -iJ);
    maquis::DMRGInterface<std::complex<double>> interfaceTCMinus(parameters2x2_RealSpace_U4_2Alpha1Beta);
    interfaceTCMinus.runTranscorrelated();
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
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_RealSpace2x2_Transcorrelated_CheckLeftVsRight, TranscorrelatedFixture)
{
  using MPSType = MPS<matrix, TwoU1>;
  MPSType mpsSS, mpsTS;
  // == TI-DMRG CALCULATION ==
  auto parametersTI = parameters2x2_RealSpace_U4_2Alpha1Beta;
  parametersTI.set("nsweeps", 20);
  parametersTI.set("max_bond_dimension", 100);
  maquis::DMRGInterface<double> interfaceTI(parametersTI);
  interfaceTI.optimize();
  auto energyTI = interfaceTI.energy();
  // The reference energy is reported per site, so we must multiply by 4.
  BOOST_CHECK_SMALL(energyTI - -1.60463*4, 1.0E-3);
  // == TC-DMRG CALCULATION ==
  // Setup of the parameter object
  auto parametersTCDMRG = parameters2x2_RealSpace_U4_2Alpha1Beta;
  parametersTCDMRG.set("max_bond_dimension", 100);
  parametersTCDMRG.set("alpha_initial", 1.0E-8);
  parametersTCDMRG.set("alpha_main", 1.0E-15);
  parametersTCDMRG.set("alpha_final", 1.0E-30);
  parametersTCDMRG.set("propagator_maxiter", 10);
  parametersTCDMRG.set("propagator_accuracy", 1.0E-10);
  parametersTCDMRG.set("TD_backpropagation", "no");
  parametersTCDMRG.set("time_units", "as");
  parametersTCDMRG.set("time_step", 10.);
  parametersTCDMRG.set("hamiltonian_units", "Hartree");
  parametersTCDMRG.set("transcorrelated_nsweeps_TI", 0);
  parametersTCDMRG.set("transcorrelated_nsweeps_TC", 30);
  // Note that we do the check for various J values
  std::vector<double> vectorOfTranscorrelationParameters = {-1., -0.5, -0.1, 0.1, 0.5, 1.};
  for (const auto& jValue: vectorOfTranscorrelationParameters) {
    parametersTCDMRG.set("J_Transcorrelated", jValue);
    // Single-site time evolution
    parametersTCDMRG.set("chkpfile", "SS.FermiHubbard.chkp.h5");
    parametersTCDMRG.set("optimizer", "singlesite");
    maquis::DMRGInterface<double> interfaceTC_SS(parametersTCDMRG);
    interfaceTC_SS.runTranscorrelated();
    auto energyTC_SS = interfaceTC_SS.energy();
    BOOST_CHECK_CLOSE(energyTC_SS, energyTI, 1.0E-8);
    // Two-site time evolution
    parametersTCDMRG.set("chkpfile", "TS.FermiHubbard.chkp.h5");
    parametersTCDMRG.set("optimizer", "twosite");
    maquis::DMRGInterface<double> interfaceTC_TS(parametersTCDMRG);
    interfaceTC_TS.runTranscorrelated();
    auto energyTC_TS = interfaceTC_TS.energy();
    BOOST_CHECK_CLOSE(energyTC_TS, energyTI, 1.0E-8);
    // == CHECKS TRANSCORRELATION PROPERTY ==
    auto latOriginal = Lattice(parameters2x2_RealSpace_U4_2Alpha1Beta);
    auto modelOriginal = Model<matrix, TwoU1>(latOriginal, parameters2x2_RealSpace_U4_2Alpha1Beta);
    auto mpoOriginal = make_mpo(latOriginal, modelOriginal);
    // Single-site
    load("SS.FermiHubbard.chkp.h5", mpsSS);
    boost::filesystem::remove_all("SS.FermiHubbard.chkp.h5");
    for (int iSite = 0; iSite < mpsSS.size(); iSite++)
      mpsSS[iSite].scaleByExponentialProductOfCharges(jValue);
    // Two-site
    load("TS.FermiHubbard.chkp.h5", mpsTS);
    boost::filesystem::remove_all("TS.FermiHubbard.chkp.h5");
    for (int iSite = 0; iSite < mpsTS.size(); iSite++)
      mpsTS[iSite].scaleByExponentialProductOfCharges(jValue);
    auto energySSFromExponential = std::real(expval(mpsSS, mpoOriginal)/overlap(mpsSS, mpsSS));
    auto energyTSFromExponential = std::real(expval(mpsTS, mpoOriginal)/overlap(mpsTS, mpsTS));
    BOOST_CHECK_CLOSE(energySSFromExponential, energyTI, 1.0E-8);
    BOOST_CHECK_CLOSE(energyTSFromExponential, energyTI, 1.0E-8);
  }
}

/**
 * @brief Compares the results of a conventional and transcorrelated DMRG calculation.
 *
 * We use as a reference the *asymmetric* two-dimensional real-space Fermi-Hubbard Hamiltonian.
 * We take a simple lattice (2x2, with 3 electrons), for tcDMRG is expected to converge to the same
 * limit with a relatively low bond dimension of m=100 for any J value.
 * We then explicitly apply e^{2*J} onto the left eigenvector, and verify that the
 * resulting MPS is equivalent to the right one.
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_AsymmetricRealSpace2x2_Transcorrelated_CheckLeftVsRight, TranscorrelatedFixture)
{
  // == J=0.5 ==
  auto parameters = parameters2x2_AsymmetricRealSpace_U4_2Alpha1Beta;
  parameters.set("transcorrelated_nsweeps_TI", 0);
  parameters.set("transcorrelated_nsweeps_TC", 20);
  parameters.set("propagator_maxiter", 10);
  parameters.set("propagator_accuracy", 1.0E-10);
  parameters.set("imaginary_time", "yes");
  parameters.set("TD_backpropagation", "no");
  parameters.set("time_units", "as");
  parameters.set("time_step", 10.);
  parameters.set("J_Transcorrelated", 0.5);
  parameters.set("hamiltonian_units", "Hartree");
  parameters.set("chkpfile", "JPlus.checkpoint.h5");
  maquis::DMRGInterface<std::complex<double>> interfaceJPlus(parameters);
  interfaceJPlus.runTranscorrelated();
  auto energyPlus = maquis::real(interfaceJPlus.energy());
  MPS<cmatrix, TwoU1> mpsJPlus;
  load("JPlus.checkpoint.h5", mpsJPlus);
  // == J=-0.5 ==
  parameters.set("J_Transcorrelated", -0.5);
  parameters.set("chkpfile", "JMinus.checkpoint.h5");
  maquis::DMRGInterface<std::complex<double>> interfaceJMinus(parameters);
  interfaceJMinus.runTranscorrelated();
  auto energyMinus = maquis::real(interfaceJMinus.energy());
  MPS<cmatrix, TwoU1> mpsJMinus;
  load("JMinus.checkpoint.h5", mpsJMinus);
  // Cleans up tmp files
  boost::filesystem::remove_all("JPlus.checkpoint.h5");
  boost::filesystem::remove_all("JMinus.checkpoint.h5");
  // First checks that the energies are the same
  BOOST_CHECK_CLOSE(energyPlus, energyMinus, 1.0E-8);
  // Now applies the correlator and checks that the final overlap is the same
  for (int iSite = 0; iSite < mpsJPlus.size(); iSite++) {
    mpsJPlus[iSite].scaleByExponentialProductOfCharges(0.5);
    mpsJMinus[iSite].scaleByExponentialProductOfCharges(-0.5);
  }
  mpsJPlus[0] /= std::sqrt(norm(mpsJPlus));
  mpsJMinus[0] /= std::sqrt(norm(mpsJMinus));
  // Final checl
  auto similiarity = overlap(mpsJMinus, mpsJPlus);
  BOOST_CHECK_CLOSE(std::abs(similiarity), 1.0, 1.0E-8);
}

/**
 * @brief iTD-DMRG calculation on the momentum-space Fermi-Hubbard model.
 *
 * Here we check that, for the 2x2 momentum-space Fermi-Hubbard Hamiltonian, iTD-DMRG
 * and TI-DMRG return the same energy.
 * This is the first test where we look into the momentun-space FH Hamiltonian.
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimevsTIFermiHubbard_Momentum2x2, TranscorrelatedFixture)
{
  // TI-DMRG
  parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("nsweeps", 10);
  parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("max_bond_dimension", 300);
  maquis::DMRGInterface<double> interfaceTI(parameters2x2_MomentumSpace_U4_2Alpha1Beta);
  interfaceTI.optimize();
  auto energyTI = interfaceTI.energy();
  // TC-DMRG
  parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("transcorrelated_nsweeps_TI", 10);
  parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("transcorrelated_nsweeps_TC", 30);
  parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("time_step", 10.);
  parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("propagator_maxiter", 10);
  parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("TD_backpropagation", "no");
  parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("time_units", "fs");
  parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("J_Transcorrelated", 0.5);
  maquis::DMRGInterface<std::complex<double>> interfaceTD(parameters2x2_MomentumSpace_U4_2Alpha1Beta);
  interfaceTD.runTranscorrelated();
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
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("transcorrelated_nsweeps_TI", 0);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("transcorrelated_nsweeps_TC", 25);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("max_bond_dimension", 50);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("init_state", "const");
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("optimization", "twosite");
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("propagator_maxiter", 10);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("imaginary_time", "yes");
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("TD_backpropagation", "no");
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("time_units", "as");
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("time_step", 100.);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("J_Transcorrelated", iJPair.first);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("hamiltonian_units", "Hartree");
    maquis::DMRGInterface<double> interfaceTCPlus(parameters2x2_MomentumSpace_U4_2Alpha1Beta);
    interfaceTCPlus.runTranscorrelated();
    auto energyTCPlus = std::real(interfaceTCPlus.energy());
    //
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("J_Transcorrelated", iJPair.second);
    maquis::DMRGInterface<double> interfaceTCMinus(parameters2x2_MomentumSpace_U4_2Alpha1Beta);
    interfaceTCMinus.runTranscorrelated();
    auto energyTCMinus = std::real(interfaceTCMinus.energy());
    BOOST_CHECK_CLOSE(energyTCPlus, energyTCMinus, 1.0E-8);
  }
}

#endif // HAVE_TwoU1 and DMRG_TC