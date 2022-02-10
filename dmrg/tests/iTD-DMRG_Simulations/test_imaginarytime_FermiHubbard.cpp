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

#include "Fixtures/TranscorrelatedFixture.h"
#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include "utils/io.hpp"
#include <iostream>
#include "maquis_dmrg.h"

/** @brief iTD-DMRG calculation on the Fermi-Hubbard model */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimevsTIFermiHubbard_2x2, TranscorrelatedFixture)
{
#ifdef HAVE_TwoU1
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("nsweeps", 10);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("max_bond_dimension", 300);
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

/** @brief iTD-DMRG calculation on the 2x2 Fermi-Hubbard model */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_2x2_Transcorrelated_J0, TranscorrelatedFixture)
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

/** @brief iTD-DMRG calculation on the 2x2 Fermi-Hubbard model for J=+/-0.1 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_2x2_Transcorrelated_J0p1, TranscorrelatedFixture)
{
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("nsweeps", 100);
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
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("J_Transcorrelated", 0.1);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("hamiltonian_units", "Hartree");
    maquis::DMRGInterface<std::complex<double>> interfaceTCPlus(parameters2x2_RealSpace_U4_2Alpha1Beta);
    interfaceTCPlus.evolve();
    auto energyTCPlus = std::real(interfaceTCPlus.energy());
    //
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("J_Transcorrelated", -0.1);
    maquis::DMRGInterface<std::complex<double>> interfaceTCMinus(parameters2x2_RealSpace_U4_2Alpha1Beta);
    interfaceTCMinus.evolve();
    auto energyTCMinus = std::real(interfaceTCMinus.energy());
    BOOST_CHECK_CLOSE(energyTCPlus, energyTCMinus, 1.0E-8);
}

/** @brief iTD-DMRG calculation on the 2x2 Fermi-Hubbard model for J=+/-1 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_2x2_Transcorrelated_J1p0, TranscorrelatedFixture)
{
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("nsweeps", 100);
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
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("J_Transcorrelated", 1.0);
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("hamiltonian_units", "Hartree");
    maquis::DMRGInterface<std::complex<double>> interfaceTCPlus(parameters2x2_RealSpace_U4_2Alpha1Beta);
    interfaceTCPlus.evolve();
    auto energyTCPlus = std::real(interfaceTCPlus.energy());
    //
    parameters2x2_RealSpace_U4_2Alpha1Beta.set("J_Transcorrelated", -1.0);
    maquis::DMRGInterface<std::complex<double>> interfaceTCMinus(parameters2x2_RealSpace_U4_2Alpha1Beta);
    interfaceTCMinus.evolve();
    auto energyTCMinus = std::real(interfaceTCMinus.energy());
    BOOST_CHECK_CLOSE(energyTCPlus, energyTCMinus, 1.0E-8);
}

/** @brief iTD-DMRG calculation on the Fermi-Hubbard model */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimevsTIFermiHubbard_2x2_Momentum, TranscorrelatedFixture)
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

/** @brief iTD-DMRG calculation on the 2x2 Fermi-Hubbard model for J=+/-1 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeFermiHubbard_2x2_Momentum_Transcorrelated_J1p0, TranscorrelatedFixture)
{
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("nsweeps", 100);
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
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("J_Transcorrelated", 1.0);
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("hamiltonian_units", "Hartree");
    maquis::DMRGInterface<std::complex<double>> interfaceTCPlus(parameters2x2_MomentumSpace_U4_2Alpha1Beta);
    interfaceTCPlus.evolve();
    auto energyTCPlus = std::real(interfaceTCPlus.energy());
    //
    parameters2x2_MomentumSpace_U4_2Alpha1Beta.set("J_Transcorrelated", -1.0);
    maquis::DMRGInterface<std::complex<double>> interfaceTCMinus(parameters2x2_MomentumSpace_U4_2Alpha1Beta);
    interfaceTCMinus.evolve();
    auto energyTCMinus = std::real(interfaceTCMinus.energy());
    BOOST_CHECK_CLOSE(energyTCPlus, energyTCMinus, 1.0E-8);
}

#endif // HAVE_TwoU1