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


#define BOOST_TEST_MODULE TranscorrelatedHelium

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
BOOST_FIXTURE_TEST_CASE(TestTCMolecular_HE_NormalOrdered, TranscorrelatedFixture)
{
    parametersHeTranscorrelated.set("truncation_initial", 1e-50);
    parametersHeTranscorrelated.set("truncation_final", 1e-50);
    parametersHeTranscorrelated.set("init_state", "default");
    parametersHeTranscorrelated.set("site_types", "0,0,0,0,0");
    parametersHeTranscorrelated.set("max_bond_dimension", 100);
    parametersHeTranscorrelated.set("propagator_accuracy", 1.0E-10);
    parametersHeTranscorrelated.set("propagator_maxiter", 10);
    parametersHeTranscorrelated.set("nsweeps", 0);
    parametersHeTranscorrelated.set("hamiltonian_units", "Hartree");
    parametersHeTranscorrelated.set("time_units", "fs");
    parametersHeTranscorrelated.set("TD_backpropagation", "no");
    parametersHeTranscorrelated.set("symmetry", "2u1");
    parametersHeTranscorrelated.set("time_step", 0.1);
    parametersHeTranscorrelated.set("optimization", "singlesite");
    parametersHeTranscorrelated.set("transcorrelated_nsweeps_TI", 0);
    parametersHeTranscorrelated.set("transcorrelated_nsweeps_TC", 20);
    parametersHeTranscorrelated.set("integral_file", CONVENTIONAL_HE_FCIDUMP_PATH);
    parametersHeTranscorrelated.set("transcorrelated_integral_file", TRANSCORRELATED_HE_FCIDUMP_PATH);

    maquis::DMRGInterface<double> interface(parametersHeTranscorrelated);
    interface.runTranscorrelated();

    auto energyDMRG = maquis::real(interface.energy());

    parametersHeTranscorrelated.set("transcorrelated_3body_max_coupling", 4);
    maquis::DMRGInterface<double> NO_interface(parametersHeTranscorrelated);
    NO_interface.runTranscorrelated();

    auto energyDMRG_NO = maquis::real(interface.energy());

    BOOST_CHECK_SMALL(std::abs(energyDMRG-energyDMRG_NO), 1.0E-10);

    maquis::cout << energyDMRG << " " << energyDMRG_NO << std::endl;
}