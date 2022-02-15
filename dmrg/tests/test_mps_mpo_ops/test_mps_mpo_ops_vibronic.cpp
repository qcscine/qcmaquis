/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
*               2022 Alberto Baiardi <abaiardi@ethz.ch>
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
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_rotate.h"
#include "Fixtures/VibronicFixture.h"
#include "dmrg/sim/matrix_types.h"

#ifdef DMRG_VIBRONIC

/** @brief Checks that the Hamonic energy is calculated correctly via [expval] */
BOOST_FIXTURE_TEST_CASE(Test_ExpVal_U1_HarmonicEnergy_Vibronic, VibronicFixture)
{
#ifdef HAVE_U1
    parametersVibronicPyrazineRedDim.set("init_state", "basis_state_generic");
    parametersVibronicPyrazineRedDim.set("init_basis_state", "1,0,0,0,0,0");
    auto vibronicLattice = Lattice(parametersVibronicPyrazineRedDim);
    auto vibronicModel = Model<matrix, U1>(vibronicLattice, parametersVibronicPyrazineRedDim);
    auto vibronicHarmonicMPO = make_mpo(vibronicLattice, vibronicModel);
    auto mps = MPS<matrix, U1>(vibronicLattice.size(), *(vibronicModel.initializer(vibronicLattice, parametersVibronicPyrazineRedDim)));
    BOOST_CHECK_CLOSE(norm(mps), 1., 1.0E-15);
    auto energy = expval(mps, vibronicHarmonicMPO)/norm(mps);
    auto referenceHarmonicEnergy = -2036.1425;
    BOOST_CHECK_CLOSE(referenceHarmonicEnergy, energy, 1e-7);
#endif // HAVE_U1
}

#ifdef HAVE_U1

/** @brief Checks that the Hamonic excitation energy is calculated correctly via [expval] */
BOOST_FIXTURE_TEST_CASE(Test_ExpVal_U1_HarmonicEnergyDifference_Vibronic, VibronicFixture)
{
    // Creates the vibrational ground state
    parametersVibronicPyrazineRedDim.set("init_state", "basis_state_generic");
    parametersVibronicPyrazineRedDim.set("init_basis_state", "1,0,0,0,0,0");
    auto vibronicLattice = Lattice(parametersVibronicPyrazineRedDim);
    auto vibronicModel = Model<matrix, U1>(vibronicLattice, parametersVibronicPyrazineRedDim);
    auto vibronicHarmonicMPO = make_mpo(vibronicLattice, vibronicModel);
    auto mpsGroundState = MPS<matrix, U1>(vibronicLattice.size(), *(vibronicModel.initializer(vibronicLattice, parametersVibronicPyrazineRedDim)));
    // Creates the vibrationally excited state
    parametersVibronicPyrazineRedDim.set("init_state", "basis_state_generic");
    parametersVibronicPyrazineRedDim.set("init_basis_state", "1,0,1,1,0,0");
    auto mpsExcitedState = MPS<matrix, U1>(vibronicLattice.size(), *(vibronicModel.initializer(vibronicLattice, parametersVibronicPyrazineRedDim)));
    // Calculates energy difference
    auto groundStateEnergy = expval(mpsGroundState, vibronicHarmonicMPO)/norm(mpsGroundState);
    auto excitedStateEnergy = expval(mpsExcitedState, vibronicHarmonicMPO)/norm(mpsExcitedState);
    // These numeric values are taken from the [VibronicFixture] file
    auto refEnergyDifference = (325.84799 + 559.34550)*2;
    BOOST_CHECK_CLOSE(refEnergyDifference, excitedStateEnergy-groundStateEnergy, 1e-7);
}

/** @brief Checks that the Hamonic excitation energy is calculated correctly via [expval] 
 *  (same as above, but for the excitonic Hamiltonian) */
BOOST_FIXTURE_TEST_CASE(Test_ExpVal_U1_HarmonicEnergyDifference_Excitonic, VibronicFixture)
{
    // General setup
    parametersExcitonicAggregate.set("vibronic_sorting", "intertwined");
    parametersExcitonicAggregate.set("init_state", "basis_state_generic");
    parametersExcitonicAggregate.set("init_basis_state", "1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
    auto excitonicLattice = Lattice(parametersExcitonicAggregate);
    auto excitonicModel = Model<matrix, U1>(excitonicLattice, parametersExcitonicAggregate);
    auto excitonicHarmonicMPO = make_mpo(excitonicLattice, excitonicModel);
    // Creates the vibrational ground state
    auto mpsGroundState = MPS<matrix, U1>(excitonicLattice.size(), *(excitonicModel.initializer(excitonicLattice, parametersExcitonicAggregate)));
    // Creates the vibrationally excited state
    parametersExcitonicAggregate.set("init_basis_state", "1,1,1,2,3,4,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
    auto mpsExcitedState = MPS<matrix, U1>(excitonicLattice.size(), *(excitonicModel.initializer(excitonicLattice, parametersExcitonicAggregate)));
    // Calculates energy difference
    auto groundStateEnergy = expval(mpsGroundState, excitonicHarmonicMPO)/norm(mpsGroundState);
    auto excitedStateEnergy = expval(mpsExcitedState, excitonicHarmonicMPO)/norm(mpsExcitedState);
    // These numeric values are taken from the [VibronicFixture] file
    auto refEnergyDifference = (103.00 + 105.50 + 2*270.00 + 3*276.00 + 4*375.50 + 662.50 + 734.50 + 814.50)*2;
    BOOST_CHECK_CLOSE(refEnergyDifference, excitedStateEnergy-groundStateEnergy, 1e-7);
}

/** @brief Checks that < bra | H | ket > == < ket | H | bra > */
BOOST_FIXTURE_TEST_CASE(Test_ExpVal_U1_BraKetHermitian_Vibronic, VibronicFixture)
{
    parametersVibronicPyrazineRedDimFull.set("init_state", "const");
    auto vibronicLattice = Lattice(parametersVibronicPyrazineRedDimFull);
    auto vibronicModel = Model<matrix, U1>(vibronicLattice, parametersVibronicPyrazineRedDimFull);
    auto vibronicMPO = make_mpo(vibronicLattice, vibronicModel);
    auto mpsConst = MPS<matrix, U1>(vibronicLattice.size(), *(vibronicModel.initializer(vibronicLattice, parametersVibronicPyrazineRedDimFull)));
    //
    parametersVibronicPyrazineRedDimFull.set("init_state", "default");
    vibronicLattice = Lattice(parametersVibronicPyrazineRedDimFull);
    auto mpsDefault = MPS<matrix, U1>(vibronicLattice.size(), *(vibronicModel.initializer(vibronicLattice, parametersVibronicPyrazineRedDimFull)));
    //
    auto energy1 = expval(mpsConst, mpsDefault, vibronicMPO);
    auto energy2 = expval(mpsDefault, mpsConst, vibronicMPO);
    BOOST_CHECK_CLOSE(energy1, energy2, 1e-7);
}

#endif // HAVE_U1

#endif // DMRG_VIBRATIONAL