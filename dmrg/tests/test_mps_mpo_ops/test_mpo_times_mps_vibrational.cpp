/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
*               2021 Alberto Baiardi <abaiardi@ethz.ch>
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

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/WatsonFixture.h"
#include "Fixtures/NModeFixture.h"

#ifdef DMRG_VIBRATIONAL

/**
 * @brief Checks that [mpo_times_mps] gives results that are coherent with expval.
 * Note that here we use 1) the harmonic MPS and 2) the harmonic Hamiltonian.
 */
BOOST_FIXTURE_TEST_CASE(Test_MPO_Times_MPS_None, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
  // Generates the HF MPS
  parametersEthyleneWatsonHarmonic.set("init_state", "basis_state_generic");
  parametersEthyleneWatsonHarmonic.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
  parametersEthyleneWatsonHarmonic.set("max_bond_dimension", 200);
  parametersEthyleneWatsonHarmonic.set("Nmax", 6);
  auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
  auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatsonHarmonic);
  auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
  auto mps = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
  // Calculates the MPS-MPO contraction
  auto traitClass = MPOTimesMPSTraitClass<matrix, TrivialGroup>(mps, watsonModel, lattice, watsonModel.total_quantum_numbers(parametersEthyleneWatsonHarmonic),
                                                                parametersEthyleneWatsonHarmonic["max_bond_dimension"]);
  auto outputMPS = traitClass.applyMPO(watsonHarmonicMPO);
  // Calculates the energy in two ways and check that the results are coherent
  auto energyFromMPSTimesMPO = overlap(mps, outputMPS)/norm(mps) + watsonHarmonicMPO.getCoreEnergy();
  auto energyFromExpVal = expval(mps, watsonHarmonicMPO)/norm(mps);
  BOOST_CHECK_CLOSE(energyFromMPSTimesMPO, energyFromExpVal, 1.E-10);
#endif
};

#ifdef HAVE_TrivialGroup

/** @brief Same as above, but with the const guess (and, therefore, m>1) */
BOOST_FIXTURE_TEST_CASE(Test_MPO_Times_MPS_Const_None, WatsonFixture)
{
  // Generates the HF MPS
  parametersEthyleneWatsonHarmonic.set("init_state", "const");
  parametersEthyleneWatsonHarmonic.set("max_bond_dimension", 200);
  parametersEthyleneWatsonHarmonic.set("Nmax", 4);
  auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
  auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatsonHarmonic);
  auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
  auto mps = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
  // Calculates the MPS-MPO contraction
  auto traitClass = MPOTimesMPSTraitClass<matrix, TrivialGroup>(mps, watsonModel, lattice, watsonModel.total_quantum_numbers(parametersEthyleneWatsonHarmonic),
                                                                parametersEthyleneWatsonHarmonic["max_bond_dimension"]);
  auto outputMPS = traitClass.applyMPO(watsonHarmonicMPO);
  // Calculates the energy in two ways and check that the results are coherent
  auto energyFromMPSTimesMPO = overlap(mps, outputMPS)/norm(mps) + watsonHarmonicMPO.getCoreEnergy();
  auto energyFromExpVal = expval(mps, watsonHarmonicMPO)/norm(mps);
  BOOST_CHECK_CLOSE(energyFromMPSTimesMPO, energyFromExpVal, 1.E-10);
};

/**
 * @brief Checks that the H^2 expectation value, caluclated via [mpo_times_mps] and expva, gives coherent results.
 */
BOOST_FIXTURE_TEST_CASE(Test_MPO_Times_MPS_Variance_None, WatsonFixture)
{
  // Generates the HF MPS
  parametersEthyleneWatsonHarmonic.set("init_state", "basis_state_generic");
  parametersEthyleneWatsonHarmonic.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
  parametersEthyleneWatsonHarmonic.set("max_bond_dimension", 10000);
  parametersEthyleneWatsonHarmonic.set("Nmax", 2);
  auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
  auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatsonHarmonic);
  auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
  auto mps = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
  // Calculates the MPS-MPO contraction
  auto traitClass = MPOTimesMPSTraitClass<matrix, TrivialGroup>(mps, watsonModel, lattice, watsonModel.total_quantum_numbers(parametersEthyleneWatsonHarmonic),
                                                                parametersEthyleneWatsonHarmonic["max_bond_dimension"]);
  auto outputMPS = traitClass.applyMPO(watsonHarmonicMPO);
  // Calculates the energy in two ways and check that the results are coherent
  auto squaredEnergy = overlap(outputMPS, outputMPS)/norm(mps) + watsonHarmonicMPO.getCoreEnergy();
  auto energyFromExpVal = expval(mps, watsonHarmonicMPO)/norm(mps);
  BOOST_CHECK_CLOSE(std::sqrt(squaredEnergy), energyFromExpVal, 1.E-10);
};

#endif

/**
 * @brief Checks that [mpo_times_mps] gives results that are coherent with expval.
 * Note that here we use the NU1 symmetry group and the default guess.
 */
BOOST_FIXTURE_TEST_CASE(Test_MPO_Times_MPS_NU1, NModeFixture)
{
#ifdef HAVE_NU1
  parametersFADTwoBody.set("init_state", "default");
  parametersFADTwoBody.set("max_bond_dimension", 10);
  auto lattice = Lattice(parametersFADTwoBody);
  auto nModeModel = Model<matrix, NU1_template<2>>(lattice, parametersFADTwoBody);
  auto nModeMPO = make_mpo(lattice, nModeModel);
  auto mps = MPS<matrix, NU1_template<2>>(lattice.size(), *(nModeModel.initializer(lattice, parametersFADTwoBody)));
  // Calculates the MPS-MPO contraction
  auto traitClass = MPOTimesMPSTraitClass<matrix, NU1_template<2>>(mps, nModeModel, lattice, nModeModel.total_quantum_numbers(parametersFADTwoBody),
                                                                   parametersFADTwoBody["max_bond_dimension"]);
  auto outputMPS = traitClass.applyMPO(nModeMPO);
  // Calculates the energy in two ways and check that the results are coherent
  auto energyFromMPSTimesMPO = overlap(mps, outputMPS)/norm(mps) + nModeMPO.getCoreEnergy();
  auto energyFromExpVal = expval(mps, nModeMPO)/norm(mps);
  BOOST_CHECK_CLOSE(energyFromMPSTimesMPO, energyFromExpVal, 1.E-10);
#endif
};

#endif // DMRG_VIBRATIONAL
