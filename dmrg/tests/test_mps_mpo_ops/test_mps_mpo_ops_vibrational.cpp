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

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_rotate.h"
#include "Fixtures/NModeFixture.h"
#include "Fixtures/WatsonFixture.h"
#include "dmrg/sim/matrix_types.h"

#ifdef DMRG_VIBRATIONAL

BOOST_FIXTURE_TEST_CASE(Test_ExpVal_None_HarmonicEnergy, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
    parametersEthyleneWatsonHarmonic.set("init_state", "basis_state_generic");
    parametersEthyleneWatsonHarmonic.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
    auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatsonHarmonic);
    auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
    auto mps = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
    BOOST_CHECK_CLOSE(norm(mps), 1., 1.0E-15);
    auto energy = expval(mps, watsonHarmonicMPO)/norm(mps);
    BOOST_CHECK_CLOSE(referenceHarmonicEnergy, energy, 1e-7);
#endif // HAVE_TrivialGroup
}

#ifdef HAVE_TrivialGroup

BOOST_FIXTURE_TEST_CASE(Test_ExpVal_None_BraKetHermitian, WatsonFixture)
{
    parametersEthyleneWatsonHarmonic.set("init_state", "const");
    auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatsonHarmonic);
    auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
    auto mpsConst = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
    //
    parametersEthyleneWatsonHarmonic.set("init_state", "default");
    lattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto mpsDefault = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
    //
    auto energy1 = expval(mpsConst, mpsDefault, watsonHarmonicMPO);
    auto energy2 = expval(mpsDefault, mpsConst, watsonHarmonicMPO);
    BOOST_CHECK_CLOSE(energy1, energy2, 1e-7);
}

#endif // HAVE_TrivialGroup

BOOST_FIXTURE_TEST_CASE(Test_ExpVal_NU1_SameBraKet, NModeFixture)
{
#ifdef HAVE_NU1
    parametersFADTwoBody.set("init_state", "default");
    auto lattice = Lattice(parametersFADTwoBody);
    auto nModeModel = Model<matrix, NU1_template<2>>(lattice, parametersFADTwoBody);
    auto nModeMPO = make_mpo(lattice, nModeModel);
    auto mps = MPS<matrix, NU1_template<2>>(lattice.size(), *(nModeModel.initializer(lattice, parametersFADTwoBody)));
    auto energyBeforeNormalization = expval(mps, nModeMPO)/norm(mps);
    // Now normalizes the MPS
    auto mpsNorm = std::sqrt(norm(mps));
    mps[0] /= mpsNorm;
    auto energyAfterNormalization = expval(mps, nModeMPO);
    BOOST_CHECK_CLOSE(energyBeforeNormalization, energyAfterNormalization, 1e-7);
#endif // HAVE_NU1
}

#endif // DMRG_VIBRATIONAL