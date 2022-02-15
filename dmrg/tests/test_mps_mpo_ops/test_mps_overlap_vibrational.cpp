/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
*               2021- Alberto Baiardi <abaiardi@ethz.ch>
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

#define BOOST_TEST_MODULE MPSOverlapVibrational

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/WatsonFixture.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

#ifdef DMRG_VIBRATIONAL


/** @brief Checks that overlap between two ONV is zero for the TrivialGroup case */
BOOST_FIXTURE_TEST_CASE( Test_MPS_Overlap_Vibrational, WatsonFixture )
{
#ifdef HAVE_TrivialGroup
    // ONV1 
    parametersEthyleneWatson.set("init_state", "basis_state_generic");
    parametersEthyleneWatson.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
    auto lattice = Lattice(parametersEthyleneWatson);
    auto model = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatson);
    auto mpsHF1 = MPS<matrix, TrivialGroup>(lattice.size(), *(model.initializer(lattice, parametersEthyleneWatson)));
    // ONV2
    parametersEthyleneWatson.set("init_state", "basis_state_generic");
    parametersEthyleneWatson.set("init_basis_state", "0,0,0,0,0,0,0,0,1,0,0,0");
    lattice = Lattice(parametersEthyleneWatson);
    model = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatson);
    auto mpsHF2 = MPS<matrix, TrivialGroup>(lattice.size(), *(model.initializer(lattice, parametersEthyleneWatson)));
    // Overlap calculationb
    double overlapHF1 = overlap(mpsHF1, mpsHF2);
    double overlapHF2 = overlap(mpsHF2, mpsHF1);
    BOOST_CHECK_CLOSE(overlapHF1, overlapHF2, 1.E-10);
    BOOST_CHECK_CLOSE(overlapHF1, 0., 1.E-10);
#endif // HAVE_TrivialGroup
}

#ifdef HAVE_TrivialGroup

/** @brief Checks Hermitianity of the overlap calculation for the real-valued TrivialGroup case */
BOOST_FIXTURE_TEST_CASE( Test_MPS_Overlap_Hermitian_Vibrational, WatsonFixture )
{
    // Global variables
    auto lattice = Lattice(parametersEthyleneWatson);
    auto model = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatson);
    // Modifies the init parameters to create two different MPSs
    parametersEthyleneWatson.set("init_state", "const");
    auto mpsConst = MPS<matrix, TrivialGroup>(lattice.size(), *(model.initializer(lattice, parametersEthyleneWatson)));
    parametersEthyleneWatson.set("init_state", "default");
    auto mpsDefault = MPS<matrix, TrivialGroup>(lattice.size(), *(model.initializer(lattice, parametersEthyleneWatson)));
    // Calculates the overlap
    double overlapOriginal = overlap(mpsConst, mpsDefault);
    double overlapHerm = overlap(mpsDefault, mpsConst);
    BOOST_CHECK_CLOSE(overlapOriginal, overlapHerm, 1.E-10);
}

#endif // HAVE_TrivialGroup

#endif // DMRG_VIBRATIONAL
