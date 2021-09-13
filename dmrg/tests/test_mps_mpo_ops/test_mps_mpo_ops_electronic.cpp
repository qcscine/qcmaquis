/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
*               2021-   Alberto Baiardi <abaiardi@ethz.ch>
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

// This files contains the tests that check the functionality of the 
// methods stored in the [mps_mpo_ops] files and that are generic, 
// meaning that they can be tested for both the TwoU1 and the SU2U1
// group.

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include <boost/mpl/list.hpp>
#include "utils/io.hpp" // has to be first include because of impi
#include <iostream>
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_rotate.h"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/BenzeneFixture.h"
#include "test_mps.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/**
 * @brief Checks that the MPO-MPS contraction fulfill the Hermitianity of the Hamiltonian 
 * Note that this test uses the "standard" implementation of the expval function, which uses the left boundary.
 */
BOOST_FIXTURE_TEST_CASE_TEMPLATE( Test_MPS_MPO_Hermitian_Left_Electronic, S, symmetries, BenzeneFixture )
{
    // Generates the HF MPS
    auto lattice = Lattice(parametersBenzene);
    auto modelHF = Model<matrix, S>(lattice, parametersBenzene);
    auto mpsHF = MPS<matrix, S>(lattice.size(), *(modelHF.initializer(lattice, parametersBenzene)));
    // Generates the const guess. Note that the bond dimension will be != than that of the HF 
    // guess, so we check that the contraction routines work also for MPS with different bond dimension.
    parametersBenzene.set("init_state", "const");
    auto modelConst = Model<matrix, S>(lattice, parametersBenzene);
    auto mpsConst = MPS<matrix, S>(lattice.size(), *(modelConst.initializer(lattice, parametersBenzene)));
    // Generates the default guess
    parametersBenzene.set("init_state", "default");
    auto modelDefault = Model<matrix, S>(lattice, parametersBenzene);
    auto mpsDefault = MPS<matrix, S>(lattice.size(), *(modelDefault.initializer(lattice, parametersBenzene)));
    // Calculates the MPO (any model would work here)
    auto mpo = make_mpo(lattice, modelConst);
    // Check 1
    double expVal1 = expval(mpsDefault, mpsConst, mpo);
    double expVal2 = expval(mpsConst, mpsDefault, mpo);
    BOOST_CHECK_CLOSE(expVal1, expVal2, 1.E-10);
    // Check 2
    expVal1 = expval(mpsHF, mpsConst, mpo);
    expVal2 = expval(mpsConst, mpsHF, mpo);
    BOOST_CHECK_CLOSE(expVal1, expVal2, 1.E-10);
    // Check 3
    expVal1 = expval(mpsDefault, mpsHF, mpo);
    expVal2 = expval(mpsHF, mpsDefault, mpo);
    BOOST_CHECK_CLOSE(expVal1, expVal2, 1.E-10);
}

/** @brief Same as above, but with the right boundaries */
BOOST_FIXTURE_TEST_CASE_TEMPLATE( Test_MPS_MPO_Hermitian_Right_Electronic, S, symmetries, BenzeneFixture )
{
    // Generates the HF MPS
    auto lattice = Lattice(parametersBenzene);
    auto modelHF = Model<matrix, S>(lattice, parametersBenzene);
    auto mpsHF = MPS<matrix, S>(lattice.size(), *(modelHF.initializer(lattice, parametersBenzene)));
    // Generates the const guess. Note that the bond dimension will be != than that of the HF 
    // guess, so we check that the contraction routines work also for MPS with different bond dimension.
    parametersBenzene.set("init_state", "const");
    auto modelConst = Model<matrix, S>(lattice, parametersBenzene);
    auto mpsConst = MPS<matrix, S>(lattice.size(), *(modelConst.initializer(lattice, parametersBenzene)));
    // Generates the default guess
    parametersBenzene.set("init_state", "default");
    auto modelDefault = Model<matrix, S>(lattice, parametersBenzene);
    auto mpsDefault = MPS<matrix, S>(lattice.size(), *(modelDefault.initializer(lattice, parametersBenzene)));
    // Calculates the MPO (any model would work here)
    auto mpo = make_mpo(lattice, modelConst);
    // Check 1
    double expVal1 = expvalFromRight(mpsDefault, mpsConst, mpo);
    double expVal2 = expvalFromRight(mpsConst, mpsDefault, mpo);
    BOOST_CHECK_CLOSE(expVal1, expVal2, 1.E-10);
    // Check 2
    expVal1 = expvalFromRight(mpsHF, mpsConst, mpo);
    expVal2 = expvalFromRight(mpsConst, mpsHF, mpo);
    BOOST_CHECK_CLOSE(expVal1, expVal2, 1.E-10);
    // Check 3
    expVal1 = expvalFromRight(mpsDefault, mpsHF, mpo);
    expVal2 = expvalFromRight(mpsHF, mpsDefault, mpo);
    BOOST_CHECK_CLOSE(expVal1, expVal2, 1.E-10);
}