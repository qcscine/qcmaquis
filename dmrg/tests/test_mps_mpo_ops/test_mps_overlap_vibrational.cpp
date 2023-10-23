/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
    parametersEthyleneWatson.set("init_type", "basis_state_generic");
    parametersEthyleneWatson.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0");
    auto lattice = Lattice(parametersEthyleneWatson);
    auto model = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatson);
    auto mpsHF1 = MPS<matrix, TrivialGroup>(lattice.size(), *(model.initializer(lattice, parametersEthyleneWatson)));
    // ONV2
    parametersEthyleneWatson.set("init_type", "basis_state_generic");
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
    parametersEthyleneWatson.set("init_type", "const");
    auto mpsConst = MPS<matrix, TrivialGroup>(lattice.size(), *(model.initializer(lattice, parametersEthyleneWatson)));
    parametersEthyleneWatson.set("init_type", "default");
    auto mpsDefault = MPS<matrix, TrivialGroup>(lattice.size(), *(model.initializer(lattice, parametersEthyleneWatson)));
    // Calculates the overlap
    double overlapOriginal = overlap(mpsConst, mpsDefault);
    double overlapHerm = overlap(mpsDefault, mpsConst);
    BOOST_CHECK_CLOSE(overlapOriginal, overlapHerm, 1.E-10);
}

#endif // HAVE_TrivialGroup

#endif // DMRG_VIBRATIONAL
