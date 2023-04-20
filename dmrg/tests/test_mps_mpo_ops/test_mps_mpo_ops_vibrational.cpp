/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_rotate.h"
#include "Fixtures/WatsonFixture.h"
#include "dmrg/sim/matrix_types.h"

#ifdef DMRG_VIBRATIONAL

BOOST_FIXTURE_TEST_CASE(Test_ExpVal_None_HarmonicEnergy, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
    parametersEthyleneWatsonHarmonic.set("init_type", "basis_state_generic");
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
    parametersEthyleneWatsonHarmonic.set("init_type", "const");
    auto lattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersEthyleneWatsonHarmonic);
    auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
    auto mpsConst = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
    //
    parametersEthyleneWatsonHarmonic.set("init_type", "default");
    lattice = Lattice(parametersEthyleneWatsonHarmonic);
    auto mpsDefault = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersEthyleneWatsonHarmonic)));
    //
    auto energy1 = expval(mpsConst, mpsDefault, watsonHarmonicMPO);
    auto energy2 = expval(mpsDefault, mpsConst, watsonHarmonicMPO);
    BOOST_CHECK_CLOSE(energy1, energy2, 1e-7);
}

BOOST_FIXTURE_TEST_CASE(Test_ExpVal_H2CO_BraKetHermitian, WatsonFixture)
{
    parametersH2COWatsonInternal.set("init_type", "const");
    auto lattice = Lattice(parametersH2COWatsonInternal);
    auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersH2COWatsonInternal);
    auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
    auto mpsConst = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersH2COWatsonInternal)));
    //
    parametersH2COWatsonInternal.set("init_type", "default");
    lattice = Lattice(parametersH2COWatsonInternal);
    auto mpsDefault = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersH2COWatsonInternal)));
    //
    auto energy1 = expval(mpsConst, mpsDefault, watsonHarmonicMPO);
    auto energy2 = expval(mpsDefault, mpsConst, watsonHarmonicMPO);
    BOOST_CHECK_CLOSE(energy1, energy2, 1.0e-7);
}

#endif // HAVE_TrivialGroup

#endif // DMRG_VIBRATIONAL
