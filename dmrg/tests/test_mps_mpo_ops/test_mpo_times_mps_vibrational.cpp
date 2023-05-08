/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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

#ifdef DMRG_VIBRATIONAL

/**
 * @brief Checks that [mpo_times_mps] gives results that are coherent with expval.
 * Note that here we use 1) the harmonic MPS and 2) the harmonic Hamiltonian.
 */
BOOST_FIXTURE_TEST_CASE(Test_MPO_Times_MPS_None, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
  // Generates the HF MPS
  parametersEthyleneWatsonHarmonic.set("init_type", "basis_state_generic");
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
  parametersEthyleneWatsonHarmonic.set("init_type", "const");
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

/** @brief Check Hermiticity of None operator with Coriolis */
BOOST_FIXTURE_TEST_CASE(Test_MPO_Times_MPS_Hermiticity_None, WatsonFixture)
{
  // Generates the HF MPS
  parametersH2COWatson.set("init_type", "const");
  parametersH2COWatson.set("max_bond_dimension", 200);
  auto lattice = Lattice(parametersH2COWatson);
  auto watsonModel = Model<matrix, TrivialGroup>(lattice, parametersH2COWatson);
  auto watsonHarmonicMPO = make_mpo(lattice, watsonModel);
  auto mps1 = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersH2COWatson)));
  // Calculates the default guess
  parametersH2COWatson.set("init_type", "default");
  auto mps2 = MPS<matrix, TrivialGroup>(lattice.size(), *(watsonModel.initializer(lattice, parametersH2COWatson)));
  auto energy1 = expval(mps1, mps2, watsonHarmonicMPO);
  auto energy2 = expval(mps2, mps1, watsonHarmonicMPO);
  BOOST_CHECK_CLOSE(energy1, energy2, 1.0E-11);
};

/**
 * @brief Checks that the H^2 expectation value, caluclated via [mpo_times_mps] and expva, gives coherent results.
 */
BOOST_FIXTURE_TEST_CASE(Test_MPO_Times_MPS_Variance_None, WatsonFixture)
{
  // Generates the HF MPS
  parametersEthyleneWatsonHarmonic.set("init_type", "basis_state_generic");
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

#endif // DMRG_VIBRATIONAL
