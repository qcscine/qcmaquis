/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#ifdef DMRG_VIBRATIONAL

#include <boost/test/included/unit_test.hpp>
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/MaleimideFixture.h"

/** @brief Watson-based calculation on maleimide - one-body potential */
BOOST_FIXTURE_TEST_CASE(Test_vDMRG_Calculation_Maleimide, MaleimideFixture)
{
#if defined(HAVE_TrivialGroup) && ORDER_NONE >= 12
  using InterfaceType = maquis::DMRGInterface<double>;
  using MPSType = MPS<matrix, TrivialGroup>;
  // Adds the final input parameters
  parametersMaleimideOneBody.set("init_type", "const");
  parametersMaleimideOneBody.set("nsweeps", 10);
  parametersMaleimideOneBody.set("ngrowsweeps", 2);
  parametersMaleimideOneBody.set("nmainsweeps", 2);
  parametersMaleimideOneBody.set("optimize", "twosite");
  parametersMaleimideOneBody.set("twosite_truncation", "heev_truncate");
  parametersMaleimideOneBody.set("alpha_initial", 1.0E-8);
  parametersMaleimideOneBody.set("alpha_main", 1.0E-15);
  parametersMaleimideOneBody.set("alpha_final", 0.);
  parametersMaleimideOneBody.set("max_bond_dimension", 20);
  // Creates the interface and calculates the corresponding energy
  InterfaceType interface(parametersMaleimideOneBody);
  interface.optimize();
  auto energyFromInterface = maquis::real(interface.energy());
  // Now calculates the energy of the HF determinant
  parametersMaleimideOneBody.set("init_type", "basis_state_generic");
  parametersMaleimideOneBody.set("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
  auto maleimideLattice = Lattice(parametersMaleimideOneBody);
  auto maleimideModel = Model<matrix, TrivialGroup>(maleimideLattice, parametersMaleimideOneBody);
  auto mpsHF = MPSType(maleimideLattice.size(), *(maleimideModel.initializer(maleimideLattice, parametersMaleimideOneBody)));
  auto mpoMaleimide = make_mpo(maleimideLattice, maleimideModel);
  auto energyFromHF = maquis::real(expval(mpsHF, mpoMaleimide)/norm(mpsHF));
  BOOST_CHECK(energyFromInterface < energyFromHF);
#endif // HAVE_TrivialGroup && ORDER_NONE >= 12
}

#endif // DMRG_VIBRATIONAL