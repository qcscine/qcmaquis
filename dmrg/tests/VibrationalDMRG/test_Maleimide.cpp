/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
 *               2022- by Nina Glaser <nglaser@ethz.ch>
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
  parametersMaleimideOneBody.set("init_state", "const");
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
  parametersMaleimideOneBody.set("init_state", "basis_state_generic");
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