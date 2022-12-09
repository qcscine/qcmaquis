/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
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

#define BOOST_TEST_MODULE mps_vibrational

#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "Fixtures/NModeFixture.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/mps_initializers_helper.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/generate_mpo.hpp"

BOOST_FIXTURE_TEST_CASE(Test_MPS_Vibrational_NU1, NModeFixture)
{
#ifdef HAVE_NU1
  using Symmetry = NU1_template<2>;
  using MPSType = MPS<matrix, Symmetry>;
  // Populates the physical indices
  parametersFADTwoBody.set("init_type", "const");
  auto lattice = Lattice(parametersFADTwoBody);
  int latticeSize = lattice.size();
  auto nModeModel = Model<matrix, Symmetry>(lattice, parametersFADTwoBody);
  auto mps = MPSType(lattice.size(), *(nModeModel.initializer(lattice, parametersFADTwoBody)));
  double normBefore = norm(mps);
  mps.scaleByScalar(3.);
  double normAfter = norm(mps);
  BOOST_CHECK_CLOSE(std::sqrt(normAfter/normBefore), 3., 5.0E-14);
#endif
}