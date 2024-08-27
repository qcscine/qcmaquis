/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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
#include "utils/fpcomparison.h"
#include "utils/io.hpp"
#include <iostream>
#include "maquis_dmrg.h"
#include "Fixtures/TimeEvolversFixture.h"

/**
 * @brief Tests that the energy is conserved along a relativistic TD-DMRG
 * propagation. The data are obtained for N2+ and the 3-21G basis set.
 */
BOOST_FIXTURE_TEST_CASE(TestImaginaryTimeRelativistic, TestTimeEvolverFixture) {
#if defined(HAVE_U1DG) and defined(DMRG_TD)
  // Two-site evolutions
  parametersRelativistic.set("optimization", "twosite");
  parametersRelativistic.set("time_step", 10.);
  parametersRelativistic.set("nsweeps", 100);
  // TD
  maquis::DMRGInterface<std::complex<double>> interfaceTD(parametersRelativistic
  );
  interfaceTD.evolve();
  auto energyTD = std::real(interfaceTD.energy());
  // TI
  parametersRelativistic.set("nsweeps", 40);
  maquis::DMRGInterface<std::complex<double>> interfaceTI(parametersRelativistic
  );
  interfaceTI.optimize();
  auto energyTI = std::real(interfaceTI.energy());
  // The threshold is here a bit looser because the iTD-DMRG convergence is
  // rather slow
  BOOST_CHECK_CLOSE(energyTD, energyTI, 1.0E-8);
#endif  // HAVE_U1DG and DMRG_TD
}