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

#define BOOST_TEST_MAIN

#ifdef DMRG_VIBRATIONAL

#include <boost/test/included/unit_test.hpp>
#include "Fixtures/WatsonFixture.h"
#include "maquis_dmrg.h"

/** @brief Checks that increasing NMax leads to a lower energy for H2CO */
BOOST_FIXTURE_TEST_CASE(Test_DMRG_, WatsonFixture)
{
#ifdef HAVE_TrivialGroup
  using InterfaceType = maquis::DMRGInterface<double>;
  // Adds the final input parameters
  parametersH2COWatsonNoCoriolis.set("init_type", "basis_state_generic");
  parametersH2COWatsonNoCoriolis.set("init_basis_state", "0,0,0,0,0,0");
  parametersH2COWatsonNoCoriolis.set("optimization", "singlesite");
  parametersH2COWatsonNoCoriolis.set("alpha_initial", 1.0E-8);
  parametersH2COWatsonNoCoriolis.set("alpha_initial", 1.0E-15);
  parametersH2COWatsonNoCoriolis.set("alpha_initial", 0.);
  parametersH2COWatsonNoCoriolis.set("nsweeps", 20);
  parametersH2COWatsonNoCoriolis.set("ngrowsweeps", 2);
  parametersH2COWatsonNoCoriolis.set("nmainsweeps", 2);
  parametersH2COWatsonNoCoriolis.set("max_bond_dimension", 50);
  parametersH2COWatsonNoCoriolis.set("MODEL", "watson");
  parametersH2COWatsonNoCoriolis.set("Nmax", "3,2,3,2,2,4");
  InterfaceType interface(parametersH2COWatsonNoCoriolis);
  interface.optimize();
  parametersH2COWatsonNoCoriolis.set("Nmax", 2);
  InterfaceType interfaceSmaller(parametersH2COWatsonNoCoriolis);
  interfaceSmaller.optimize();
  BOOST_TEST(interface.energy() < interfaceSmaller.energy());
#endif // HAVE_TrivialGroup
}

#endif // DMRG_VIBRATIONAL