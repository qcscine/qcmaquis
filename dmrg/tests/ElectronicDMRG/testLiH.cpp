/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022 by Alberto Baiardi <abaiardi@ethz.ch>
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

#define BOOST_TEST_MODULE TestLiH

#include <boost/mpl/list.hpp>
#include <boost/test/included/unit_test.hpp>
#include "Fixtures/LiHFixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/** @brief Test FEAST for electronic calculations (we target the ground and first excited state) */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_LiH_DMRG_SSvsTS, S, symmetries, LiHFixture)
{
  // Generic parameters
  parametersLiH.set("max_bond_dimension", 50);
  parametersLiH.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
  parametersLiH.set("init_type", "const");
  parametersLiH.set("nsweeps", 20);
  parametersLiH.set("ngrowsweeps", 2);
  parametersLiH.set("nmainsweeps", 5);
  parametersLiH.set("optimizer", "singlesite");
  // Constructs the interface for the SS optimizer
  maquis::DMRGInterface<double> ssOptimizer(parametersLiH);
  ssOptimizer.optimize();
  auto energyFromSS = ssOptimizer.energy();
  // Now runs the TS optimizer
  parametersLiH.set("optimizer", "twosite");
  maquis::DMRGInterface<double> tsOptimizer(parametersLiH);
  tsOptimizer.optimize();
  auto energyFromTS = tsOptimizer.energy();
  BOOST_CHECK_CLOSE(energyFromSS, energyFromTS, 1.0E-8);
}

///** @brief Test conventional DMRG with dumping the boundaries to File */
//BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_LiH_DMRG_BoundaryStorage, S, symmetries, LiHFixture)
//{
//  // Generic parameters
//  parametersLiH.set("max_bond_dimension", 50);
//  parametersLiH.set("init_type", "const");
//  parametersLiH.set("seed", 42);
//  parametersLiH.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
//  parametersLiH.set("nsweeps", 20);
//  parametersLiH.set("optimization", "twosite");
//  parametersLiH.set("storagedir", "tmp");
//  maquis::DMRGInterface<double> optimizer(parametersLiH);
//  optimizer.optimize();
//  // auto energyFromTS = optimizer.energy();
//  // BOOST_CHECK_CLOSE(energyFromSS, energyFromTS, 1.0E-8);
//}