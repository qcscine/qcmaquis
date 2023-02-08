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

#define BOOST_TEST_MODULE SweepOptimizationTraits

#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/SweepBasedAlgorithms/SweepOptimizationTypeTrait.h"

/** @brief Test of traits for single-site optimizer */
BOOST_AUTO_TEST_CASE(Test_SweepOptimizationTrait_SS)
{
  using SweepTraitClass = SweepOptimizationTypeTrait<SweepOptimizationType::SingleSite>;
  int latticeSize = 10;
  auto numberOfMicroiterations = SweepTraitClass::getNumberOfMicroiterations(latticeSize);
  BOOST_CHECK_EQUAL(numberOfMicroiterations, 18);
  // Checks which site is visited, and how many times
  std::vector<int> visitedSites;
  std::vector<SweepDirectionType> directions;
  for (int iSite = 0; iSite < numberOfMicroiterations; iSite++) {
    visitedSites.push_back(SweepTraitClass::convertMicroIterationToSite(latticeSize, iSite));
    directions.push_back(SweepTraitClass::getSweepDirection(latticeSize, iSite));
  }
  // Checks the number of visits
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 0), 1);
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 3), 2);
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 5), 2);
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 9), 1);
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 10), 0);
  // Does the same for the directions
  BOOST_CHECK_EQUAL(std::count(directions.begin(), directions.end(), SweepDirectionType::Forward), 10);
  BOOST_CHECK_EQUAL(std::count(directions.begin(), directions.end(), SweepDirectionType::Backward), 8);
}

/** @brief Test of traits for two-site optimizer */
BOOST_AUTO_TEST_CASE(Test_SweepOptimizationTrait_TS)
{
  using SweepTraitClass = SweepOptimizationTypeTrait<SweepOptimizationType::TwoSite>;
  int latticeSize = 10;
  auto numberOfMicroiterations = SweepTraitClass::getNumberOfMicroiterations(latticeSize);
  BOOST_CHECK_EQUAL(numberOfMicroiterations, 16);
  // Checks which site is visited, and how many times
  std::vector<int> visitedSites;
  std::vector<SweepDirectionType> directions;
  for (int iSite = 0; iSite < numberOfMicroiterations; iSite++) {
    visitedSites.push_back(SweepTraitClass::convertMicroIterationToSite(latticeSize, iSite));
    directions.push_back(SweepTraitClass::getSweepDirection(latticeSize, iSite));
  }
  // Checks the number of visits
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 0), 1);
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 2), 2);
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 6), 2);
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 8), 1);
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 9), 0);
  BOOST_CHECK_EQUAL(std::count(visitedSites.begin(), visitedSites.end(), 10), 0);
  // Does the same for the directions
  BOOST_CHECK_EQUAL(std::count(directions.begin(), directions.end(), SweepDirectionType::Forward), 9);
  BOOST_CHECK_EQUAL(std::count(directions.begin(), directions.end(), SweepDirectionType::Backward), 7);
}