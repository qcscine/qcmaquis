/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#ifdef DMRG_VIBRATIONAL

#include <boost/test/included/unit_test.hpp>
#include "Fixtures/WatsonFixture.h"
#include "maquis_dmrg.h"

/** @brief Checks that increasing NMax leads to a lower energy for H2CO */
BOOST_FIXTURE_TEST_CASE(Test_DMRG_H2CO, WatsonFixture)
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
  // parametersH2COWatsonNoCoriolis.set("Nmax", "3,2,3,2,2,4");
  parametersH2COWatsonNoCoriolis.set("Nmax", "3");
  InterfaceType interface(parametersH2COWatsonNoCoriolis);
  interface.optimize();
  parametersH2COWatsonNoCoriolis.set("Nmax", 2);
  InterfaceType interfaceSmaller(parametersH2COWatsonNoCoriolis);
  interfaceSmaller.optimize();
  BOOST_TEST(interface.energy() < interfaceSmaller.energy());
#endif // HAVE_TrivialGroup
}

#endif // DMRG_VIBRATIONAL