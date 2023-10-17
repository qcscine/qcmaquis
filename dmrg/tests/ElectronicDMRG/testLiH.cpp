/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE TestLiH

#include <boost/mpl/list.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/filesystem/operations.hpp>
#include "Fixtures/LiHFixture.h"

//TODO NOTE THAT THE 2U1 VERSION OF THESE TESTS GETS STUCK IN A LOCAL MINIMUM -- TO BE CHECKED

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/** @brief Test FEAST for electronic calculations, both with SS and TS optimization */
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

/** @brief Test conventional DMRG with dumping the boundaries to File */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_LiH_DMRG_BoundaryStorage, S, symmetries, LiHFixture)
{
  // Generic parameters
  parametersLiH.set("max_bond_dimension", 50);
  parametersLiH.set("init_type", "default");
  parametersLiH.set("seed", 98789);
  parametersLiH.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
  parametersLiH.set("nsweeps", 20);
  parametersLiH.set("ngrowsweeps", 2);
  parametersLiH.set("nmainsweeps", 5);
  parametersLiH.set("optimization", "singlesite");
  parametersLiH.set("alpha_initial", 1.0E-8);
  parametersLiH.set("alpha_main", 1.0E-15);
  parametersLiH.set("alpha_final", 0.);
  parametersLiH.set("storagedir", "tmpDMRGSS");
  maquis::DMRGInterface<double> optimizerSS(parametersLiH);
  optimizerSS.optimize();
  auto energyFromSS = optimizerSS.energy();
  // Runs the same with the two-site optimizer
  parametersLiH.set("optimization", "twosite");
  parametersLiH.set("storagedir", "tmpDMRGTS");
  maquis::DMRGInterface<double> optimizerTS(parametersLiH);
  optimizerTS.optimize();
  auto energyFromTS = optimizerTS.energy();
  BOOST_CHECK_CLOSE(energyFromSS, energyFromTS, 1.0E-8);
  // This reference is taken from test2.cpp
  BOOST_CHECK_CLOSE(energyFromSS, referenceEnergy, 1.0e-7);
  boost::filesystem::remove_all("tmpDMRGSS");
  boost::filesystem::remove_all("tmpDMRGTS");
}

/** @brief Test DMRG-IPI with dumping the boundaries to File */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_LiH_IPI_BoundaryStorage, S, symmetries, LiHFixture)
{
  // Generic parameters
  parametersLiH.set("max_bond_dimension", 50);
  parametersLiH.set("init_type", "hf");
  parametersLiH.set("hf_occ", "4,1,1,1");
  parametersLiH.set("orbital_order", "4,1,2,3");
  parametersLiH.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
  parametersLiH.set("nsweeps", 20);
  parametersLiH.set("ngrowsweeps", 2);
  parametersLiH.set("nmainsweeps", 5);
  parametersLiH.set("optimization", "twosite");
  // IPI-specific parametrs
  parametersLiH.set("ipi_sweep_energy_threshold", 1.0E-10);
  parametersLiH.set("ipi_sweep_overlap_threshold", 1.0E-10);
  parametersLiH.set("ipi_sweeps_per_system", 5);
  parametersLiH.set("ipi_iterations", 5);
  parametersLiH.set("ipi_shift", -7.905);
  // DMRG-IPI calculation via interface without storing
  maquis::DMRGInterface<double> interface(parametersLiH);
  interface.runInversePowerIteration();
  auto energy = interface.energy();
  // DMRG-IPI calculation via interface with storing
  parametersLiH.set("storagedir", "tmpIPI");
  maquis::DMRGInterface<double> interfaceStorage(parametersLiH);
  interfaceStorage.runInversePowerIteration();
  auto energyStorage = interfaceStorage.energy();
  // Consistency check
  BOOST_CHECK_CLOSE(energy, energyStorage, 1.0e-10);
  BOOST_CHECK_CLOSE(energy, referenceEnergy, 1.0e-10);
  boost::filesystem::remove_all("tmpIPI");
}