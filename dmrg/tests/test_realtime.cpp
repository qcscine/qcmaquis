/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include "utils/io.hpp"
#include <iostream>
#include "maquis_dmrg.h"
#include "Fixtures/H2Fixture.h"
#include "Fixtures/BenzeneFixture.h"
#include "Fixtures/TimeEvolversFixture.h"
#include "Fixtures/PreBOTimeEvolversFixture.h"
#include "Fixtures/VibronicFixture.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/sim/matrix_types.h"

/**
 * @brief Tests that the energy is conserved along a TD-DMRG propagation.
 *
 * The data are obtained for CAS(6, 6) and based on the cc-pVDZ basis set.
 * The data are stored in the BenzeneFixture class.
 */
BOOST_FIXTURE_TEST_CASE(TestRealTime, BenzeneFixture)
{
  std::vector<std::string> symmetries;
  #ifdef HAVE_SU2U1PG
  symmetries.push_back("su2u1pg");
  #endif
  #ifdef HAVE_SU2U1
  symmetries.push_back("su2u1");
  #endif
  #ifdef HAVE_TwoU1PG
  symmetries.push_back("2u1pg");
  #endif
  #ifdef HAVE_TwoU1
  symmetries.push_back("2u1");
  #endif
  parametersBenzeneRealTime.set("storagedir", "tmpTD");
  for (auto&& s: symmetries) {
    parametersBenzeneRealTime.set("symmetry", s);
    // Single-site evolution
    maquis::cout << "Running SS real-time evolution test for symmetry " << s << std::endl;
    parametersBenzeneRealTime.set("optimization", "singlesite");
    maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersBenzeneRealTime);
    auto initialEnergy = std::real(interfaceSS.energy());
    interfaceSS.evolve();
    auto finalEnergy = std::real(interfaceSS.energy());
    boost::filesystem::remove_all("tmpTD");
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-8);
    // Two-site evolution
    maquis::cout << "Running TS real-time evolution test for symmetry " << s << std::endl;
    parametersBenzeneRealTime.set("optimization", "twosite");
    maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersBenzeneRealTime);
    initialEnergy = std::real(interfaceTS.energy());
    interfaceTS.evolve();
    finalEnergy = std::real(interfaceTS.energy());
    boost::filesystem::remove_all("tmpTD");
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
  }
}


/**
 * @brief Tests non-trivial dynamics in H2.
 *
 * Here, we test that, if one takes the Hamiltonian of H2 in a minimal basis
 * and prepares the initial state to be the HF state, then the dynamics "populates"
 * the doubly-excited state, and the dynamics follows a Rabi oscillation pattern.
 */
#ifdef HAVE_TwoU1PG
BOOST_FIXTURE_TEST_CASE(TestRabiOscillation, H2Fixture)
{
  // General initialization
  auto latticeH2 = Lattice(parametersH2);
  auto modelHF = Model<matrix, TwoU1PG>(latticeH2, parametersH2);
  auto mpo = make_mpo(latticeH2, modelHF);
  // Check first that the HF state and the doubly excited state of H2 do not interact
  // with the single excitation (for the HF determinant, they should because of the Brillouin
  // theorem. For the doubly-excited state, they should because of symmetry).
  parametersH2.set("init_type", "basis_state_generic");
  parametersH2.set("init_basis_state", "4,1");
  auto mpsHF = MPS<matrix, TwoU1PG>(latticeH2.size(), *(modelHF.initializer(latticeH2, parametersH2)));
  // Single excitations
  parametersH2.set("init_basis_state", "2,3");
  auto mpsSingle1 = MPS<matrix, TwoU1PG>(latticeH2.size(), *(modelHF.initializer(latticeH2, parametersH2)));
  parametersH2.set("init_basis_state", "3,2");
  auto mpsSingle2 = MPS<matrix, TwoU1PG>(latticeH2.size(), *(modelHF.initializer(latticeH2, parametersH2)));
  // Double excitation
  parametersH2.set("init_basis_state", "1,4");
  auto mpsDouble = MPS<matrix, TwoU1PG>(latticeH2.size(), *(modelHF.initializer(latticeH2, parametersH2)));
  // Does the actual checks
  BOOST_CHECK_SMALL(expval(mpsHF, mpsSingle1, mpo), 1.0E-12);
  BOOST_CHECK_SMALL(expval(mpsHF, mpsSingle2, mpo), 1.0E-12);
  BOOST_CHECK_SMALL(expval(mpsDouble, mpsSingle1, mpo), 1.0E-12);
  BOOST_CHECK_SMALL(expval(mpsDouble, mpsSingle2, mpo), 1.0E-12);
}
#endif // HAVE_TwoU1PG

#ifdef HAVE_U1DG

/**
 * @brief Tests that the energy is conserved along a relativistic TD-DMRG propagation.
 * The data are obtained for N2+ and the 3-21G basis set.
 */
BOOST_FIXTURE_TEST_CASE(TestRealTimeRelativistic, TestTimeEvolverFixture)
{
  // Generic settings
  parametersRelativistic.set("imaginary_time", "no");
  parametersRelativistic.set("TD_backpropagation", "yes");
  parametersRelativistic.set("time_units", "as");
  parametersRelativistic.set("time_step", 1.0E+00);
  // Single-site evolution
  parametersRelativistic.set("optimization", "singlesite");
  maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersRelativistic);
  auto initialEnergy = std::real(interfaceSS.energy());
  interfaceSS.evolve();
  auto finalEnergy = std::real(interfaceSS.energy());
  BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
  // Two-site evolutions
  parametersRelativistic.set("optimization", "twosite");
  maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersRelativistic);
  initialEnergy = std::real(interfaceTS.energy());
  interfaceTS.evolve();
  finalEnergy = std::real(interfaceTS.energy());
  BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
}

#endif

#ifdef DMRG_PREBO

/**
 * @brief Tests that the energy is conserved along a PreBO TD-DMRG propagation.
 *
 * Note that this example uses the BO hamiltonian for H2 in a minimal basis set - i.e.,
 * only electrons are included.
 */
BOOST_FIXTURE_TEST_CASE( TestRealTimePreBOWithBOHamiltonian, PreBOTestTimeEvolverFixture )
{
  // Generic settings
  parametersBOComplex.set("imaginary_time", "no");
  parametersBOComplex.set("TD_backpropagation", "yes");
  parametersBOComplex.set("time_units", "as");
  parametersBOComplex.set("time_step", 1.0E+00);
  // Single-site evolution
  maquis::cout << "Running SS real-time evolution test for PreBO model " << std::endl;
  parametersBOComplex.set("optimization", "singlesite");
  maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersBOComplex);
  auto initialEnergy = std::real(interfaceSS.energy());
  interfaceSS.evolve();
  auto finalEnergy = std::real(interfaceSS.energy());
  BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
  // Two-site evolutions
  maquis::cout << "Running TS real-time evolution test for PreBO model" << std::endl;
  parametersBOComplex.set("optimization", "twosite");
  maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersBOComplex);
  initialEnergy = std::real(interfaceTS.energy());
  interfaceTS.evolve();
  finalEnergy = std::real(interfaceTS.energy());
  BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
}

/**
 * @brief Tests that the energy is conserved along a "true" PreBO TD-DMRG propagation.
 */
BOOST_FIXTURE_TEST_CASE( TestRealTimePreBO, PreBOTestTimeEvolverFixture )
{
  // Generic settings
  parametersPreBOComplex.set("imaginary_time", "no");
  parametersPreBOComplex.set("TD_backpropagation", "yes");
  parametersPreBOComplex.set("time_units", "as");
  parametersPreBOComplex.set("time_step", 1.0E+00);
  // Single-site evolution
  maquis::cout << "Running SS real-time evolution test for PreBO model " << std::endl;
  parametersPreBOComplex.set("optimization", "singlesite");
  maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersPreBOComplex);
  auto initialEnergy = std::real(interfaceSS.energy());
  interfaceSS.evolve();
  auto finalEnergy = std::real(interfaceSS.energy());
  BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
  // Two-site evolutions
  maquis::cout << "Running TS real-time evolution test for PreBO model" << std::endl;
  parametersPreBOComplex.set("optimization", "twosite");
  maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersPreBOComplex);
  initialEnergy = std::real(interfaceTS.energy());
  interfaceTS.evolve();
  finalEnergy = std::real(interfaceTS.energy());
  BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
}

#endif // DMRG_PREBO

#ifdef DMRG_VIBRONIC

/** @brief Tests that the energy is conserved along an Excitonic TD-DMRG propagation. */
BOOST_FIXTURE_TEST_CASE(TestRealTimeExcitonic, VibronicFixture)
{
#ifdef HAVE_U1
  // Generic settings
  parametersExcitonicAggregateTwoSites.set("imaginary_time", "no");
  parametersExcitonicAggregateTwoSites.set("TD_backpropagation", "yes");
  parametersExcitonicAggregateTwoSites.set("TD_noise", "yes");
  parametersExcitonicAggregateTwoSites.set("propagator_maxiter", 30);
  parametersExcitonicAggregateTwoSites.set("propagator_accuracy", 1.0E-10);
  parametersExcitonicAggregateTwoSites.set("time_units", "as");
  parametersExcitonicAggregateTwoSites.set("time_step", 0.01);
  parametersExcitonicAggregateTwoSites.set("nsweeps", 2);
  parametersExcitonicAggregateTwoSites.set("max_bond_dimension", 100);
  parametersExcitonicAggregateTwoSites.set("init_type", "basis_state_generic");
  parametersExcitonicAggregateTwoSites.set("init_basis_state", "1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
  parametersExcitonicAggregateTwoSites.set("alpha_initial", 1.0E-15);
  parametersExcitonicAggregateTwoSites.set("alpha_main", 1.0E-20);
  parametersExcitonicAggregateTwoSites.set("alpha_final", 0.);
  // Single-site evolution
  maquis::cout << "Running SS real-time evolution test for Excitonic model " << std::endl;
  parametersExcitonicAggregateTwoSites.set("optimization", "singlesite");
  maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersExcitonicAggregateTwoSites);
  auto initialEnergy = std::real(interfaceSS.energy());
  interfaceSS.evolve();
  auto finalEnergy = std::real(interfaceSS.energy());
  BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
  // Two-site evolutions
  maquis::cout << "Running TS real-time evolution test for Excitonic model" << std::endl;
  parametersExcitonicAggregateTwoSites.set("optimization", "twosite");
  maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersExcitonicAggregateTwoSites);
  initialEnergy = std::real(interfaceTS.energy());
  interfaceTS.evolve();
  finalEnergy = std::real(interfaceTS.energy());
  BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
#endif // HAVE_U1
}

#endif // DMRG_PREBO

#ifdef DMRG_VIBRONIC

/** @brief Tests that the energy is conserved along an Excitonic TD-DMRG propagation. */
BOOST_FIXTURE_TEST_CASE( TestRealTimeExcitonic, VibronicFixture )
{
#ifdef HAVE_U1
    // Generic settings
    parametersExcitonicAggregateTwoSites.set("imaginary_time", "no");
    parametersExcitonicAggregateTwoSites.set("TD_backpropagation", "yes");
    parametersExcitonicAggregateTwoSites.set("propagator_maxiter", 30);
    parametersExcitonicAggregateTwoSites.set("propagator_accuracy", 1.0E-10);
    parametersExcitonicAggregateTwoSites.set("time_units", "as");
    parametersExcitonicAggregateTwoSites.set("time_step", 0.01);
    parametersExcitonicAggregateTwoSites.set("nsweeps", 2);
    parametersExcitonicAggregateTwoSites.set("max_bond_dimension", 50);
    parametersExcitonicAggregateTwoSites.set("init_state", "basis_state_generic");
    parametersExcitonicAggregateTwoSites.set("init_basis_state", "1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
    // Single-site evolution
    maquis::cout << "Running SS real-time evolution test for Excitonic model " << std::endl;
    parametersExcitonicAggregateTwoSites.set("optimization", "singlesite");
    maquis::DMRGInterface<std::complex<double>> interfaceSS(parametersExcitonicAggregateTwoSites);
    auto initialEnergy = std::real(interfaceSS.energy());
    interfaceSS.evolve();
    auto finalEnergy = std::real(interfaceSS.energy());
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
    // Two-site evolutions
    maquis::cout << "Running TS real-time evolution test for Excitonic model" << std::endl;
    parametersExcitonicAggregateTwoSites.set("optimization", "twosite");
    maquis::DMRGInterface<std::complex<double>> interfaceTS(parametersExcitonicAggregateTwoSites);
    initialEnergy = std::real(interfaceTS.energy());
    interfaceTS.evolve();
    finalEnergy = std::real(interfaceTS.energy());
    BOOST_CHECK_CLOSE(initialEnergy, finalEnergy, 1.0E-10);
#endif // HAVE_U1
}

#endif // DMRG_VIBRONIC
