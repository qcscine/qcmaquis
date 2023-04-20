/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE ElectronicFEAST

#include <boost/mpl/list.hpp>
#include <boost/test/included/unit_test.hpp>
#include "dmrg/MetaSweepSimulations/FEASTLauncher.h"
#include "dmrg/MetaSweepSimulations/FEASTSimulator.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/sim/matrix_types.h"
#include "Fixtures/H2Fixture.h"
#include "Fixtures/LiHFixture.h"
#include "Fixtures/BenzeneFixture.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
> symmetries;

/** @brief Test FEAST for electronic calculations (we target the ground and first excited state) */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_FEAST_Electronic_H2, S, symmetries, H2Fixture)
{
  using FEASTSimulatorType = FEASTSimulator<S>;
  using ModelType = Model<cmatrix, S>;
  //
  parametersH2.set("max_bond_dimension", 10);
  parametersH2.set("init_type", "default");
  parametersH2.set("seed", 19893003);
  parametersH2.set("optimization", "twosite");
  parametersH2.set("nsweeps", 5);
  parametersH2.set("chkpfile", "GS.H2.FEAST.chkp.h5");
  parametersH2.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
  maquis::DMRGInterface<double> interfaceOptimizerGS(parametersH2);
  interfaceOptimizerGS.optimize();
  auto energyFromOptimizerGS = interfaceOptimizerGS.energy();
  // Excited-state calculation
  parametersH2.set("n_ortho_states", 1);
  parametersH2.set("ortho_states", "GS.H2.FEAST.chkp.h5");
  parametersH2.set("chkpfile", "ES.H2.FEAST.chkp.h5");
  maquis::DMRGInterface<double> interfaceOptimizerES(parametersH2);
  interfaceOptimizerES.optimize();
  auto energyFromOptimizerES = interfaceOptimizerES.energy();
  // Cleans up stuff
  boost::filesystem::remove_all("GS.H2.FEAST.chkp.h5");
  boost::filesystem::remove_all("ES.H2.FEAST.chkp.h5");
  // FEAST for ground state
  auto eMin = energyFromOptimizerGS - (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  auto eMax = energyFromOptimizerGS + (energyFromOptimizerES-energyFromOptimizerGS)/10.;
  parametersH2.set("feast_num_states", 1);
  parametersH2.set("feast_max_iter", 1);
  parametersH2.set("feast_emin", eMin);
  parametersH2.set("feast_emax", eMax);
  parametersH2.set("feast_num_points", 8);
  parametersH2.set("init_type", "default");
  parametersH2.set("linsystem_precond", "no");
  parametersH2.set("linsystem_krylov_dim", 50);
  parametersH2.set("linsystem_tol", 1.0E-5);
  parametersH2.set("linsystem_init", "last");
  auto H2Lattice = Lattice(parametersH2);
  auto H2Model = ModelType(H2Lattice, parametersH2);
  auto H2MPO = make_mpo(H2Lattice, H2Model);
  auto feastSimulator = FEASTSimulatorType(parametersH2, H2Model, H2Lattice, H2MPO);
  feastSimulator.runFEAST();
  auto feastEnergy = feastSimulator.getEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromOptimizerGS, 1.0E-6);
  // FEAST for excited state
  eMin = energyFromOptimizerES - 0.0001;
  eMax = energyFromOptimizerES + 0.0001;
  parametersH2.set("feast_max_iter", 5);
  parametersH2.set("feast_emin", eMin);
  parametersH2.set("feast_emax", eMax);
  auto feastSimulatorES = FEASTSimulatorType(parametersH2, H2Model, H2Lattice, H2MPO);
  feastSimulatorES.runFEAST();
  feastEnergy = feastSimulatorES.getEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromOptimizerES, 1.0E-6);
}

/** @brief Test FEAST for electronic calculations (we target the ground and first excited state) */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_FEAST_Electronic_H2_Interface, S, symmetries, H2Fixture)
{
  using ComplexType = std::complex<double>;
  using FEASTSimulatorType = FEASTSimulator<S>;
  using ModelType = Model<cmatrix, S>;
  //
  parametersH2.set("max_bond_dimension", 10);
  parametersH2.set("optimization", "twosite");
  parametersH2.set("nsweeps", 5);
  parametersH2.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
  parametersH2.set("feast_num_states", 1);
  parametersH2.set("feast_max_iter", 1);
  parametersH2.set("feast_emin", -0.981);
  parametersH2.set("feast_emax", -0.98);
  parametersH2.set("feast_num_points", 8);
  parametersH2.set("init_type", "const");
  parametersH2.set("linsystem_precond", "no");
  parametersH2.set("linsystem_krylov_dim", 50);
  parametersH2.set("linsystem_tol", 1.0E-10);
  parametersH2.set("linsystem_init", "last");
  auto H2Lattice = Lattice(parametersH2);
  auto H2Model = ModelType(H2Lattice, parametersH2);
  auto H2MPO = make_mpo(H2Lattice, H2Model);
  auto feastSimulator = FEASTSimulatorType(parametersH2, H2Model, H2Lattice, H2MPO);
  feastSimulator.runFEAST();
  auto feastEnergy = maquis::real(feastSimulator.getEnergy(0));
  // Gets the energy via the interface
  maquis::DMRGInterface<ComplexType> interfaceFEAST(parametersH2);
  interfaceFEAST.runFEAST();
  auto feastEnergyInterface = maquis::real(interfaceFEAST.energyFEAST(0));
  // Final checks
  BOOST_CHECK_CLOSE(feastEnergy, feastEnergyInterface, 1.0E-6);
}

/** @brief Test that running FEAST with real-valued parameters raises an exception */
BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_FEAST_Electronic_H2_Real, S, symmetries, H2Fixture)
{
  parametersH2.set("max_bond_dimension", 10);
  parametersH2.set("init_type", "default");
  parametersH2.set("seed", 19893003);
  parametersH2.set("optimization", "twosite");
  parametersH2.set("nsweeps", 5);
  parametersH2.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
  maquis::DMRGInterface<double> interfaceFEASTReal(parametersH2);
  BOOST_CHECK_THROW(interfaceFEASTReal.runFEAST(), FEASTException);
}

 #ifdef HAVE_TwoU1PG

/** @brief Test FEAST wit two guesses */
BOOST_FIXTURE_TEST_CASE(Test_FEAST_Electronic_LiH_TwoGuesses, LiHFixture)
{
  // Generic data
  using ComplexType = std::complex<double>;
  auto refEnergyGS = -7.904357504731657;
  // Generic parameters
  parametersLiH.set("max_bond_dimension", 50);
  parametersLiH.set("optimization", "twosite");
  parametersLiH.set("symmetry", "2u1pg");
  parametersLiH.set("seed", 4);
  parametersLiH.set("nsweeps", 5);
  parametersLiH.set("nmainsweeps", 2);
  parametersLiH.set("ngrowsweeps", 2);
  parametersLiH.set("truncation_initial", 1.0E-30);
  parametersLiH.set("truncation_final", 1.0E-30);
  parametersLiH.set("linsystem_exact_error", "yes");
  // Linear system parameters (note that the same set of parameters is used also for DMRG[FEAST])
  // parametersLiH.set("linsystem_precond", "yes");
  parametersLiH.set("linsystem_init", "last");
  parametersLiH.set("linsystem_max_it", 1);
  parametersLiH.set("linsystem_tol", 1.0E-10);
  parametersLiH.set("linsystem_krylov_dim", 50);
  parametersLiH.set("linsystem_exact_error", "yes");
  parametersLiH.set("linsystem_verbose", "no");
  // FEAST parameters
  parametersLiH.set("feast_max_iter", 2);
  parametersLiH.set("feast_num_points", 8);
  parametersLiH.set("init_type", "default");
  parametersLiH.set("feast_overlap_convergence_threshold", 1.0E-5);
  parametersLiH.set("feast_energy_convergence_threshold", 1.0E-6);
  // Note that the interval includes two states, but we use three guesses.
  parametersLiH.set("feast_emin", refEnergyGS-0.001);
  parametersLiH.set("feast_emax", refEnergyGS+0.001);
  parametersLiH.set("feast_num_states", 2);
  parametersLiH.set("feast_calculate_standard_deviation", "yes");
  // Constructs the interface and runs FEAST
  maquis::DMRGInterface<ComplexType> interfaceFEASTGS(parametersLiH);
  interfaceFEASTGS.runFEAST();
  auto energyGS = maquis::real(interfaceFEASTGS.energyFEAST(0));
  BOOST_CHECK_CLOSE(energyGS, refEnergyGS, 1.0E-10);
}

#endif // HAVE_TwoU1PG

/** 
 * @brief Test FEAST on Benzene with two guesses and an interval comprising 1 state 
 * Note that the threshold used in the test is rather high (1.0E-5) because we use 
 * a bond dimension (50) which does not allow to match precisely the reference energy,
 * which was obtained from a fully converged DMRG calculation.
 */
BOOST_FIXTURE_TEST_CASE(Test_FEAST_Electronic_Benzene_TwoGuesses, BenzeneFixture)
{
  // Generic data
  using ComplexType = std::complex<double>;
  auto referenceEnergy = -230.7555388354673;
  // Generic parameters
  parametersBenzene.set("max_bond_dimension", 50);
  parametersBenzene.set("optimization", "twosite");
  parametersBenzene.set("seed", 42);
  parametersBenzene.set("nsweeps", 2);
  parametersBenzene.set("truncation_initial", 1.0E-30);
  parametersBenzene.set("truncation_final", 1.0E-30);
  // parametersBenzene.set("linsystem_exact_error", "yes");
  // Linear system parameters (note that the same set of parameters is used also for DMRG[FEAST])
  parametersBenzene.set("linsystem_init", "last");
  parametersBenzene.set("linsystem_max_it", 1);
  parametersBenzene.set("linsystem_tol", 1.0E-8);
  parametersBenzene.set("linsystem_krylov_dim", 10);
  // FEAST parameters
  parametersBenzene.set("feast_max_iter", 2);
  parametersBenzene.set("feast_num_points", 8);
  parametersBenzene.set("init_type", "default");
  parametersBenzene.set("feast_overlap_convergence_threshold", 1.0E-5);
  parametersBenzene.set("feast_energy_convergence_threshold", 1.0E-6);
  // Note that the interval includes two states, but we use three guesses.
  parametersBenzene.set("feast_emin", -230.8);
  parametersBenzene.set("feast_emax", -230.7);
  parametersBenzene.set("feast_num_states", 2);
#ifdef HAVE_TwoU1
  // Constructs the interface and runs FEAST.
  // Only one state is included in the interval.
  parametersBenzene.set("symmetry", "2u1");
  parametersBenzene.set("linsystem_exact_error", "yes");
  parametersBenzene.set("feast_calculate_standard_eviation", "no");
  // Constructs the interface and runs FEAST
  maquis::DMRGInterface<ComplexType> interfaceBenzene(parametersBenzene);
  interfaceBenzene.runFEAST();
  auto energy1 = maquis::real(interfaceBenzene.energyFEAST(0));
  BOOST_CHECK_CLOSE(energy1, referenceEnergy, 1.0E-5);
#endif // HAVE_TwoU1
  // Does the same for the point group case
#ifdef HAVE_TwoU1PG
  parametersBenzene.set("symmetry", "2u1pg");
  maquis::DMRGInterface<ComplexType> interfaceBenzenePointGroup(parametersBenzene);
  interfaceBenzenePointGroup.runFEAST();
  auto energy2 = maquis::real(interfaceBenzenePointGroup.energyFEAST(0));
  BOOST_CHECK_CLOSE(energy2, referenceEnergy, 1.0E-5);
#endif // HAVE_TwoU1PG
}

/**
 * @brief Test FEAST on Benzene with 3 guesses and an interval comprising 2 states.
 * Note that, to differentiate the test compared to the previous test case, we use
 * the symmetry-adapted code.
 * Note that, as for the test above, the check threshold is low because -- to keep
 * the test short -- we use a relatively low value for the bond dimension.
 */
BOOST_FIXTURE_TEST_CASE(Test_FEAST_Electronic_Benzene_ThreeGuesses, BenzeneFixture)
{
  // Generic data
  using ComplexType = std::complex<double>;
  auto referenceEnergy0 = -230.7555388354674;
  auto referenceEnergy1 = -230.4936880170761;
  // Generic parameters
  parametersBenzene.set("max_bond_dimension", 50);
  parametersBenzene.set("optimization", "twosite");
  parametersBenzene.set("seed", 42);
  parametersBenzene.set("nsweeps", 5);
  parametersBenzene.set("truncation_initial", 1.0E-30);
  parametersBenzene.set("truncation_final", 1.0E-30);
  parametersBenzene.set("linsystem_exact_error", "no");
  // Linear system parameters (note that the same set of parameters is used also for DMRG[FEAST])
  parametersBenzene.set("linsystem_init", "last");
  parametersBenzene.set("linsystem_max_it", 1);
  parametersBenzene.set("linsystem_tol", 1.0E-8);
  parametersBenzene.set("linsystem_krylov_dim", 10);
  // FEAST parameters
  parametersBenzene.set("feast_max_iter", 3);
  parametersBenzene.set("feast_num_points", 8);
  parametersBenzene.set("init_type", "default");
  // parametersBenzene.set("init_basis_state", "4,4,4,1,1,1|4,4,1,4,1,1|4,1,4,1,4,1");
  parametersBenzene.set("feast_overlap_convergence_threshold", 1.0E-5);
  parametersBenzene.set("feast_energy_convergence_threshold", 1.0E-6);
  // Note that the interval includes two states, but we use three guesses.
  parametersBenzene.set("feast_emin", -230.80);
  parametersBenzene.set("feast_emax", -230.48);
  parametersBenzene.set("feast_num_states", 3);
#ifdef HAVE_SU2U1
  // Constructs the interface and runs FEAST.
  // Only one state is included in the interval.
  parametersBenzene.set("symmetry", "su2u1");
  maquis::DMRGInterface<ComplexType> interfaceBenzene(parametersBenzene);
  interfaceBenzene.runFEAST();
  std::vector<double> vectorOFEnergies = {maquis::real(interfaceBenzene.energyFEAST(0)),
                                          maquis::real(interfaceBenzene.energyFEAST(1))};
  std::sort(vectorOFEnergies.begin(), vectorOFEnergies.end());
  // Now sets the intervals to be smaller.
  parametersBenzene.set("feast_emin", vectorOFEnergies[0]-0.001);
  parametersBenzene.set("feast_emax", vectorOFEnergies[0]+0.001);
  maquis::DMRGInterface<ComplexType> interfaceBenzeneSmaller(parametersBenzene);
  interfaceBenzene.runFEAST();
  BOOST_CHECK_CLOSE(maquis::real(interfaceBenzene.energyFEAST(0)), referenceEnergy0, 1.0E-8);
  // BOOST_CHECK_CLOSE(vectorOFEnergies[1], referenceEnergy1, 1.0E-8);
#endif // HAVE_SU2U1
  // Does the same for the point group case
#ifdef HAVE_SU2U1PG
  parametersBenzene.set("symmetry", "su2u1pg");
  // We reset the upper and lower bounds to the "loose" values
  parametersBenzene.set("feast_emin", -230.80);
  parametersBenzene.set("feast_emax", -230.48);
  maquis::DMRGInterface<ComplexType> interfaceBenzenePG(parametersBenzene);
  interfaceBenzenePG.runFEAST();
  std::vector<double> vectorOFEnergiesPG = {maquis::real(interfaceBenzenePG.energyFEAST(0)),
                                            maquis::real(interfaceBenzenePG.energyFEAST(1))};
  std::sort(vectorOFEnergiesPG.begin(), vectorOFEnergiesPG.end());
  // Now sets the intervals to be smaller.
  parametersBenzene.set("feast_emin", vectorOFEnergiesPG[1]-0.001);
  parametersBenzene.set("feast_emax", vectorOFEnergiesPG[1]+0.001);
  maquis::DMRGInterface<ComplexType> interfaceBenzeneSmallerPG(parametersBenzene);
  interfaceBenzeneSmallerPG.runFEAST();
  BOOST_CHECK_CLOSE(maquis::real(interfaceBenzeneSmallerPG.energyFEAST(0)), referenceEnergy1, 1.0E-8);
#endif // HAVE_SU2U1PG
}

#ifdef HAVE_SU2U1PG

/**
 * @brief Test FEAST for electronic calculations (we target the first excited state)
 * Note that, unlike the other test-cases, here we benchmark against DMRG[IPI] and also
 * repeat the FEAST iterations.
 */
BOOST_FIXTURE_TEST_CASE(Test_FEAST_Electronic_LiH, LiHFixture)
{
  using FEASTSimulatorType = FEASTSimulator<SU2U1PG>;
  using ModelType = Model<cmatrix, SU2U1PG>;
  // Generic parameters
  parametersLiH.set("max_bond_dimension", 50);
  parametersLiH.set("init_type", "hf");
  parametersLiH.set("hf_occ", "4,1,1,1");
  parametersLiH.set("optimization", "twosite");
  parametersLiH.set("symmetry", "su2u1pg");
  parametersLiH.set("init_type", "const");
  // IPI-specific parameters
  parametersLiH.set("ipi_sweep_overlap_threshold", 1.0E-5);
  parametersLiH.set("ipi_sweep_energy_threshold", 1.0E-5);
  parametersLiH.set("ipi_sweeps_per_system", 5);
  parametersLiH.set("ipi_iterations", 8);
  // Linear system parameters (note that the same set of parameters is used also for DMRG[FEAST])
  parametersLiH.set("linsystem_precond", "no");
  parametersLiH.set("linsystem_init", "last");
  parametersLiH.set("linsystem_max_it", 1);
  parametersLiH.set("linsystem_tol", 1.0E-10);
  parametersLiH.set("linsystem_krylov_dim", 20);
  // == IPI SIMULATION ==
  // Ground state, reference taken from test2.cpp == -7.90436
  parametersLiH.set("ipi_shift", -8.);
  maquis::DMRGInterface<double> optimizerIpiGS(parametersLiH);
  optimizerIpiGS.runInversePowerIteration();
  auto energyFromIpiGS = optimizerIpiGS.energy();
  BOOST_CHECK_SMALL(std::abs(energyFromIpiGS - -7.90436), 1.0E-4);
  // First excited-state, reference taken from test2.cpp == -7.77349
  parametersLiH.set("ipi_shift", -7.8);
  maquis::DMRGInterface<double> optimizerIpiES1(parametersLiH);
  optimizerIpiES1.runInversePowerIteration();
  auto energyFromIpiES1 = optimizerIpiES1.energy();
  BOOST_CHECK_SMALL(std::abs(energyFromIpiES1 - -7.77349), 1.0E-4);
  // Second excited state (no reference from test2.cpp, but DMRG data == -7.275314952)
  parametersLiH.set("ipi_shift", -7.3);
  maquis::DMRGInterface<double> optimizerIpiES2(parametersLiH);
  optimizerIpiES2.runInversePowerIteration();
  auto energyFromIpiES2 = optimizerIpiES2.energy();
  BOOST_CHECK_SMALL(std::abs(energyFromIpiES2 - -7.275314952), 1.0E-4);
  // FEAST for ground state
  parametersLiH.set("nsweeps", 5);
  parametersLiH.set("feast_num_states", 1);
  parametersLiH.set("feast_max_iter", 3);
  parametersLiH.set("feast_num_points", 8);
  parametersLiH.set("init_type", "default");
  parametersLiH.set("feast_overlap_convergence_threshold", 1.0E-5);
  parametersLiH.set("feast_energy_convergence_threshold", 1.0E-6);
  //
  parametersLiH.set("linsystem_precond", "no");
  parametersLiH.set("linsystem_krylov_dim", 50);
  parametersLiH.set("linsystem_tol", 1.0E-5);
  parametersLiH.set("linsystem_init", "last");
  // Ground state
  auto eMin = energyFromIpiGS - 0.001;
  auto eMax = energyFromIpiGS + 0.001;
  parametersLiH.set("feast_emin", eMin);
  parametersLiH.set("feast_emax", eMax);
  auto LiHLattice = Lattice(parametersLiH);
  auto LiHModel = ModelType(LiHLattice, parametersLiH);
  auto LiHMPO = make_mpo(LiHLattice, LiHModel);
  auto feastSimulator = FEASTSimulatorType(parametersLiH, LiHModel, LiHLattice, LiHMPO);
  feastSimulator.runFEAST();
  auto feastEnergy = feastSimulator.getEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromIpiGS, 1.0E-5);
  // First excited state
  eMin = energyFromIpiES1 - 0.001;
  eMax = energyFromIpiES1 + 0.001;
  parametersLiH.set("feast_emin", eMin);
  parametersLiH.set("feast_emax", eMax);
  auto feastSimulatorES1 = FEASTSimulatorType(parametersLiH, LiHModel, LiHLattice, LiHMPO);
  feastSimulatorES1.runFEAST();
  feastEnergy = feastSimulatorES1.getEnergy(0);
  BOOST_CHECK_CLOSE(feastEnergy, energyFromIpiES1, 1.0E-5);
  // Simultaneous calculation on ground and excited state (so, 2 roots)
  eMin = energyFromIpiGS - 0.001;
  eMax = energyFromIpiES1 + 0.001;
  parametersLiH.set("feast_emin", eMin);
  parametersLiH.set("feast_emax", eMax);
  parametersLiH.set("feast_num_states", 2);
  auto feastSimulatorTwoStates = FEASTSimulatorType(parametersLiH, LiHModel, LiHLattice, LiHMPO);
  feastSimulatorTwoStates.runFEAST();
  auto feastEnergyGS = feastSimulatorTwoStates.getEnergy(0);
  auto feastEnergyES = feastSimulatorTwoStates.getEnergy(1);
  BOOST_CHECK_CLOSE(feastEnergyGS+feastEnergyES, energyFromIpiGS+energyFromIpiES1, 1.0E-5);
}

#endif // HAVE_SU2U1PG
