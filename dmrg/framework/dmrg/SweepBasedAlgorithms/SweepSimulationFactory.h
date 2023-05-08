/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SWEEP_SIMULATION_FACTORY
#define SWEEP_SIMULATION_FACTORY

#include <memory>
#include "GenericSweepSimulation.h"
#include "SweepBasedEnergyMinimization.h"
#include "SweepBasedLinearSystem.h"
#ifdef DMRG_TD
#include "SweepBasedTimeEvolution.h"
#endif // DMRG_TD
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/utils/BaseParameters.h"

template<class Matrix, class SymmGroup, class Storage>
class SweepSimulationFactory {
  // Types declaration
  using GenericSSSimulationType = GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepOptimizationType::SingleSite>;
  using GenericTSSimulationType = GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepOptimizationType::TwoSite>;
  using OptimizationSSSimulationType = SweepBasedEnergyMinimization<Matrix, SymmGroup, Storage, SweepOptimizationType::SingleSite>;
  using OptimizationTSSimulationType = SweepBasedEnergyMinimization<Matrix, SymmGroup, Storage, SweepOptimizationType::TwoSite>;
  using LinearSystemSSSimulationType = SweepBasedLinearSystem<Matrix, SymmGroup, Storage, SweepOptimizationType::SingleSite>;
  using LinearSystemTSSimulationType = SweepBasedLinearSystem<Matrix, SymmGroup, Storage, SweepOptimizationType::TwoSite>;
#ifdef DMRG_TD
  using EvolverSSSimulationType = SweepBasedTimeEvolution<Matrix, SymmGroup, Storage, SweepOptimizationType::SingleSite>;
  using EvolverTSSimulationType = SweepBasedTimeEvolution<Matrix, SymmGroup, Storage, SweepOptimizationType::TwoSite>;
#endif // DMRG_TD
  using PointerToSSSimulatorType = std::unique_ptr<GenericSSSimulationType>;
  using PointerToTSSimulatorType = std::unique_ptr<GenericTSSimulationType>;
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using ModelType = Model<Matrix, SymmGroup>;

public:
  SweepSimulationFactory(std::string simulationName, SweepOptimizationType sweepType, MPSType& mps, const MPOType& mpo,
                         BaseParameters& parms, const ModelType& model, const Lattice& lattice)
    : sweepType_(sweepType)
  {
    // Optimization
    if (simulationName == "optimize") {
      if (sweepType_ == SweepOptimizationType::SingleSite) {
        ssSimulator_ = std::make_unique<OptimizationSSSimulationType>(mps, mpo, parms, model, lattice, true);
      }
      else if (sweepType_ == SweepOptimizationType::TwoSite) {
        tsSimulator_ = std::make_unique<OptimizationTSSimulationType>(mps, mpo, parms, model, lattice, true);
      }
    }
    // Solution of a linear system
    else if (simulationName == "linear_system") {
      bool verbose = parms["linsystem_verbose"] == "yes";
      if (sweepType_ == SweepOptimizationType::SingleSite) {
        ssSimulator_ = std::make_unique<LinearSystemSSSimulationType>(mps, mpo, parms, model, lattice, verbose);
      }
      else if (sweepType_ == SweepOptimizationType::TwoSite) {
        tsSimulator_ = std::make_unique<LinearSystemTSSimulationType>(mps, mpo, parms, model, lattice, verbose);
      }
    }
#ifdef DMRG_TD
    else if (simulationName == "evolve") {
      if (sweepType_ == SweepOptimizationType::SingleSite) {
        ssSimulator_ = std::make_unique<EvolverSSSimulationType>(mps, mpo, parms, model, lattice, true);
      }
      else if (sweepType_ == SweepOptimizationType::TwoSite) {
        tsSimulator_ = std::make_unique<EvolverTSSimulationType>(mps, mpo, parms, model, lattice, true);
      }
    }
#endif // DMRG_TD
    else {
      throw std::runtime_error("Sweep-based simulation type not recognized");
    }
    // Final check
    if (!ssSimulator_ && !tsSimulator_)
      throw std::runtime_error("Error in parameters for [SweepSimulationFactory] object");
  };

  /** @brief Runs a complete sweep-based optimization */
  void runSweepSimulation() {
    if (ssSimulator_)
      ssSimulator_->runSweepSimulation();
    else
      tsSimulator_->runSweepSimulation();
  }

  /** @brief Runs a single sweep (back and forth) */
  void runSingleSweep(int iSweep) {
    if (ssSimulator_)
      ssSimulator_->runSingleSweep(iSweep);
    else
      tsSimulator_->runSingleSweep(iSweep);
  }

  /** @brief Retrieves simulation results */
  auto getIterationResults() {
    if (ssSimulator_)
      return ssSimulator_->iteration_results();
    else
      return tsSimulator_->iteration_results();
  }

private:
  PointerToSSSimulatorType ssSimulator_;
  PointerToTSSimulatorType tsSimulator_;
  SweepOptimizationType sweepType_;
};

#endif // SWEEP_SIMULATION_FACTORY
