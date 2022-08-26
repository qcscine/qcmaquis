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

#ifndef SWEEP_SIMULATION_FACTORY
#define SWEEP_SIMULATION_FACTORY

#include <memory>
#include "GenericSweepSimulation.h"
#include "SweepBasedEnergyMinimization.h"
#include "SweepBasedLinearSystem.h"
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
  using LinearSystemSSSimulationType = SweepBasedEnergyMinimization<Matrix, SymmGroup, Storage, SweepOptimizationType::SingleSite>;
  using LinearSystemTSSimulationType = SweepBasedEnergyMinimization<Matrix, SymmGroup, Storage, SweepOptimizationType::TwoSite>;
  using PointerToSSSimulatorType = std::unique_ptr<GenericSSSimulationType>;
  using PointerToTSSimulatorType = std::unique_ptr<GenericTSSimulationType>;
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;

public:
  SweepSimulationFactory(std::string simulationName, SweepOptimizationType sweepType,
                         MPSType& mps, const MPOType& mpo, BaseParameters& parms, int initSite)
    : sweepType_(sweepType)
  {
    // Optimization
    if (simulationName == "optimize")
      if (sweepType_ == SweepOptimizationType::SingleSite)
        ssSimulator_ = std::make_unique<OptimizationSSSimulationType>(mps, mpo, parms, initSite);
      else if (sweepType_ == SweepOptimizationType::TwoSite)
        tsSimulator_ = std::make_unique<OptimizationTSSimulationType>(mps, mpo, parms, initSite);
    // Solution of a linear system
    if (simulationName == "linear_system")
      if (sweepType_ == SweepOptimizationType::SingleSite)
        ssSimulator_ = std::make_unique<LinearSystemSSSimulationType>(mps, mpo, parms, initSite);
      else if (sweepType_ == SweepOptimizationType::TwoSite)
        tsSimulator_ = std::make_unique<LinearSystemTSSimulationType>(mps, mpo, parms, initSite);
    if (!ssSimulator_ && !tsSimulator_)
      throw std::runtime_error("Error in parameters for [SweepSimulationFactory] object");
  };

  void runSweepSimulation() {
    if (ssSimulator_)
      ssSimulator_->runSweepSimulation();
    else
      tsSimulator_->runSweepSimulation();
  }


  void runSingleSweep(int iSweep) {
    if (ssSimulator_)
      ssSimulator_->runSingleSweep(iSweep);
    else
      tsSimulator_->runSingleSweep(iSweep);
  }

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
