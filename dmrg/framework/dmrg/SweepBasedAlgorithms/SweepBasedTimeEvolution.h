/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SWEEP_BASED_TIME_EVOLUTION_H
#define SWEEP_BASED_TIME_EVOLUTION_H

#include "GenericSweepSimulation.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/evolve/TimeEvolvers/timeevolver.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/utils/storage.h"
#include "dmrg/utils/time_limit_exception.h"
#include "dmrg/utils/checks.h"
#include "BoundaryPropagator.h"

#ifdef DMRG_TD

template<class Matrix, class SymmGroup, class Storage, SweepOptimizationType SweepType>
class SweepBasedTimeEvolution : public GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepType> {
public:
  using Base = GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepType>;
  using SweepTraitClass = SweepOptimizationTypeTrait<SweepType>;
  using SiteProblemType = SiteProblem<Matrix, SymmGroup>;
  using MPSType = typename Base::MPSType;
  using ModelType = typename Base::ModelType;
  using MPOType = typename Base::MPOType;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using TimeEvolverType = TimeEvolver<Matrix, SymmGroup, BaseParameters>;
  using ValueType = typename MPSTensorType::value_type;
  //
  using Base::boundaryPropagator_;
  using Base::getAlpha;
  using Base::get_cutoff;
  using Base::get_Mmax;
  using Base::getSpecificResult;
  using Base::indexOfMicroIteration_;
  using Base::isTerminal;
  using Base::iterationResults_;
  using Base::L_;
  using Base::mps_;
  using Base::mpsContainer_;
  using Base::mpsUpdater_;
  using Base::mpoContainer_;
  using Base::parms_;
  using Base::runSweepSimulation;
  using Base::siteLeft_;
  using Base::siteRight_;

  /** @brief Class constructor */
  SweepBasedTimeEvolution(MPSType& mps, const MPOType& mpo, BaseParameters& parms, const ModelType& model,
                          const Lattice& lattice, bool verbose)
    : Base(mps, mpo, parms, model, lattice, verbose, std::string("Time Evolution")), perturbMPS_(false), doBackpropagation_(true),
      isImaginaryTime_(false)
  {
    // Generate classes needed for propagation
    timeEvolver_ = std::make_shared<TimeEvolverType>(parms_);
    timeStep_ = timeEvolver_->get_time();
    if (parms_["TD_backpropagation"] == "no")
      doBackpropagation_ = false;
    if (parms_["TD_noise"] == "yes")
      perturbMPS_ = true;
    if (parms_["imaginary_time"] == "yes")
      isImaginaryTime_ = true;
    // == TODO Move this in the [GenericSweepBasedSimulation] part ==
    // perturber_ = std::make_shared< PerturberType >(left_, right_, mpo_, parms_);
  }

  /** @brief Method called at the beginning of each sweep */
  void prepareSweep() override final {
    iterationResults_.clear();
  }

  /** @brief Method called before each microiteration */
  void prepareMicroiteration() override final {
    siteProblem_ = std::make_unique<SiteProblemType>(boundaryPropagator_->getLeftBoundary(siteLeft_), boundaryPropagator_->getRightBoundary(siteRight_),
                                                     mpoContainer_.getMPOTensor(siteLeft_));
  }

  /** @brief General method for performing the back-propagation step (specialized later) */
  void performBackPropagation(GrowBoundaryModality boundaryGrowthModality) override final {
    if (doBackpropagation_)
      mpsUpdater_->performBackPropagation(boundaryGrowthModality, siteLeft_, siteRight_, timeEvolver_);
  }

  /** @brief Propagation of the MPS for a given site */
  MPSTensorType solveLocalProblem() override final {
    MPSTensorType mpsToPropagate = mpsContainer_.getMPSTensor(siteLeft_);
    timeEvolver_->evolve(*(siteProblem_.get()), mpsToPropagate, true, isTerminal());
    // Note that the energy is calculated only once per sweep -- it will (or should) anyways be conserved,
    // so it's not useful to print the value of the energy for each microiteration.
    if (siteLeft_ == 0 || isImaginaryTime_) {
      auto energy = ietl::get_energy(*(siteProblem_.get()), mpsToPropagate) + maquis::real(mpoContainer_.getMPO().getCoreEnergy());
      resultOfLocalSiteProblem_.first = energy;
      maquis::cout << std::setprecision(10) << " Energy = " << std::setprecision(16) << resultOfLocalSiteProblem_.first << std::endl;
    }
    iterationResults_["Energy"] << resultOfLocalSiteProblem_.first;
    resultOfLocalSiteProblem_.second = mpsToPropagate;
    return resultOfLocalSiteProblem_.second;
    // // Loads the final results
    // auto energy = resultOfLocalSiteProblem_.first + mpoContainer_.getMPO().getCoreEnergy();
    // maquis::cout << std::setprecision(10) << " Energy = " << std::setprecision(16) << energy << std::endl;
    // iterationResults_["Energy"] << energy;
    // return resultOfLocalSiteProblem_.second;
  }

  /** @brief Operations to be executed at the end of a microiteration */
  void finalizeMicroIteration(const truncation_results& trunc) override final {
    iterationResults_["BondDimension"]   << trunc.bond_dimension;
    iterationResults_["TruncatedWeight"] << trunc.truncated_weight;
    iterationResults_["SmallestEV"]      << trunc.smallest_ev;
    // If has reached the end of the chain, increases the time
    if (SweepTraitClass::changeDirectionNextMicroiteration(L_, indexOfMicroIteration_) ||
        indexOfMicroIteration_ == SweepTraitClass::getNumberOfMicroiterations(L_))
    timeEvolver_->add_to_current_time(timeStep_);
  }

  /** @brief Operations to be executed at the end of the sweep */
  void finalizeSweep() override final { }

  /** @brief Whether to normalize the MPS at the end of a half-sweep */
  bool normalizeAtEnd() override final {
    return true;
  }

  /** @brief Whether to appy the noise-based perturbation */
  bool activatePerturbation() override final {
    return perturbMPS_;
  }

private:
  // Class members
  std::shared_ptr< TimeEvolver< Matrix, SymmGroup, BaseParameters > > timeEvolver_;
  std::unique_ptr<SiteProblemType> siteProblem_;
  std::pair<double, MPSTensorType > resultOfLocalSiteProblem_;
  bool doBackpropagation_, perturbMPS_, isImaginaryTime_;
  double timeStep_;
  // std::shared_ptr< SiteShifter< Matrix, SymmGroup, TimeEvolver<Matrix, SymmGroup, BaseParameters>, PerturberType > > site_shifter_;
  // std::shared_ptr<PerturberType> perturber_;
};

#endif // DMRG_TD

#endif // SWEEP_BASED_TIME_EVOLUTION_H
