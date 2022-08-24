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

#ifndef GENERIC_SWEEPS_SIMULATION_H
#define GENERIC_SWEEPS_SIMULATION_H

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/results_collector.h"
#include "SweepMPSContainer.h"
#include "SweepMPOContainer.h"
#include "SweepMPSUpdater.h"
#include "SweepOptimizationTypeTrait.h"

/**
 * @brief Class representing a generic sweep-based simulation.
 *
 * @tparam Matrix matrix class underlying the matrix storage.
 * @tparam SymmGroup symmetry group of the Hamiltonian.
 */

template<class Matrix, class SymmGroup, class Storage, SweepOptimizationType SweepType>
class GenericSweepSimulation {
public:
  // Types declaration
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using MPSContainerType = SweepMPSContainer<Matrix, SymmGroup, SweepType>;
  using MPOContainerType = SweepMPOContainer<Matrix, SymmGroup, SweepType>;
  using SweepMPSUpdaterType = SweepMPSUpdater<Matrix, SymmGroup, Storage, SweepType>;
  using BoundaryPropagatorType = BoundaryPropagator<Matrix, SymmGroup, Storage>;
  using SweepTraitClass = SweepOptimizationTypeTrait<SweepType>;

  /** @brief Class constructor */
  GenericSweepSimulation(MPSType& mps, const MPOType& mpo, BaseParameters& parms,
                         int initSite=0)
    : mps_(mps), mpo_(mpo), parms_(parms), L_(mps_.length()),
      mpoContainer_(mpo_, mps_), mpsContainer_(mps)
  {
    lastSite_ = SweepTraitClass::getLastSite(L_);
    boundaryPropagator_ = std::make_shared<BoundaryPropagatorType>(mps_, mpo_);
    mpsUpdater_ = std::make_unique<SweepMPSUpdaterType>(mpo_, mps_, boundaryPropagator_, parms_);
  };

  /**
   * @brief Execution of a generic sweep-based optimization algorithm.
   *
   * Note that we delegate every action to the derived class, with the exception of the
   * memory management, which is done here to ensure that
   *
   */
  void runSweepSimulation() {
    // Operations performed at the beginning of the 
    int maxNumberOfSweeps = parms_["nsweeps"];
    // == LOOP OVER THE SWEEPS ==
    for (int iSweep = 0; iSweep < maxNumberOfSweeps; iSweep++) {
      // Preparatory operations.
      maquis::cout << " == SWEEP NUMBER = " << maxNumberOfSweeps << " ==" << std::endl;
      maquis::cout << std::endl;
      this->prepareSweep();
      indexOfMicroIteration_ = 0;
      currentSite_ = SweepTraitClass::convertMicroIterationToSite(L_, indexOfMicroIteration_);
      // Prefetches the boundaries that will be needed for the first sweep
      Storage::prefetch(boundaryPropagator_->getLeftBoundary(SweepTraitClass::getIndexOfLeftBoundary(currentSite_, SweepDirectionType::Forward)));
      Storage::prefetch(boundaryPropagator_->getRightBoundary(SweepTraitClass::getIndexOfRightBoundary(currentSite_, SweepDirectionType::Forward)));
      // == LOOP OVER THE MICROITERATIONS ==
      while (indexOfMicroIteration_ < 2*lastSite_) {
        maquis::cout << " -- Microiteration number = " << indexOfMicroIteration_ << " --" << std::endl;
        maquis::cout << std::endl;
        // Calculates the relevant indices on the DMRG lattice.
        auto sweepType = (indexOfMicroIteration_ < lastSite_) ? SweepDirectionType::Forward : SweepDirectionType::Backward;
        currentSite_ = SweepTraitClass::convertMicroIterationToSite(L_, indexOfMicroIteration_);
        siteLeft_ = SweepTraitClass::getIndexOfLeftBoundary(currentSite_, sweepType);
        siteRight_ = SweepTraitClass::getIndexOfRightBoundary(currentSite_, sweepType);
        mpoContainer_.updatePlacements(indexOfMicroIteration_, siteLeft_);
        // We must be careful here because, for the two-site case, there is the risk of fetching twice the boundaries.
        // In fact, we run the optimization of sites (L-1, L) twice consequently
        //TODO ALB THIS SHOULD BE FIXED!
        if (SweepTraitClass::countEndSiteTwice_ && indexOfMicroIteration_ != lastSite_) {
          Storage::fetch(boundaryPropagator_->getLeftBoundary(siteLeft_));
          Storage::fetch(boundaryPropagator_->getRightBoundary(siteRight_));
        }
        // Starts prefetching what will be needed in the following microiteration.
        // Note that, for instance, we don't prefetch the left boundary for the l2r sweep because
        // this will be taken care in the boundary propagation (in other words, there is no 
        // need to prefetch the left boundary since it will be anyways modified by the boundary
        // propagation)
        if (sweepType == SweepDirectionType::Forward) {
          auto nextIndex = SweepTraitClass::getIndexOfNextRightBoundary(currentSite_, sweepType);
          if (nextIndex <= L_)
            Storage::prefetch(boundaryPropagator_->getRightBoundary(nextIndex));
        }
        else if (sweepType == SweepDirectionType::Backward) {
          auto nextIndex = SweepTraitClass::getIndexOfNextLeftBoundary(currentSite_, sweepType);
          if (nextIndex >= 0)
            Storage::prefetch(boundaryPropagator_->getLeftBoundary(nextIndex));
        }
        // == SOLUTION OF THE LOCAL PROBLEM ==
        this->prepareMicroiteration();
        auto outputTensor = this->solveLocalProblem();
        // == MPS UPDATE ==
        auto truncationResults = mpsUpdater_->updateMPS(siteLeft_, siteRight_, sweepType, outputTensor, this->getAlpha(iSweep),
                                                        this->get_cutoff(iSweep), this->get_Mmax(iSweep));
        // == BOUNDARY PROPAGATION ==
        this->propagateBoundaries();
        // After the boundary propagation we can do two operations at the memory level.
        // 1) we can drop the memory of the right boundary (in the case of a l2r sweep).
        //    In fact, that memory will anyways be overwritten by the r2l sweep that
        //    will follow.
        // 2) We can write to file the left boundary that we are "leaving behind"
        if (sweepType == SweepDirectionType::Forward && siteLeft_ != L_-1) {
          Storage::drop(boundaryPropagator_->getRightBoundary(siteRight_));
          Storage::StoreToFile(boundaryPropagator_->getLeftBoundary(siteLeft_));
        }
        else if (sweepType == SweepDirectionType::Backward && siteLeft_ != 0) {
          Storage::drop(boundaryPropagator_->getLeftBoundary(siteLeft_));
          Storage::StoreToFile(boundaryPropagator_->getRightBoundary(siteRight_));
        }
        this->finalizeMicroIteration(truncationResults);
        indexOfMicroIteration_ += 1;
        maquis::cout << std::endl;
      }
    }
    this->finalizeSweep();
  }

  /** @brief Gets the container with the results of each iteration */
  const auto& iteration_results() const { return iterationResults_; }

  /** @brief Gets a specific value of the iteration result */
  template<class CastType>
  CastType getSpecificResult(std::string resultName) {
    if (!iterationResults_.has(resultName))
      throw std::runtime_error("Trying to access non-existing simulation result");
    return boost::any_cast<CastType>(iterationResults_[resultName].get()[0]);
  }

protected:

  /**
   * @brief Collects the operation to be done before a sweep.
   * Note that these operations are done only once per sweep, i.e. they are not
   * repeated at each microiteration.
   */
  virtual void prepareSweep() = 0;

  /** @brief Collects the operation to be performed before a microiteration */
  virtual void prepareMicroiteration() = 0;

  /** @brief Runs the actual sweep simulation */
  virtual MPSTensorType solveLocalProblem() = 0;

  /** @brief Boundary propagation method */
  virtual void propagateBoundaries() = 0;

  /** @brief Collects the operation to be run at the end of a micro iteration */
  virtual void finalizeMicroIteration(const truncation_results& trunc) = 0;

  /** @brief Collects the operation to be run at the end of the simulation */
  virtual void finalizeSweep() = 0;

  /** @brief Simple utility function for a logarithmic interpolation */
  static double log_interpolate(double y0, double y1, int N, int i)
  {
    double ret;
    if (N < 2)
      return y1;
    if (y0 == 0)
      return 0;
    double x = log(y1/y0)/(N-1);
    return y0*exp(x*i);
  }

  /** @brief Method to get the truncation threshold for a given sweep */
  double get_cutoff(int sweep) const {
    return (sweep >= parms_.template get<int>("ngrowsweeps")) ? parms_.template get<double>("truncation_final")
        : log_interpolate(parms_.template get<double>("truncation_initial"),
                          parms_.template get<double>("truncation_final"),
                          parms_.template get<int>("ngrowsweeps"), sweep);
  }

  /** @brief Method to get the proper noise parameter */
  auto getAlpha(int iSweep) const {
    double alpha;
    int ngs = parms_.template get<int>("ngrowsweeps");
    int nms = parms_.template get<int>("nmainsweeps");
    if (iSweep < ngs)
      alpha = parms_.template get<double>("alpha_initial");
    else if (iSweep < ngs + nms)
      alpha = parms_.template get<double>("alpha_main");
    else
      alpha = parms_.template get<double>("alpha_final");
    return alpha;
  }

  /** @brief Calculats the maximum value of the bond dimension for a given sweep */
  int get_Mmax(int sweep) const {
    std::size_t Mmax;
    if (parms_.is_set("sweep_bond_dimensions")) {
      std::vector<std::size_t> ssizes = parms_.template get<std::vector<std::size_t> >("sweep_bond_dimensions");
      if (sweep >= ssizes.size())
        Mmax = *ssizes.rbegin();
      else
        Mmax = ssizes[sweep];
    } else {
      Mmax = parms_.template get<std::size_t>("max_bond_dimension");
    }
    return Mmax;
  }

protected:
  MPSType& mps_;
  const MPOType& mpo_;
  MPOContainerType mpoContainer_;
  MPSContainerType mpsContainer_;
  std::unique_ptr<SweepMPSUpdaterType> mpsUpdater_;
  int initSite_, L_, indexOfMicroIteration_, currentSite_, lastSite_, siteLeft_, siteRight_;
  BaseParameters& parms_;
  results_collector iterationResults_;
  std::shared_ptr<BoundaryPropagatorType> boundaryPropagator_;
};

#endif // GENERIC_SWEEPS_SIMULATION_H
