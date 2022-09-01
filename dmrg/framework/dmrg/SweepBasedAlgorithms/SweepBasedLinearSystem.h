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

#ifndef SWEEP_BASED_LINEAR_SYSTEM_H
#define SWEEP_BASED_LINEAR_SYSTEM_H

#include "GenericSweepSimulation.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/LinearSystem/linsolver.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/utils/storage.h"
#include "dmrg/utils/time_limit_exception.h"
#include "dmrg/utils/parallel/placement.hpp"
#include "dmrg/utils/checks.h"
#include "BoundaryPropagator.h"
#include "OverlapPropagator.h"

template<class Matrix, class SymmGroup, class Storage, SweepOptimizationType SweepType>
class SweepBasedLinearSystem : public GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepType> {
public:
  using Base = GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepType>;
  using OverlapPropagatorType = OverlapPropagator<Matrix, SymmGroup, Storage>;
  using SweepTraitClass = SweepOptimizationTypeTrait<SweepType>;
  using SiteProblemType = SiteProblem<Matrix, SymmGroup>;
  using LinearSolverType = LinSolver<Matrix, SymmGroup>;
  using MPSType = typename Base::MPSType;
  using MPOType = typename Base::MPOType;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using BlockMatrixType = block_matrix<Matrix, SymmGroup>;
  using ValueType = typename MPSTensorType::value_type;
  //
  using Base::boundaryPropagator_;
  using Base::getSpecificResult;
  using Base::indexOfMicroIteration_;
  using Base::initSite_;
  using Base::iterationResults_;
  using Base::lastSite_;
  using Base::L_;
  using Base::mps_;
  using Base::mpsContainer_;
  using Base::mpoContainer_;
  using Base::mpo_;
  using Base::parms_;
  using Base::runSweepSimulation;
  using Base::siteLeft_;
  using Base::siteRight_;

  /** @brief Class constructor */
  SweepBasedLinearSystem(MPSType& mps, const MPOType& mpo, BaseParameters& parms, int initSite=0)
    : Base(mps, mpo, parms, std::string("Linear system solver"), initSite), adaptiveBondDimension_(false),
      shiftParameter_(0.), isPrecond_(false)
  {
    /* // Folded simulation --> To be reactivated when implementing the folded operator 
    if (parms["pI_folded"] == "yes") {
        maquis::cout << " Activating folded treatment " << std::endl;
        isSquared = true;
    } */
    rhsMps_ = mps;
    overlapPropagator_ = std::make_unique<OverlapPropagatorType>(mps_, rhsMps_, initSite_);
    /* To be reactivated when implementing the folded operator 
    if (isSquared) {
      leftSquared_.resize(mpo.length()+1);
      rightSquared_.resize(mpo.length()+1);
      leftCross_.resize(mpo.length()+1);
      rightCross_.resize(mpo.length()+1);
    } */
    // Adaptive m
    if (parms.is_set("linsystem_truncation_ratio")) {
      adaptiveBondDimension_ = true;
      truncationRatio_ = parms["linsystem_truncation_ratio"].as<double>();
    }
    // Note that we subtract the core energy to the shift parameter (the SiteProblem object
    // does not include that contribution)
    if (parms_.is_set("ipi_shift"))
      shiftParameter_ = parms["ipi_shift"].as<ValueType>()-mpo_.getCoreEnergy();
    if (parms_["linsystem_precond"] == "yes")
      isPrecond_ = true;
  }

  /** @brief Setter for the shift */
  void setShift(ValueType newShift) {
    shiftParameter_ = newShift;
  }

  /** @brief Method called at the beginning of each sweep */
  void prepareSweep() override final {
    iterationResults_.clear();
  }

  /** @brief Method called before each microiteration */
  void prepareMicroiteration() override final {
    siteProblem_ = std::make_unique<SiteProblemType>(boundaryPropagator_->getLeftBoundary(siteLeft_), boundaryPropagator_->getRightBoundary(siteRight_),
                                                     mpoContainer_.getMPOTensor(siteLeft_));
    rhs_ = overlapPropagator_->template getOrthogonalVector<SweepType>(siteLeft_, siteRight_);
    if (isPrecond_)
      preconditioner_ = std::make_unique<BlockMatrixType>(contraction::diagonal_hamiltonian(boundaryPropagator_->getLeftBoundary(siteLeft_),
                                                                                            boundaryPropagator_->getRightBoundary(siteRight_),
                                                                                            mpoContainer_.getMPOTensor(siteLeft_),
                                                                                            mpsContainer_.getMPSTensor(siteLeft_)));
  }

  /** @brief Solution of the site-centered problem */
  MPSTensorType solveLocalProblem() override final {
    auto& mpsToOptimize = mpsContainer_.getMPSTensor(siteLeft_);
    LinearSolverType ls(siteProblem_, mpsToOptimize, rhs_, shiftParameter_, parms_, preconditioner_);
    resultOfLocalSiteProblem_ = ls.res();
    // mps[site] = res.second;
    iterationResults_["Energy"] << resultOfLocalSiteProblem_.first + maquis::real(mpo_.getCoreEnergy());
    return resultOfLocalSiteProblem_.second;
  }

  /** @brief Propagates the boundaries */
  void propagateBoundaries() override final {
    auto sweepType = (indexOfMicroIteration_ < lastSite_) ? SweepDirectionType::Forward : SweepDirectionType::Backward;
    // Boundary propagation
    if (sweepType == SweepDirectionType::Forward) {
      rhsMps_.move_normalization_l2r(siteLeft_, siteLeft_+1, DefaultSolver());
      boundaryPropagator_->updateLeftBoundary(siteLeft_+1);
      if (overlapPropagator_)
        overlapPropagator_->updateLeftOverlapBoundaries(siteLeft_+1);
    }
    else if (sweepType == SweepDirectionType::Backward) {
      rhsMps_.move_normalization_r2l(siteRight_-1, siteRight_-2);
      boundaryPropagator_->updateRightBoundary(siteRight_-1);
      if (overlapPropagator_)
        overlapPropagator_->updateRightOverlapBoundaries(siteRight_-1);
    }
  }

  /** @brief Operations to be executed at the end of a microiteration */
  void finalizeMicroIteration(const truncation_results& trunc) override final {
    iterationResults_["BondDimension"]   << trunc.bond_dimension;
    iterationResults_["TruncatedWeight"] << trunc.truncated_weight;
    iterationResults_["SmallestEV"]      << trunc.smallest_ev;
  }

  /** @brief Operations to be executed at the end of the sweep */
  void finalizeSweep() override final {
    initSite_ = -1;
  }

  /** @brief Whether to normalize the MPS at the end of a half-sweep */
  bool normalizeAtEnd() override final {
    return false;
  }

private:
  // Class members
  MPSType rhsMps_;                                                // RHS for the solution of the linear system.
  std::shared_ptr<BlockMatrixType> preconditioner_;               // If needed, stores the preconditioner.
  bool adaptiveBondDimension_;                                    // Whether to dynamically adapt the bond dimension.
  bool isPrecond_;                                                // If true, activates the preconditioning.
  double truncationRatio_;                                        // Parameter for a DBSS-like solution of the linear system.
  ValueType shiftParameter_;                                      // Shift parameter for the linear system
  MPSTensorType rhs_;                                             // RHS of the local linear system (updated at each microiteration).
  std::unique_ptr<OverlapPropagatorType> overlapPropagator_;      // Object needed to store the partial MPS/MPS contraction
  std::shared_ptr<SiteProblemType> siteProblem_;                  // Site problem associated with the solution of the linear system.
  std::pair<ValueType, MPSTensorType > resultOfLocalSiteProblem_; // TO CHECK IF NEEDED
};

#endif // SWEEP_BASED_LINEAR_SYSTEM_H
