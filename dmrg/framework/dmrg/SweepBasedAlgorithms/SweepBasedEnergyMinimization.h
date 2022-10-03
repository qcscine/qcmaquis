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

#ifndef SWEEP_BASED_ENERGY_MINIMIZATION_H
#define SWEEP_BASED_ENERGY_MINIMIZATION_H

#include "GenericSweepSimulation.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/optimize/ietl_jacobi_davidson.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/utils/storage.h"
#include "dmrg/utils/time_limit_exception.h"
#include "dmrg/utils/parallel/placement.hpp"
#include "dmrg/utils/checks.h"
#include "BoundaryPropagator.h"
#include "OverlapPropagator.h"

template<class Matrix, class SymmGroup, class Storage, SweepOptimizationType SweepType>
class SweepBasedEnergyMinimization : public GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepType> {
public:
  using Base = GenericSweepSimulation<Matrix, SymmGroup, Storage, SweepType>;
  using OverlapPropagatorType = OverlapPropagator<Matrix, SymmGroup, Storage>;
  using SweepTraitClass = SweepOptimizationTypeTrait<SweepType>;
  using SiteProblemType = SiteProblem<Matrix, SymmGroup>;
  using MPSType = typename Base::MPSType;
  using MPOType = typename Base::MPOType;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
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
  SweepBasedEnergyMinimization(MPSType& mps, const MPOType& mpo, BaseParameters& parms, int initSite=0)
    : Base(mps, mpo, parms, std::string("Optimization"), initSite), nOrtho_(0)
  {
    // mps_.canonize(initSite_);
    if (parms_.is_set("ortho_states") && parms_["ortho_states"] != "") {
      files_ = parms_["ortho_states"].str();
      std::vector<std::string> files;
      boost::split(files, files_, boost::is_any_of(", "));
      if (!parms_.is_set("n_ortho_states"))
        throw std::runtime_error("Please set [n_ortho_states]");
      else
        nOrtho_ = parms_["n_ortho_states"];
      overlapPropagator_ = std::make_unique<OverlapPropagatorType>(mps_, files, parms_);
      if (nOrtho_ != overlapPropagator_->getNumberOfOverlapMPSs())
        throw std::runtime_error("Nuber of chkp files not coherent with [n_ortho_states] parameter");
      orthoLocal_.resize(nOrtho_);
      maquis::cout << "Running a constrained optimization with respect to " << nOrtho_ << " states." << std::endl;
    }
  }

  /** @brief Method called at the beginning of each sweep */
  void prepareSweep() override final {
    iterationResults_.clear();
  }

  /** @brief Method called before each microiteration */
  void prepareMicroiteration() override final {
    siteProblem_ = std::make_unique<SiteProblemType>(boundaryPropagator_->getLeftBoundary(siteLeft_), boundaryPropagator_->getRightBoundary(siteRight_),
                                                     mpoContainer_.getMPOTensor(siteLeft_));
    if (overlapPropagator_)
      for (int iState = 0; iState < nOrtho_; iState++)
        orthoLocal_[iState] = overlapPropagator_->template getOrthogonalVector<SweepType>(iState, siteLeft_, siteRight_);
  }

  /** @brief Solution of the site-centered problem */
  MPSTensorType solveLocalProblem() override final {
    auto& mpsToOptimize = mpsContainer_.getMPSTensor(siteLeft_);
    if (parms_["eigensolver"] == std::string("IETL"))
      resultOfLocalSiteProblem_ = solve_ietl_lanczos(*(siteProblem_.get()), mpsToOptimize, parms_);
    else if (parms_["eigensolver"] == std::string("IETL_JCD"))
      resultOfLocalSiteProblem_ = solve_ietl_jcd(*(siteProblem_.get()), mpsToOptimize, parms_, orthoLocal_);
    else if (parms_["eigensolver"] == std::string("IETL_DAVIDSON"))
      resultOfLocalSiteProblem_ = solve_ietl_jcd(*(siteProblem_.get()), mpsToOptimize, parms_, orthoLocal_);
    else
      throw std::runtime_error("I don't know this eigensolver.");
    // Loads the final results
    auto energy = resultOfLocalSiteProblem_.first + mpo_.getCoreEnergy();
    maquis::cout << std::setprecision(10) << " Energy = " << energy << std::endl;
    iterationResults_["Energy"] << energy;
    return resultOfLocalSiteProblem_.second;
  }

  /** @brief Propagates the boundaries */
  void propagateBoundaries() override final {
    auto sweepType = (indexOfMicroIteration_ < lastSite_) ? SweepDirectionType::Forward : SweepDirectionType::Backward;
    // Boundary propagation
    if (sweepType == SweepDirectionType::Forward) {
      boundaryPropagator_->updateLeftBoundary(siteLeft_+1);
      if (overlapPropagator_)
        overlapPropagator_->updateLeftOverlapBoundaries(siteLeft_+1);
    }
    else if (sweepType == SweepDirectionType::Backward) {
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
    return true;
  }

private:
  // Class members
  int nOrtho_;
  std::vector<MPSTensorType> orthoLocal_;
  std::unique_ptr<OverlapPropagatorType> overlapPropagator_;
  std::unique_ptr<SiteProblemType> siteProblem_;
  std::string files_;
  std::pair<ValueType, MPSTensorType > resultOfLocalSiteProblem_;
};

#endif // SWEEP_BASED_ENERGY_MINIMIZATION_H
