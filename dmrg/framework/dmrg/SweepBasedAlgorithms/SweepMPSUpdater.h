/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiardi <alberto.baiardi@phys.chem.ethz.ch>
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

#ifndef SWEEP_MPS_UPDATER
#define SWEEP_MPS_UPDATER

#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/ts_ops.h"
#include "BoundaryPropagator.h"
#include "SweepOptimizationTypeTrait.h"

/**
 * @brief Class that manages the update of an MPS at the end of the sweep.
 *
 * Also in this case, the class hides the detail of how the
 *
 * @tparam SweepType can be SingleSite or TwoSite.
 */
template<class Matrix, class SymmGroup, class Storage, SweepOptimizationType SweepType>
class SweepMPSUpdater {};

/** @brief Specialization for the single-site case */
template<class Matrix, class SymmGroup, class Storage>
class SweepMPSUpdater<Matrix, SymmGroup, Storage, SweepOptimizationType::SingleSite> {
public:
  // Type declaration
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using BoundaryPropagatorType = BoundaryPropagator<Matrix, SymmGroup, Storage>;

  /** @brief Class constructor */
  SweepMPSUpdater(const MPOType& mpo, MPSType& mps, std::shared_ptr<BoundaryPropagatorType> boundaryPropagator,
                  BaseParameters& parms)
    : mpo_(mpo), mps_(mps), boundaryPropagator_(boundaryPropagator), parms_(parms)
  {
    L_ = mps_.size();
  };

  /** @brief Method to update the MPS for a given site */
  auto updateMPS(int siteLeft, int siteRight, SweepDirectionType sweepDirection, const MPSTensorType& inputMPS,
                 double alpha, double cutoff, double mMax, bool normalizeEnd) {
    // Printing
    mps_[siteLeft] = inputMPS;
    truncation_results truncationOutput;
    // Forward sweep case
    if (sweepDirection == SweepDirectionType::Forward) {
      if (siteLeft < L_-1) {
        truncationOutput = mps_.grow_l2r_sweep(mpo_[siteLeft], boundaryPropagator_->getLeftBoundary(siteLeft),
                                               boundaryPropagator_->getRightBoundary(siteRight), siteLeft, alpha,
                                               cutoff, mMax);
      }
      else if (normalizeEnd) {
        mps_[siteLeft].leftNormalize(DefaultSolver());
      }
    }
    // Backward case
    else if (sweepDirection == SweepDirectionType::Backward) {
      if (siteLeft > 0) {
        truncationOutput = mps_.grow_r2l_sweep(mpo_[siteLeft], boundaryPropagator_->getLeftBoundary(siteLeft),
                                               boundaryPropagator_->getRightBoundary(siteRight), siteLeft,
                                               alpha, cutoff, mMax);
      }
      else if (normalizeEnd) {
        mps_[siteLeft].rightNormalize(DefaultSolver());
      }
    }
    return truncationOutput;
  }

private:
  std::shared_ptr<BoundaryPropagatorType> boundaryPropagator_;
  const MPOType& mpo_;
  MPSType& mps_;
  BaseParameters& parms_;
  int L_;
};

/** @brief Specialization for the two-site case */
template<class Matrix, class SymmGroup, class Storage>
class SweepMPSUpdater<Matrix, SymmGroup, Storage, SweepOptimizationType::TwoSite> {
public:
  // Type declaration
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using TwoSiteTensorType = TwoSiteTensor<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using BoundaryPropagatorType = BoundaryPropagator<Matrix, SymmGroup, Storage>;

  /** @brief Class constructor */
  SweepMPSUpdater(const MPOType& mpo, MPSType& mps, std::shared_ptr<BoundaryPropagatorType> boundaryPropagator,
                  BaseParameters& parms)
    : mpo_(mpo), mps_(mps), boundaryPropagator_(boundaryPropagator), parms_(parms)
  {
    L_ = mps_.size();
  };

  /** @brief Method to update the MPS for a given site */
  auto updateMPS(int siteLeft, int siteRight, SweepDirectionType sweepDirection, const MPSTensorType& inputMPS,
                 double alpha, double cutoff, double mMax, bool normalizeEnd)
  {
    // Converts back the MPS into the two-site tensor
    TwoSiteTensorType tst(mps_[siteLeft], mps_[siteLeft+1]);
    tst << inputMPS;
    truncation_results truncationOutput;
    // Actual truncation
    if (sweepDirection == SweepDirectionType::Forward) {
      // Write back result from optimization
      if (parms_["twosite_truncation"] == "svd")
        boost::tie(mps_[siteLeft], mps_[siteLeft+1], truncationOutput) = tst.split_mps_l2r(mMax, cutoff);
      else
        boost::tie(mps_[siteLeft], mps_[siteLeft+1], truncationOutput) = tst.predict_split_l2r(mMax, cutoff, alpha, boundaryPropagator_->getLeftBoundary(siteLeft),
                                                                                               mpo_[siteLeft]);
      if (siteRight < L_) {
        auto t = mps_[siteLeft+1].leftNormalizeAndReturn(DefaultSolver());
        mps_[siteRight].multiply_from_left(t);
      }
      else if (normalizeEnd) {
        mps_[siteLeft+1].leftNormalize(DefaultSolver());
      }

    }
    else if (sweepDirection == SweepDirectionType::Backward) {
      if (parms_["twosite_truncation"] == "svd")
        boost::tie(mps_[siteLeft], mps_[siteLeft+1], truncationOutput) = tst.split_mps_r2l(mMax, cutoff);
      else
        boost::tie(mps_[siteLeft], mps_[siteLeft+1], truncationOutput) = tst.predict_split_r2l(mMax, cutoff, alpha, boundaryPropagator_->getRightBoundary(siteRight),
                                                                                               mpo_[siteLeft+1]);
      if (siteLeft > 0) {
        auto t = mps_[siteLeft].rightNormalizeAndReturn(DefaultSolver());
        mps_[siteLeft-1].multiply_from_right(t);
      }
      else if (normalizeEnd) {
        mps_[siteLeft].rightNormalize(DefaultSolver());
      }
    }
    return truncationOutput;
  }

private:
  std::shared_ptr<BoundaryPropagatorType> boundaryPropagator_;
  const MPOType& mpo_;
  MPSType& mps_;
  BaseParameters& parms_;
  int L_;
};

#endif // SWEEP_MPS_UPDATER
