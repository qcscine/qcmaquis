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
#include "dmrg/evolve/TimeEvolvers/timeevolver.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/siteproblem.h"
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
  using BlockMatrixType = block_matrix<Matrix, SymmGroup>;
  using BoundaryPropagatorType = BoundaryPropagator<Matrix, SymmGroup, Storage>;
  using Contractor = typename contraction::Engine<Matrix, Matrix, SymmGroup>;
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using TimeEvolverType = TimeEvolver<Matrix, SymmGroup, BaseParameters>;
  using ZeroSiteProblemType = ZeroSiteProblem<Matrix, SymmGroup>;

  /** @brief Class constructor */
  SweepMPSUpdater(const MPOType& mpo, MPSType& mps, std::shared_ptr<BoundaryPropagatorType> boundaryPropagator,
                  BaseParameters& parms, bool verbose)
    : mpo_(mpo), mps_(mps), boundaryPropagator_(boundaryPropagator), parms_(parms), verbose_(verbose),
      loadedUnitaryFactor_(false)
  {
    L_ = mps_.size();
  }

  /** @brief Method to perform the truncated SVD the MPS for a given site */
  auto generateUnitaryFactor(int siteLeft, int siteRight, GrowBoundaryModality boundaryModality, const MPSTensorType& inputMPS,
                             double alpha, double cutoff, double mMax, bool normalizeEnd)
  {
    loadedUnitaryFactor_ = true;
    mps_[siteLeft] = inputMPS;
    truncation_results truncationOutput;
    bool perturbDM = true;
    MPSTensorType unitaryFactor;
    // Forward sweep case
    if (boundaryModality == GrowBoundaryModality::LeftToRight) {
      if (siteLeft < L_-1) {
        /*
         truncationOutput = mps_.grow_l2r_sweep(mpo_[siteLeft], boundaryPropagator_->getLeftBoundary(siteLeft),
                                                boundaryPropagator_->getRightBoundary(siteRight), siteLeft, alpha,
                                                cutoff, mMax, true, verbose_);
        */
        boost::tie(unitaryFactor, truncationOutput) = Contractor::predict_new_state_l2r_sweep(mps_[siteLeft], mpo_[siteLeft], boundaryPropagator_->getLeftBoundary(siteLeft),
                                                                                              boundaryPropagator_->getRightBoundary(siteRight), alpha, cutoff, mMax,
                                                                                              perturbDM, verbose_);
        zeroSiteTensor_ = Contractor::getZeroSiteTensorL2R(mps_[siteLeft+1], mps_[siteLeft], unitaryFactor);
        mps_[siteLeft] = unitaryFactor;
      }
      else if (normalizeEnd) {
        mps_[siteLeft].leftNormalize(DefaultSolver());
      }
    }
    // Backward case
    else if (boundaryModality == GrowBoundaryModality::RightToLeft) {
      if (siteLeft > 0) {
        /*
          truncationOutput = mps_.grow_r2l_sweep(mpo_[siteLeft], boundaryPropagator_->getLeftBoundary(siteLeft),
                                                 boundaryPropagator_->getRightBoundary(siteRight), siteLeft, alpha,
                                                 cutoff, mMax, true, verbose_);
        */
        boost::tie(unitaryFactor, truncationOutput) = Contractor::predict_new_state_r2l_sweep(mps_[siteLeft], mpo_[siteLeft], boundaryPropagator_->getLeftBoundary(siteLeft),
                                                                                              boundaryPropagator_->getRightBoundary(siteRight), alpha, cutoff, mMax,
                                                                                              perturbDM, verbose_);
        zeroSiteTensor_ = Contractor::getZeroSiteTensorR2L(mps_[siteLeft-1], mps_[siteLeft], unitaryFactor);
        mps_[siteLeft] = unitaryFactor;
      }
      else if (normalizeEnd) {
        mps_[siteLeft].rightNormalize(DefaultSolver());
      }
    }
    return truncationOutput;
  }

  /** @brief Method to perform the back-propagation step */
  void performBackPropagation(GrowBoundaryModality boundaryModality, int siteLeft, int siteRight, std::shared_ptr<TimeEvolverType> timeEvolver) {
    int site = siteLeft;
    if (loadedUnitaryFactor_) {
      if (boundaryModality == GrowBoundaryModality::LeftToRight) {
        if (site < L_-1) {
          auto zsp = ZeroSiteProblemType(mpo_[site], mpo_[site+1], boundaryPropagator_.getLeftBoundary(site),
                                         boundaryPropagator_.getRightBoundary(site+1));
          timeEvolver->evolve(zsp, zeroSiteTensor_, true);
        }
      }
      else if (boundaryModality == GrowBoundaryModality::RightToLeft) {
        if (site > 0) {
          auto zsp = ZeroSiteProblemType(mpo_[site-1], mpo_[site], boundaryPropagator_.getLeftBoundary(site),
                                         boundaryPropagator_.getRightBoundary(site));
          timeEvolver->evolve(zsp, zeroSiteTensor_, true);
        }
      }
    }
    else {
      throw std::runtime_error("ERROR: trying to back-propagate before performing the SVD");
    }
  }

  /** @brief Moves the normalization to the previous/next site (depending whether it's l2r or r2l sweep) */
  void mergeUnitaryFactor(GrowBoundaryModality boundaryModality, int siteLeft, int siteRight, bool normalizeEnd) {
    int site = siteLeft;
    if (boundaryModality == GrowBoundaryModality::LeftToRight) {
      if (site < L_-1) {
        // mps_[site+1] = Contractor::predict_lanczos_l2r_sweep(mps_[site+1], mps_[site], unitaryFactor_);
        // mps_[site] = unitaryFactor_;
        mps_[site+1].multiply_from_left(zeroSiteTensor_);
      }
    }
    else if (boundaryModality == GrowBoundaryModality::RightToLeft) {
      if (site > 0) {
        // mps_[site-1] = Contractor::predict_lanczos_r2l_sweep(mps_[site-1], mps_[site], unitaryFactor_);
        // mps_[site] = unitaryFactor_;
        mps_[site-1].multiply_from_right(zeroSiteTensor_);
      }
    }
    loadedUnitaryFactor_ = false;
  }

private:
  std::shared_ptr<BoundaryPropagatorType> boundaryPropagator_;
  const MPOType& mpo_;
  MPSType& mps_;
  BlockMatrixType zeroSiteTensor_;
  BaseParameters& parms_;
  int L_;
  bool verbose_, loadedUnitaryFactor_;
};

/** @brief Specialization for the two-site case */
template<class Matrix, class SymmGroup, class Storage>
class SweepMPSUpdater<Matrix, SymmGroup, Storage, SweepOptimizationType::TwoSite> {
public:
  // Type declaration
  using BoundaryPropagatorType = BoundaryPropagator<Matrix, SymmGroup, Storage>;
  using Contractor = typename contraction::Engine<Matrix, Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using SiteProblemType = SiteProblem<Matrix, SymmGroup>;
  using TimeEvolverType = TimeEvolver<Matrix, SymmGroup, BaseParameters>;
  using TwoSiteTensorType = TwoSiteTensor<Matrix, SymmGroup>;

  /** @brief Class constructor */
  SweepMPSUpdater(const MPOType& mpo, MPSType& mps, std::shared_ptr<BoundaryPropagatorType> boundaryPropagator,
                  BaseParameters& parms, bool verbose)
    : mpo_(mpo), mps_(mps), boundaryPropagator_(boundaryPropagator), parms_(parms), verbose_(verbose),
      loadedUnitaryFactor_(false)
  {
    L_ = mps_.size();
  }

  /** @brief Method to perform the truncated SVD the MPS for a given site */
  auto generateUnitaryFactor(int siteLeft, int siteRight, GrowBoundaryModality boundaryModality, const MPSTensorType& inputMPS,
                             double alpha, double cutoff, double mMax, bool normalizeEnd)
  {
    // Converts back the MPS into the two-site tensor. Note that here the tst is *not* the contraction of
    // mps_[siteLeft] and mps_[siteLeft+1], since tst << inputMPS overwrites this contraction. mps_[siteLeft]
    // and mps_[siteLeft+1] are used only to have the correct indices
    TwoSiteTensorType tst(mps_[siteLeft], mps_[siteLeft+1]);
    tst << inputMPS;
    truncation_results truncationOutput;
    // Actual truncation
    if (boundaryModality == GrowBoundaryModality::LeftToRight) {
      // Write back result from optimization
      if (parms_["twosite_truncation"] == "svd")
        boost::tie(mps_[siteLeft], mps_[siteLeft+1], truncationOutput) = tst.split_mps_l2r(mMax, cutoff);
      else
        boost::tie(mps_[siteLeft], mps_[siteLeft+1], truncationOutput) = tst.predict_split_l2r(mMax, cutoff, alpha, boundaryPropagator_->getLeftBoundary(siteLeft),
                                                                                               mpo_[siteLeft]);
    }
    else if (boundaryModality == GrowBoundaryModality::RightToLeft) {
      if (parms_["twosite_truncation"] == "svd")
        boost::tie(mps_[siteLeft], mps_[siteLeft+1], truncationOutput) = tst.split_mps_r2l(mMax, cutoff);
      else
        boost::tie(mps_[siteLeft], mps_[siteLeft+1], truncationOutput) = tst.predict_split_r2l(mMax, cutoff, alpha, boundaryPropagator_->getRightBoundary(siteRight),
                                                                                               mpo_[siteLeft+1]);
    }
    loadedUnitaryFactor_= true;
    return truncationOutput;
  }

  /** @brief Back-propagates the tensor */
  void performBackPropagation(GrowBoundaryModality boundaryModality, int siteLeft, int siteRight, std::shared_ptr<TimeEvolverType> timeEvolver) {
    if (loadedUnitaryFactor_) {
      if (siteRight != L_-1) {
        SiteProblemType sp2(boundaryPropagator_->getLeftBoundary(siteRight), boundaryPropagator_->getLeftBoundary(siteRight+1), mpo_[siteRight]);
        timeEvolver->evolve(sp2, mps_[siteRight], true);
      }
    }
    else {
      throw std::runtime_error("ERROR: trying to back-propagate before performing SVD");
    }
  }

  /**
   * @brief Final merging of the unitary factor.
   * Remember that, in the two-site case, siteRight = siteLeft+2.
   */
  void mergeUnitaryFactor(GrowBoundaryModality boundaryModality, int siteLeft, int siteRight, bool normalizeEnd) {
    if (boundaryModality == GrowBoundaryModality::LeftToRight) {
      // TODO Check if this is really needed
      if (siteRight < L_) {
        auto t = mps_[siteLeft+1].leftNormalizeAndReturn(DefaultSolver());
        mps_[siteRight].multiply_from_left(t);
      }
      else if (normalizeEnd) {
        mps_[siteLeft+1].leftNormalize(DefaultSolver());
      }
    }
    else if (boundaryModality == GrowBoundaryModality::RightToLeft) {
      // TODO Check if this is really needed
      if (siteLeft > 0) {
        auto t = mps_[siteLeft].rightNormalizeAndReturn(DefaultSolver());
        mps_[siteLeft-1].multiply_from_right(t);
      }
      else if (normalizeEnd) {
        mps_[siteLeft].rightNormalize(DefaultSolver());
      }
    }
  }

private:
  std::shared_ptr<BoundaryPropagatorType> boundaryPropagator_;
  const MPOType& mpo_;
  MPSType& mps_;
  MPSTensorType unitaryFactor_;
  BaseParameters& parms_;
  int L_;
  bool verbose_, loadedUnitaryFactor_;
};

#endif // SWEEP_MPS_UPDATER
