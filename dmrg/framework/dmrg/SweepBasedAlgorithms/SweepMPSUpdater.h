/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SWEEP_MPS_UPDATER
#define SWEEP_MPS_UPDATER

#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/evolve/TimeEvolvers/timeevolver.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/mp_tensors/ts_ops.h"
#include "dmrg/mp_tensors/zerositeproblem.h"
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
                             double alpha, double cutoff, double mMax, bool normalizeEnd, bool perturbDM)
  {
    loadedUnitaryFactor_ = true;
    mps_[siteLeft] = inputMPS;
    truncation_results truncationOutput;
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

#ifdef DMRG_TD
  using TimeEvolverType = TimeEvolver<Matrix, SymmGroup, BaseParameters>;
  /** @brief Method to perform the back-propagation step */
  void performBackPropagation(GrowBoundaryModality boundaryModality, int siteLeft, int siteRight, std::shared_ptr<TimeEvolverType> timeEvolver) {
    int site = siteLeft;
    if (loadedUnitaryFactor_) {
      if (boundaryModality == GrowBoundaryModality::LeftToRight) {
        if (site < L_-1) {
          auto zsp = ZeroSiteProblemType(mpo_[site], mpo_[site+1], boundaryPropagator_->getLeftBoundary(site+1),
                                         boundaryPropagator_->getRightBoundary(site+1));
          timeEvolver->evolve(zsp, zeroSiteTensor_, false, false);
        }
      }
      else if (boundaryModality == GrowBoundaryModality::RightToLeft) {
        if (site > 0) {
          auto zsp = ZeroSiteProblemType(mpo_[site-1], mpo_[site], boundaryPropagator_->getLeftBoundary(site),
                                         boundaryPropagator_->getRightBoundary(site));
          timeEvolver->evolve(zsp, zeroSiteTensor_, false, false);
        }
      }
    }
    else {
      throw std::runtime_error("ERROR: trying to back-propagate before performing the SVD");
    }
  }
#endif // DMRG_TD

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
                             double alpha, double cutoff, double mMax, bool normalizeEnd, bool perturbDM)
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
                                                                                               mpo_[siteLeft], perturbDM);
    }
    else if (boundaryModality == GrowBoundaryModality::RightToLeft) {
      if (parms_["twosite_truncation"] == "svd")
        boost::tie(mps_[siteLeft], mps_[siteLeft+1], truncationOutput) = tst.split_mps_r2l(mMax, cutoff);
      else
        boost::tie(mps_[siteLeft], mps_[siteLeft+1], truncationOutput) = tst.predict_split_r2l(mMax, cutoff, alpha, boundaryPropagator_->getRightBoundary(siteRight),
                                                                                               mpo_[siteLeft+1], perturbDM);
    }
    loadedUnitaryFactor_= true;
    return truncationOutput;
  }

  /** @brief Back-propagates the tensor */
#ifdef DMRG_TD
  using TimeEvolverType = TimeEvolver<Matrix, SymmGroup, BaseParameters>;
  void performBackPropagation(GrowBoundaryModality boundaryModality, int siteLeft, int siteRight, std::shared_ptr<TimeEvolverType> timeEvolver) {
    if (loadedUnitaryFactor_) {
      if (boundaryModality == GrowBoundaryModality::LeftToRight) {
        if (siteRight != L_-1) {
          SiteProblemType sp2(boundaryPropagator_->getLeftBoundary(siteRight-1), boundaryPropagator_->getRightBoundary(siteRight), mpo_[siteRight-1]);
          timeEvolver->evolve(sp2, mps_[siteRight-1], false, false);
        }
      }
      else if (boundaryModality == GrowBoundaryModality::RightToLeft) {
        if (siteLeft != 0) {
          SiteProblemType sp2(boundaryPropagator_->getLeftBoundary(siteLeft), boundaryPropagator_->getRightBoundary(siteLeft+1), mpo_[siteLeft]);
          timeEvolver->evolve(sp2, mps_[siteLeft], false, false);
        }
      }
    }
    else {
      throw std::runtime_error("ERROR: trying to back-propagate before performing SVD");
    }
  }
#endif // DMRG_TD

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
