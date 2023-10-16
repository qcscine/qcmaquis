/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SWEEP_MPS_CONTAINER
#define SWEEP_MPS_CONTAINER

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/twositetensor.h"
#include "SweepOptimizationTypeTrait.h"

/**
 * @brief Class wrapping around an MPS, to be used in conjunction with sweep-based algorithms.
 *
 * Similarly to [SweepMPOContainer], this class serves as a wrapper to hide how the MPS
 * for a given site (or, for the TS case, for a pair of neighbouring sites) is constructed.
 *
 * @tparam SweepType can be OneSite or TwoSite.
 */
template<class Matrix, class SymmGroup, SweepOptimizationType SweepType>
class SweepMPSContainer {};

/** @brief Specialization for the single-site case */
template<class Matrix, class SymmGroup>
class SweepMPSContainer<Matrix, SymmGroup, SweepOptimizationType::SingleSite> {
public:
  // Type declaration
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;

  /** @brief Class constructor */
  explicit SweepMPSContainer(const MPSType& mps) : mps_(mps) {};

  /** @brief Const getter for the MPSTensor of a given site */
  const MPSTensorType& getMPSTensor(int siteLeft) const { return mps_[siteLeft]; }

  /** @brief Gets the MPS */
  const auto& getMPS() { return mps_; }

private:
  const MPSType& mps_;
};

/** @brief Specialization for the two-site case */
template<class Matrix, class SymmGroup>
class SweepMPSContainer<Matrix, SymmGroup, SweepOptimizationType::TwoSite> {
private:
  using TwoSiteTensorType = TwoSiteTensor<Matrix, SymmGroup>;

public:
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;

  /** @brief Class constructor */
  explicit SweepMPSContainer(const MPSType& mps) : mps_(mps) { };

  /** @brief Const getter for the MPSTensor */
  const MPSTensorType& getMPSTensor(int siteLeft) const {
    generateTwoSiteTensor(siteLeft);
    return localMPS_;
  }

  /** @brief Gets the MPS */
  const auto& getMPS() { return mps_; }

private:

  /** @brief Generates the two-site tensor */
  void generateTwoSiteTensor(int site) const {
    TwoSiteTensor<Matrix, SymmGroup> tst(mps_[site], mps_[site+1]);
    localMPS_ = tst.make_mps();
  }

  // Class members
  const MPSType& mps_;
  mutable MPSTensorType localMPS_;
};

#endif // SWEEP_MPO_CONTAINER
