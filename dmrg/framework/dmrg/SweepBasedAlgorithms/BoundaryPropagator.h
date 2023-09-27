/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef BOUNDARY_PROPAGATOR_H
#define BOUNDARY_PROPAGATOR_H

#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"

/**
 * @brief This class serves as a wrapper around the boundary propagation routine.
 *
 * The boundary propagation is represented by the following contraction patters:
 *
 *         o--o--o--o--o--o
 *         |  |  |  |  |  |
 *         +--+--+--+--+--+
 *         |  |  |  |  |  |
 *         o--o--o--o--o--o
 *
 * And the result of the contraction between the MPS and the MPO are stored in
 * so-called boundaries.
 * The (i)-th left boundary collects the partial MPS/MPO contraction up to site (i)
 * from the left, with the (i)-th site *excluded*.
 * The (i)-th right boundary collects instead the partial MPS/MPO contraction up
 * to site (i) included.
 * When solving the local problem on a given site j one, therefore, needs, the j-th
 * left boundaries and the (j+1)-th right boundary (for the single-site case).
 */
template<class Matrix, class SymmGroup, class Storage>
class BoundaryPropagator {
public:
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using BoundariesType = std::vector<Boundary<typename storage::constrained<Matrix>::type, SymmGroup> >;
  using Contraction = contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>;

  /** @brief Class constructor */
  BoundaryPropagator(const MPSType& mps, const MPOType& mpo, int initSite=0)
    : initSite_(initSite), L_(mps.length()), mps_(mps), mpo_(mpo)
  {
    left_.resize(L_+1);
    right_.resize(L_+1);
    generateLeftBoundary();
    generateRightBoundary();
  }

  /** @brief Getter for the left boundary */
  auto& getLeftBoundary(int iSite) {
    assert(iSite >= 0 && iSite <= L_);
    return left_[iSite];
  }

  /** @brief Const getter for the left boundary */
  const auto& getLeftBoundary(int iSite) const {
    assert(iSite >= 0 && iSite <= L_);
    return left_[iSite];
  }

  /** @brief Getter for the right boundary */
  auto& getRightBoundary(int iSite) {
    assert(iSite >= 0 && iSite <= L_);
    return right_[iSite];
  }

  /** @brief Const getter for the right boundary */
  const auto& getRightBoundary(int iSite) const {
    assert(iSite >= 0 && iSite <= L_);
    return right_[iSite];
  }

  /** @brief Getter for the left boundary vector */
  auto& getLeftBoundaries() { return left_; }

  /** @brief Const getter for the left boundary vector */
  const auto& getLeftBoundaries() const { return left_; }

  /** @brief Getter for the right boundary */
  auto& getRightBoundaries() { return right_; }

  /** @brief Const getter for the right boundary */
  const auto& getRightBoundary() const { return right_; }

  /**
   * @brief Propagation of the left boundary.
   *
   * Updates the left boundary element that is sitting on [iSite].
   * Therefore, it contracts left[iSite-1] with mps[iSite-1].
   *
   * @param iSite site on which the left boundary must be updated.
   */
  inline void updateLeftBoundary(int iSite) {
    if (iSite > 0 && iSite <= L_) {
      Storage::drop(left_[iSite]);
      left_[iSite] = Contraction::overlap_mpo_left_step(mps_[iSite-1], mps_[iSite-1],
                                                        left_[iSite-1], mpo_[iSite-1]);
      Storage::StoreToFile(left_[iSite-1]);
    }
  }

  /**
   * @brief Propagation of the right boundary.
   *
   * Updates the right boundary element that is sitting on [iSite].
   * Therefore, contracts right[iSite+1] with mps_[iSite] to yield right[iSite].
   *
   * @param iSite site on which the right boundary must be updated.
   */
  inline void updateRightBoundary(int iSite) {
    if (iSite >= 0 && iSite < L_) {
      Storage::drop(right_[iSite]);
      right_[iSite] = Contraction::overlap_mpo_right_step(mps_[iSite], mps_[iSite],
                                                          right_[iSite+1], mpo_[iSite]);
      Storage::StoreToFile(right_[iSite+1]);
    }
  }

private:

  /** @brief Generates the left boundary */
  void generateLeftBoundary() {
    Storage::drop(left_[0]);
    left_[0] = mps_.left_boundary();
    for (int iSite = 1; iSite <= initSite_; iSite++)
      updateLeftBoundary(iSite);
    Storage::StoreToFile(left_[initSite_]);
  }

  /** @brief Generates the right boundary */
  void generateRightBoundary() {
    Storage::drop(right_[L_]);
    right_[L_] = mps_.right_boundary();
    for (int iSite = L_-1; iSite > initSite_; iSite--)
      updateRightBoundary(iSite);
    Storage::StoreToFile(right_[initSite_+1]);
  }

  // Class members
  int L_;
  BoundariesType left_, right_;
  int initSite_;
  const MPOType& mpo_;
  const MPSType& mps_;
};

#endif
