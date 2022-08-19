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

#ifndef BOUNDARY_PROPAGATOR_H
#define BOUNDARY_PROPAGATOR_H

#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/utils/parallel/placement.hpp"

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
    parallel::construct_placements(mpo_);
    left_.resize(L_+1);
    right_.resize(L_+1);
    // Generation of the left boundary
    generateLeftBoundary();
    maquis::cout << "Boundaries are partially initialized...\n";
    generateRightBoundary();
    maquis::cout << "Boundaries are fully initialized...\n";
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

  /**
   * @brief Propagation of the left boundary.
   * 
   * Note that this function assumes that the MPS at site [siteInitial]
   * was changed, and that the optimization is now moved up to site [siteFinal].
   * Therefore, all boundaries ranging from left_[siteInitial+1] up to left_[siteFinal]
   * are changed.
   * 
   * @param siteInitial starting site of the propagation.
   * @param siteFInal last site of the boundary propagation.
   */
  inline void propagateLeftBoundary(int siteInitial, int siteFinal) {
    if (siteInitial < siteFinal) {
      for (int iSite = siteInitial; iSite < siteFinal; iSite++) {
        // The incoming boundary can be dropped - it's anyway overwritten
        Storage::drop(left_[iSite+1]);
        left_[iSite+1] = Contraction::overlap_mpo_left_step(mps_[iSite], mps_[iSite],
                                                            left_[iSite], mpo_[iSite]);
        // We start writing the boundary that has been just used
        Storage::evict(left_[iSite]);
        parallel::sync();
      }
    }
  }

  /**
   * @brief Propagation of the right boundary.
   * 
   * Analogously to [propagateLeftBoundary], this function assumes that the MPS at 
   * site [siteInitial] was changed, and that the optimization is now moved to the *right*
   * up to site [siteFinal].
   * Therefore, all boundaries ranging from right_[siteInitial] up to right_[siteFinal+1]
   * are changed.
   * Note that it makes sense to call this function with siteFinal == -1, in order to calculate
   * right_[0] (which should just contain the energy).
   * 
   * @param siteInitial starting site of the propagation.
   * @param siteFInal last site of the boundary propagation.
   */
  inline void propagateRightBoundary(int siteInitial, int siteFinal) {
    if (siteInitial > siteFinal) {
      for (int iSite = siteInitial; iSite > siteFinal; iSite--) {
        Storage::drop(right_[iSite]);
        right_[iSite] = Contraction::overlap_mpo_right_step(mps_[iSite], mps_[iSite],
                                                            right_[iSite+1], mpo_[iSite]);
        Storage::evict(right_[iSite+1]);
        parallel::sync();
      }
    }
  }

private:

  /** @brief Generates the left boundary */
  void generateLeftBoundary() {
    Storage::drop(left_[0]);
    left_[0] = mps_.left_boundary();
    propagateLeftBoundary(0, initSite_);
    Storage::evict(left_[initSite_]);
  }

  /** @brief Generates the right boundary */
  void generateRightBoundary() {
    Storage::drop(right_[L_]);
    right_[L_] = mps_.right_boundary();
    propagateRightBoundary(L_-1, initSite_);
    Storage::evict(right_[initSite_]);
  }

  // Class members
  int L_;
  BoundariesType left_, right_;
  int initSite_;
  const MPOType& mpo_;
  const MPSType& mps_;
};

#endif 