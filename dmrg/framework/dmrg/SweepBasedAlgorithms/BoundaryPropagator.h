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

private:

  /** @brief Generates the left boundary */
  void generateLeftBoundary() {
    Storage::drop(left_[0]);
    left_[0] = mps_.left_boundary();
    Storage::pin(left_[0]);
    for (int i = 0; i < initSite_; i++) {
      Storage::drop(left_[i+1]);
      propagateLeftBoundary(i);
      Storage::evict(left_[i]);
      parallel::sync();
    }
    Storage::evict(left_[initSite_]);
  }

  /** @brief Generates the right boundary */
  void generateRightBoundary() {
    Storage::drop(right_[L_]);
    right_[L_] = mps_.right_boundary();
    Storage::pin(right_[L_]);
    for (int i = L_-1; i >= initSite_; i--) {
      Storage::drop(right_[i]);
      propagateRightBoundary(i);
      Storage::evict(right_[i+1]);
      parallel::sync();
    }
    Storage::evict(right_[initSite_]);
  }

  /** @brief Propagation of the left boundary */
  inline void propagateLeftBoundary(int referenceSite) {
    assert(siteInitial <= siteFinal);
    assert(siteInitial >= 0 && siteInitial <= L_);
    assert(siteFinal >= 0 && siteFinal <= L_);
    left_[referenceSite+1] = Contraction::overlap_mpo_left_step(mps_[referenceSite], mps_[referenceSite],
                                                                left_[referenceSite], mpo_[referenceSite]);
    Storage::pin(left_[referenceSite]);
  }

  /** @brief Propagation of the right boundary */
  inline void propagateRightBoundary(int referenceSite) {
    assert(siteInitial >= siteFinal);
    assert(siteInitial >= 0 && siteInitial <= L_);
    assert(siteFinal >= 0 && siteFinal <= L_);
    right_[referenceSite] = Contraction::overlap_mpo_right_step(mps_[referenceSite], mps_[referenceSite],
                                                                right_[referenceSite+1], mpo_[referenceSite]);
    Storage::pin(right_[referenceSite]);
  }

  // Class members
  int L_;
  BoundariesType left_, right_;
  int initSite_;
  const MPOType& mpo_;
  const MPSType& mps_;
};

#endif 