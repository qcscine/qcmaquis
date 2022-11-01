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

#ifndef SWEEP_MPO_CONTAINER
#define SWEEP_MPO_CONTAINER

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/ts_ops.h"
#include "dmrg/utils/parallel/placement.hpp"
#include "SweepOptimizationTypeTrait.h"

/**
 * @brief Class wrapping around an MPO, to be used in conjunction with sweep-based algorithms.
 * 
 * The scope of this function is to hide the details of how an MPOTensor is constructed
 * in a sweep-based algorithm. For the single-site case, nothing needs to be done but, for 
 * the two-site case, the two-site MPO must be constructed locally. This class enables not
 * to care about this intermediate step
 * 
 * @tparam SweepType can be OneSite or TwoSite.
 */
template<class Matrix, class SymmGroup, SweepOptimizationType SweepType>
class SweepMPOContainer {};

/** @brief Specialization for the single-site case */
template<class Matrix, class SymmGroup>
class SweepMPOContainer<Matrix, SymmGroup, SweepOptimizationType::SingleSite> {
public:
  // Type declaration
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;

  /** @brief Class constructor */
  SweepMPOContainer(const MPOType& mpo, const MPSType& mps) : mpo_(mpo) {};

  /** @brief Getter for the MPOTensor */
  const auto& getMPOTensor(int site) const {
    return mpo_[site];
  }

  /** @brief Getter for the MPO */
  const auto& getMPO() const {
    return mpo_;
  }

private:
  const MPOType& mpo_;
};

/** @brief Specialization for the two-site case */
template<class Matrix, class SymmGroup>
class SweepMPOContainer<Matrix, SymmGroup, SweepOptimizationType::TwoSite> {
public:
  using MPOType = MPO<Matrix, SymmGroup>;
  using MPSType = MPS<Matrix, SymmGroup>;

  /** @brief Class constructor */
  SweepMPOContainer(const MPOType& mpo, const MPSType& mps) : mpo_(mpo) {
    make_ts_cache_mpo(mpo, twoSiteMPOCache_, mps);
  }

  /** @brief Getter for the MPOTensor */
  const auto& getMPOTensor(int site) const {
    return twoSiteMPOCache_[site];
  }

  /** @brief Getter for the MPO */
  const auto& getMPO() const {
    return mpo_;
  }

private:
  const MPOType& mpo_;
  MPOType twoSiteMPOCache_;
};

#endif // SWEEP_MPO_CONTAINER