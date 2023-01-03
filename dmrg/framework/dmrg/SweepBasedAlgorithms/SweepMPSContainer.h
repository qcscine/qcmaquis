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
