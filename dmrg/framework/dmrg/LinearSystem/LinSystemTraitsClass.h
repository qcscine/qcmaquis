/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef LINSYSTEM_TRAIT_CLASS_H
#define LINSYSTEM_TRAIT_CLASS_H

#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"

template<class Matrix, class SymmGroup>
class LinSystemTraitClass {
public:
  using ModelType = Model<Matrix, SymmGroup>;
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using ValueType = typename Matrix::value_type;
  using ChargeType = typename SymmGroup::charge;

  /** @brief Calculates the error associated with the solution of a linear system */
  static double calculateError(const MPSType& mps, const MPSType& rhsMps, const MPOType& mpo, ValueType shift,
                               const ModelType& model, const Lattice& lattice, ChargeType totalQN, int mMax)
  {
    auto traitClass = MPOTimesMPSTraitClass<Matrix, SymmGroup>(mps, model, lattice, totalQN, mMax);
    auto outputMPS = traitClass.applyMPO(mpo);
    double res = 0.;
    // First term, (H - zI)^2. Note that we do not include the core contribution in H since it is
    // not included in the solution of the linear system as well.
    res += norm(outputMPS);
    res -= 2.*maquis::real(shift*overlap(mps, outputMPS));
    res += std::norm(shift)*norm(mps);
    // Second term
    res -= 2.*maquis::real(overlap(rhsMps, outputMPS));
    res += 2.*maquis::real(shift*overlap(rhsMps, mps));
    // Third term
    res += norm(rhsMps);
    return res;
  }
};

/** @brief Specialization for the spin-adpated case */
template<class Matrix>
class LinSystemTraitClass<Matrix, SU2U1> {
public:
  using ModelType = Model<Matrix, SU2U1>;
  using MPSType = MPS<Matrix, SU2U1>;
  using MPOType = MPO<Matrix, SU2U1>;
  using ValueType = typename Matrix::value_type;
  using ChargeType = typename SU2U1::charge;

  /** @brief Calculates the error associated with the solution of a linear system */
  static double calculateError(const MPSType& mps, const MPSType& rhsMps, const MPOType& mpo, ValueType shift,
                               const ModelType& model, const Lattice& lattice, ChargeType totalQN, int mMax)
  {
    throw std::runtime_error("Variance calculation for SU2U1 NYI");
  }
};

template<class Matrix>
class LinSystemTraitClass<Matrix, SU2U1PG> {
public:
  using ModelType = Model<Matrix, SU2U1PG>;
  using MPSType = MPS<Matrix, SU2U1PG>;
  using MPOType = MPO<Matrix, SU2U1PG>;
  using ValueType = typename Matrix::value_type;
  using ChargeType = typename SU2U1PG::charge;

  /** @brief Calculates the error associated with the solution of a linear system */
  static double calculateError(const MPSType& mps, const MPSType& rhsMps, const MPOType& mpo, ValueType shift,
                               const ModelType& model, const Lattice& lattice, ChargeType totalQN, int mMax)
  {
    throw std::runtime_error("Variance calculation for SU2U1PG NYI");
  }
};

#endif