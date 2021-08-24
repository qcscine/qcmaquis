/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2019 by Alberto Baiardi <abaiardi@ethz.ch>
 *
 * This software is part of the ALPS libraries, published under the ALPS
 * Library License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Library License along with
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
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

#ifndef MAQUIS_DMRG_TIMEEVOLUTIONALGORITHM_H
#define MAQUIS_DMRG_TIMEEVOLUTIONALGORITHM_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/mp_tensors/zerositeproblem.h"

/**
 * @brief Virtual class representing an algorithm to propagate an MPS
 * 
 * This is the base class for all site-centered time-evolution algorithms.
 * It must implement an evolve method for both SiteProblem and ZeroSiteProblem
 * objects.
 */

template<class Matrix, class SymmGroup>
class TimeEvolutionAlgorithm {
 public:
  /** Types definition */
  using time_type = double;
  using scalar_type = typename MPSTensor<Matrix, SymmGroup>::scalar_type;

  /**
   * @brief Class constructor
   * @param time_step Time-step (in au) for the propagation.
   * @param has_td_part Whether the MPO has a time-dependent part (NYI).
   * @param is_imag true if the propagation is an imaginary-time propagation.
   */
  explicit TimeEvolutionAlgorithm(time_type time_step, bool has_td_part, bool is_imag)
    : time_step_(time_step), has_td_part_(has_td_part), is_imag_(is_imag) {};
  
  /** Virtual destructor */
  virtual ~TimeEvolutionAlgorithm() = default;

  /** Interface to the time-evolution methods */
  virtual void evolve(SiteProblem<Matrix, SymmGroup> const& site_problem, MPSTensor<Matrix, SymmGroup>& matrix,
                      bool is_forward, time_type time_current, time_type time_step) const = 0;
  virtual void evolve(ZeroSiteProblem<Matrix, SymmGroup> const& site_problem, block_matrix<Matrix, SymmGroup>& matrix,
                      bool is_forward, time_type time_current, time_type time_step) const = 0;

  /** Wrapper around the method to propagate a MatrixType (SiteProblem or ZeroSiteProblem) objects */
  template<class MatrixType, class SiteProblem>
  MatrixType apply_hamiltonian(MatrixType const& matrix, SiteProblem const& site_problem) const {
    return site_problem.apply(matrix);
  };
 protected:
  /** Class members */
  bool has_td_part_, is_imag_;
  time_type time_step_;
};

#endif
