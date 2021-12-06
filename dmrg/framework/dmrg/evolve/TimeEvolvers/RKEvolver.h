/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2021 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef MAQUIS_DMRG_RKEVOLVER_H
#define MAQUIS_DMRG_RKEVOLVER_H

#ifdef DMRG_TD

#include "TimeEvolutionAlgorithm.h"

/**
 * @brief Class implementing a fourth-order RK propagation scheme for MPSs.
 */
template<class Matrix, class SymmGroup>
class RKEvolver : public TimeEvolutionAlgorithm<Matrix, SymmGroup> {
 public:
  /* Types definition */
  using base = TimeEvolutionAlgorithm<Matrix, SymmGroup>;
  using scalar_type = typename MPSTensor<Matrix, SymmGroup>::scalar_type;
  using time_type = typename base::time_type;
  using base::has_td_part_;
  using base::is_imag_;
  using base::time_step_;
  using base::apply_hamiltonian;

  /* Class constructor */
  RKEvolver(time_type time_step, bool has_td_part, bool is_imag) : base(time_step, has_td_part, is_imag) {};

  /* Time-evolution method. Interfaces calling the main kernel routine */
  void evolve(SiteProblem<Matrix, SymmGroup> const& site_problem, MPSTensor<Matrix, SymmGroup>& matrix,
              bool is_forward, time_type time_current, time_type time_step) const
  {
    this->evolve_kernel(site_problem, matrix, is_forward, time_current, time_step);
  };
  void evolve(ZeroSiteProblem<Matrix, SymmGroup> const& site_problem, block_matrix<Matrix, SymmGroup>& matrix,
              bool is_forward, time_type time_current, time_type time_step) const
  {
    this->evolve_kernel(site_problem, matrix, is_forward, time_current, time_step);
  };
 private:

  /* Kernel for the time evolution */
  template<class SiteProblem, class MatrixType>
  void evolve_kernel(SiteProblem const& site_problem, MatrixType& matrix, bool is_forward, time_type time_current,
                     time_type time_step) const;

  /* Routine for the rescaling of the matrix to obtain, in the end, a complex value */
  template< class ArgType, class MatrixType, typename std::enable_if< std::is_same<double, ArgType >::value>::type * = nullptr >
  void final_rescale(MatrixType& input) const { };
  template< class ArgType, class MatrixType, typename std::enable_if< std::is_same< typename std::complex<double>, ArgType >::value>::type * = nullptr >
  void final_rescale(MatrixType& input) const { input *= std::complex<double>(0., 1.) ; };
  template< class ArgType, class MatrixType, typename std::enable_if< std::is_same<double, ArgType >::value>::type * = nullptr >
  void inverse_final_rescale(MatrixType& input) const { };
  template< class ArgType, class MatrixType, typename std::enable_if< std::is_same< typename std::complex<double>, ArgType >::value>::type * = nullptr >
  void inverse_final_rescale(MatrixType& input) const { input *= std::complex<double>(0., -1.) ; };

  /** Routine to add the scaling factor to the RK vectors. */
  template<class MatrixType>
  void FinalRescaling(MatrixType& input_mat, bool is_forward) const;
};

#include "RKEvolver.cpp"

#endif // DMRG_TD

#endif //MAQUIS_DMRG_RKEVOLVER_H
