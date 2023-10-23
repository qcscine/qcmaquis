/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
  using TimeType = typename base::time_type;
  using base::has_td_part_;
  using base::is_imag_;
  using base::time_step_;
  using base::apply_hamiltonian;

  /** @brief Class constructor */
  RKEvolver(TimeType time_step, bool has_td_part, bool is_imag) : base(time_step, has_td_part, is_imag) {};

  /** @brief Evolution of a MPSTensor */
  void evolve(SiteProblem<Matrix, SymmGroup> const& site_problem, MPSTensor<Matrix, SymmGroup>& matrix,
              bool is_forward, TimeType time_current, TimeType time_step) const override final {
    this->evolve_kernel(site_problem, matrix, is_forward, time_current, time_step);
  }

  /** @brief Evolution of a ZeroSite Tensor */
  void evolve(ZeroSiteProblem<Matrix, SymmGroup> const& site_problem, block_matrix<Matrix, SymmGroup>& matrix,
              bool is_forward, TimeType time_current, TimeType time_step) const override final {
    this->evolve_kernel(site_problem, matrix, is_forward, time_current, time_step);
  }

 private:
  /* Kernel for the time evolution */
  template<class SiteProblem, class MatrixType>
  void evolve_kernel(SiteProblem const& site_problem, MatrixType& matrix, bool is_forward, TimeType time_current,
                     TimeType time_step) const;

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
