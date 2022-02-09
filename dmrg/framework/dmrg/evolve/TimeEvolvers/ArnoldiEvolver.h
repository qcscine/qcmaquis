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

#ifndef MAQUIS_DMRG_ARNOLDIEVOLVER_H
#define MAQUIS_DMRG_ARNOLDIEVOLVER_H

#ifdef DMRG_TD

#include <vector>
//#include <Eigen/Core>
//#include <unsupported/Eigen/MatrixFunctions>
#include "TimeEvolutionAlgorithm.h"

template<class Matrix, class SymmGroup>
class ArnoldiEvolver : public TimeEvolutionAlgorithm<Matrix, SymmGroup> {

  /** Types definition */
  using complex_type = std::complex<double>;
  using base = typename TimeEvolutionAlgorithm<Matrix, SymmGroup>::TimeEvolutionAlgorithm;
  using scalar_type = typename MPSTensor<Matrix, SymmGroup>::scalar_type;
  using time_type = typename base::time_type;
  using matrix_complex = alps::numeric::matrix< complex_type >;
  using vector_complex = std::vector< complex_type >;

  /** Types inheritance */
  using base::is_imag_;
  using base::time_step_;
  using base::apply_hamiltonian;
public:
  /* Class constructor */
  ArnoldiEvolver(time_type time_step, bool has_td, bool is_imag, double threshold, std::size_t max_iter) 
    : base(time_step, has_td, is_imag), threshold_(threshold), max_iter_(max_iter) {}

  /* Time evolution method */
  void evolve(SiteProblem<Matrix, SymmGroup> const& site_problem, MPSTensor<Matrix, SymmGroup>& matrix,
              bool is_forward, time_type time_current, time_type time_step) const 
  {
    evolve_kernel(site_problem, matrix, is_forward, time_current, time_step);
  }

  void evolve(ZeroSiteProblem<Matrix, SymmGroup> const& site_problem, block_matrix<Matrix, SymmGroup>& matrix,
              bool is_forward, time_type time_current, time_type time_step) const 
  {
    evolve_kernel(site_problem, matrix, is_forward, time_current, time_step);
  }

 private:

  /* Kernel for the time evolution part */
  template<class SiteProblem, class MatrixType>
  void evolve_kernel(SiteProblem const& site_problem, MatrixType& matrix, bool is_forward, time_type time_current,
                     time_type time_step) const;

  /* Private method interfacing to Eigen matrix exponential calculator */
  template<class MatrixType, class VectorType>
  void apply_exponential(MatrixType& hamiltonian_matrix, VectorType& ret, size_t local_dim_) const;
  template<class SiteProblem, class MatrixType>
  MatrixType applyOperator(const MatrixType& inputVec, const SiteProblem& site_problem, int idExp, time_type time_current, bool is_forward) const;

  /* Real --> Complex conversion routines */
  template< class ArgType, typename std::enable_if< std::is_same<double, ArgType >::value>::type * = nullptr >
  complex_type initial_convert(const ArgType& input) const { return std::complex<double>(input, 0.) ; };
  template< class ArgType, typename std::enable_if< std::is_same< typename std::complex<double>, ArgType >::value>::type * = nullptr >
  complex_type initial_convert(const ArgType& input) const { return input ; };
  template< class ArgType, typename std::enable_if< std::is_same<double, ArgType >::value>::type * = nullptr >
  ArgType final_convert(const complex_type& input) const { return std::real(input) ; };
  template< class ArgType, typename std::enable_if< std::is_same< typename std::complex<double>, ArgType >::value>::type * = nullptr >
  ArgType final_convert(const complex_type& input) const { return input ; };

  /* Methods to print the results of the Arnoldi algorithm. */
  void print_header() const {
    print_line() ;
    std::cout << "  ITERATION  |   ERROR " << std::endl ;
    print_line() ;
  }
  void print_line() const {
    std::cout << "+------------+------------+" << std::endl ;
  }
  void print_data(std::size_t n_iter, time_type error) const {
    char buf[100] ;
    int n = sprintf(buf, "   %2d        |  %1.4E ", static_cast<int>(n_iter), error) ;
    std::cout << buf << std::endl ;
  };

  /* Class members */
  std::size_t max_iter_;
  double threshold_;
};

#include "ArnoldiEvolver.cpp"

#endif // DMRG_TD

#endif // MAQUIS_DMRG_ARNOLDIEVOLVER_H