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

#include "ArnoldiEvolver.h"
#include "alps/numeric/matrix/matrix.hpp"
#include <complex>
#include <cstddef>
#include <vector>

template<class Matrix, class SymmGroup>
template<class SiteProblem, class MatrixType>
void ArnoldiEvolver<Matrix, SymmGroup>::evolve_kernel(const SiteProblem& site_problem, MatrixType& matrix,
                                                      bool is_forward, time_type time_current, time_type time_step) const
{
  // Types definition
  //using matrix_type = Eigen::Matrix< typename MatrixType::scalar_type, Eigen::Dynamic, Eigen::Dynamic >;
  //using vector_type = Eigen::Matrix< typename MatrixType::scalar_type, Eigen::Dynamic, 1 >;
  using matrix_type = alps::numeric::matrix< typename MatrixType::scalar_type >;
  using vector_type = std::vector< typename MatrixType::scalar_type >;
  // Initialization.
  typename MatrixType::real_type norm_local;
  typename MatrixType::scalar_type alpha, beta;
  typename std::vector<MatrixType> lanczos_space;
  std::size_t local_dim = 1;
  MatrixType buffer_vector;
  double error;
  matrix_type matrix_representation(max_iter_, max_iter_, 0.);
  vector_type result_vector(1);
  // First step of the Arnoldi iteration
  if (verbose_) {
    print_header();
  }
  lanczos_space.reserve(max_iter_);
  lanczos_space.push_back(matrix/ietl::two_norm(matrix));
  // ==  MAIN LOOP ==
  for (int idx = 0; idx < max_iter_; idx++) {
    // Generation of the new vector
    buffer_vector = apply_hamiltonian(lanczos_space[idx], site_problem);
    for (int idx2 = 0; idx2 < idx+1; idx2++) {
      matrix_representation(idx2, idx) = ietl::dot(lanczos_space[idx2], buffer_vector);
      buffer_vector -= matrix_representation(idx2, idx)*lanczos_space[idx2];
    }
    norm_local = ietl::two_norm(buffer_vector);
    apply_exponential(matrix_representation, result_vector, local_dim);
    matrix = result_vector[0]*lanczos_space[0];
    for (int i = 1; i < local_dim; i++)
      matrix += result_vector[i]*lanczos_space[i];
    if (is_imag_)
      matrix /= ietl::two_norm(matrix);
    // Check if the norm of the new Arnoldi vector is non-zero
    if (local_dim == 1)
      error = norm_local;
    else
      error = std::norm(result_vector[local_dim-1])*norm_local;
    // Temporary representation of the matrix
    local_dim++;
    if (verbose_) {
      print_data(local_dim-1, error);
    }
    if (norm_local < 1.0E-20 || error < threshold_ || idx == max_iter_-1) {
      if (verbose_) {
        print_line();
      }
      break;
    }
    else {
      // Update of the vector space
      matrix_representation(idx+1, idx) = norm_local;
      buffer_vector /= norm_local;
      lanczos_space.push_back(buffer_vector);
    }
  }
};

template<class Matrix, class SymmGroup>
template<class SiteProblem, class MatrixType>
MatrixType ArnoldiEvolver<Matrix, SymmGroup>::applyOperator(const MatrixType& inputVec, const SiteProblem& site_problem,
                                                            int idExp, time_type time_current, bool is_forward) const
{
  MatrixType outputVec;
  outputVec = apply_hamiltonian(inputVec, site_problem);
  return outputVec;
}

template<class Matrix, class SymmGroup>
template<class MatrixType, class VectorType>
void ArnoldiEvolver<Matrix, SymmGroup>::apply_exponential(MatrixType& hamiltonian_matrix, VectorType& ret, std::size_t local_dim_) const
{
  // Types definition
  //using ArgType = typename MatrixType::Scalar;
  using ArgType = typename MatrixType::value_type;
  matrix_complex H_hess(local_dim_, local_dim_);
  //
  for (int i = 0; i < local_dim_; i++) {
    for (int j = 0; j < local_dim_; j++) {
      if (i <= j+1)
        H_hess(i,j) = initial_convert<ArgType>(hamiltonian_matrix(i,j));
      else
        H_hess(i,j) = std::complex<double>(0., 0.);
    }
  }
  //
  matrix_complex expM;
  auto coeff = (is_imag_) ? -std::complex<double>(time_step_, 0.) : std::complex<double>(0., time_step_);
  expM = alps::numeric::exp(H_hess, coeff, is_imag_);
  ret.resize(local_dim_);
  for (std::size_t idx = 0; idx < local_dim_; idx++)
    ret[idx] = final_convert<ArgType>(expM(idx, 0));
}
