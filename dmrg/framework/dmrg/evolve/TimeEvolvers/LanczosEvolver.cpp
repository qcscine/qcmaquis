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

template<class Matrix, class SymmGroup, TimeStepDistributor TimeStepDistributorClass>
template<class SiteProblem, class MatrixType>
void LanczosEvolver<Matrix, SymmGroup, TimeStepDistributorClass>::evolve_kernel(const SiteProblem& site_problem, MatrixType& matrix,
                                                                                bool is_forward, time_type time_current, time_type time_step) const
{
  // Types definition
  //using matrix_type = Eigen::Matrix< typename MatrixType::scalar_type, Eigen::Dynamic, Eigen::Dynamic >;
  //using vector_type = Eigen::Matrix< typename MatrixType::scalar_type, Eigen::Dynamic, 1 >;
  using matrix_type = alps::numeric::matrix< typename MatrixType::scalar_type >;
  using vector_type = std::vector< typename MatrixType::scalar_type >;
  // Initialization.
  MatrixType buffer_vector;
  double error;
  typename MatrixType::real_type norm_local;
  typename MatrixType::scalar_type alpha, beta;
  typename std::vector<MatrixType> lanczos_space;
  // -- Outer loop --
  // This is the loop over the exponential factors. 
  for (int iExp = 0; iExp < numberOfExponentials; iExp++) {
    std::size_t local_dim = 1;
    matrix_type matrix_representation(max_iter_, max_iter_);
    vector_type result_vector(1);
    lanczos_space.resize(0);
    lanczos_space.reserve(max_iter_);
    // First step of the Lanczos iteration
    print_header();
    lanczos_space.push_back(matrix);
    buffer_vector = applyOperator(lanczos_space[0], site_problem, iExp, time_current, is_forward);
    alpha = ietl::dot(matrix, buffer_vector);
    if (is_forward)
      matrix_representation(0,0) = -alpha;
    else
      matrix_representation(0,0) =  alpha;
    buffer_vector -= alpha*lanczos_space[0];
    // +-----------+
    //   MAIN LOOP
    // +-----------+
    for (std::size_t idx = 1; idx < max_iter_; idx++) {
      // -- ACTUAL CALCULATION OF THE EXPONENTIAL --
      apply_exponential(matrix_representation, result_vector, local_dim);
      matrix = result_vector[0]*lanczos_space[0];
      for (std::size_t i = 1; i < local_dim; i++)
        matrix += result_vector[i]*lanczos_space[i];
      if (is_imag_)
        matrix /= ietl::two_norm(matrix);
      // Check if the norm of the new Lanczos vector is non-zero
      norm_local = ietl::two_norm(buffer_vector);
      if (local_dim == 1)
        error = norm_local;
      else
        error = std::norm(result_vector[local_dim-1])*norm_local;
      // Temporary representation of the matrix
      local_dim++;
      print_data(local_dim-1, error);
      if (norm_local < 1.0E-20 || error < threshold_) {
        print_line();
        break;
      } else {
        // Update of the vector space
        buffer_vector /= norm_local;
        lanczos_space.push_back(buffer_vector);
        // Generation of the new vector
        buffer_vector = applyOperator(lanczos_space[idx], site_problem, iExp, time_current, is_forward);
        alpha = ietl::dot(buffer_vector, lanczos_space[idx]);
        buffer_vector -= alpha*lanczos_space[idx] + norm_local*lanczos_space[idx-1];
        //matrix_representation.conservativeResize(idx+1, idx+1);
        if (is_forward) {
          matrix_representation(idx, idx) = -alpha;
          matrix_representation(idx-1, idx) = -norm_local;
          matrix_representation(idx, idx-1) = -norm_local;
        } else {
          matrix_representation(idx, idx) = alpha;
          matrix_representation(idx-1, idx) = norm_local;
          matrix_representation(idx, idx-1) = norm_local;
        }
      }
      if (idx == max_iter_-1)
        print_line();
    }
  }
};

template<class Matrix, class SymmGroup, TimeStepDistributor TimeStepDistributorClass>
template<class SiteProblem, class MatrixType>
MatrixType LanczosEvolver<Matrix, SymmGroup, TimeStepDistributorClass>::applyOperator(const MatrixType& inputVec, const SiteProblem& site_problem,
                                                                                      int idExp, time_type time_current, bool is_forward) const {
  MatrixType outputVec;
  if (has_td_part_) {
    for (int iExp = 0; iExp < numberOfFactorsPerExponential; iExp++) {
      if (iExp == 0) {
          outputVec = factorsAndSteps[idExp][iExp].first*apply_hamiltonian(inputVec, site_problem);
      }
      else {
          outputVec += factorsAndSteps[idExp][iExp].first*apply_hamiltonian(inputVec, site_problem);
      }
      // Adds the TD part
      //outputVec += factorsAndSteps[idExp][iExp].first*site_problem.apply_td_part(inputVec, time_current+time_step_*factorsAndSteps[idExp][iExp].second);
    }
  }
  else {
    outputVec = apply_hamiltonian(inputVec, site_problem);
  }
  return outputVec;
}

template<class Matrix, class SymmGroup, TimeStepDistributor TimeStepDistributorClass>
template<class MatrixType, class VectorType>
void LanczosEvolver<Matrix, SymmGroup, TimeStepDistributorClass>::apply_exponential(MatrixType& hamiltonian_matrix, VectorType& ret, std::size_t local_dim_) const
{
  // Types definition
  //using ArgType = typename MatrixType::Scalar;
  using ArgType = typename MatrixType::value_type;
  matrix_complex H_hess(local_dim_, local_dim_);
  //
  for (std::size_t i = 0; i < local_dim_; i++) {
    for (std::size_t j = 0; j < local_dim_; j++) {
      if (i==j || i==j+1 || i==j-1) {
        //if (!is_imag_) {
        H_hess(i, j) = initial_convert<ArgType>(hamiltonian_matrix(i,j));
        //}
        //else {
        //  H_hess(i, j) = -initial_convert<ArgType>(hamiltonian_matrix(i,j))*std::complex<double>(time_step_, 0.);
        //}
      }
      else {
        H_hess(i, j) = std::complex<double>(0., 0.);
      }
    }
  }
  matrix_complex expM;
  auto coeff = (is_imag_) ? -std::complex<double>(time_step_, 0.) : std::complex<double>(0., time_step_);
  expM = alps::numeric::exp_hermitian(H_hess, coeff, is_imag_);
  /*
  if (is_imag_) {
    // For imaginary-time propagations, we might encounter large numbers, shifts the diagonal.
    Eigen::ComplexEigenSolver<matrix_complex> es;
    es.compute(H_hess);
    matrix_complex exponentialOfEigenValues = matrix_complex(es.eigenvalues().asDiagonal())
      - matrix_complex::Identity(local_dim_, local_dim_)*es.eigenvalues().real().maxCoeff();
    expM = es.eigenvectors()*exponentialOfEigenValues.exp()*es.eigenvectors().transpose();
  }
  else {
    expM = H_hess.exp();
  }
  */
  ret.resize(local_dim_);
  for (std::size_t idx = 0; idx < local_dim_; idx++)
    ret[idx] = final_convert<ArgType>(expM(0, idx));
}