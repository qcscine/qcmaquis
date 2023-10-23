/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

template<class Matrix, class SymmGroup>
template<class SiteProblem, class MatrixType>
void RKEvolver<Matrix, SymmGroup>::evolve_kernel(SiteProblem const& site_problem, MatrixType& matrix, bool is_forward,
                                                 TimeType time_current, TimeType time_step) const
{
  // Typedef
  typedef typename MatrixType::value_type ArgType;
  // Variable initialization
  MatrixType mps_k1, mps_k2, mps_k3, mps_k4, buffer;
  // First RK vector
  mps_k1 = apply_hamiltonian(matrix, site_problem);
  //if (has_td_part_)
  //  mps_k1 += site_problem.apply_td_part(matrix, time_current);
  mps_k1 *= time_step;
  FinalRescaling(mps_k1, is_forward);
  // Second RK vector
  buffer = matrix + mps_k1/2.;
  mps_k2 = apply_hamiltonian(buffer, site_problem); 
  //if (has_td_part_)
  //  mps_k2 += site_problem.apply_td_part(buffer, time_current+time_step/2.);
  mps_k2 *= time_step;
  FinalRescaling(mps_k2, is_forward);
  // Third RK vector
  buffer = matrix + mps_k2/2.;
  mps_k3 = apply_hamiltonian(buffer, site_problem);
  //if (has_td_part_)
  //  mps_k3 += site_problem.apply_td_part(buffer, time_current+time_step/2.);
  mps_k3 *= time_step;
  FinalRescaling(mps_k3, is_forward);
  // Fourth RK vector
  buffer = matrix + mps_k3;
  mps_k4 = apply_hamiltonian(buffer, site_problem);
  //if (has_td_part_)
  //  mps_k4 += site_problem.apply_td_part(buffer, time_current+time_step);
  mps_k4 *= time_step;
  FinalRescaling(mps_k4, is_forward);
  // Final retrieval of the output
  matrix += (mps_k1 + 2.*mps_k2 + 2.*mps_k3 + mps_k4)/6.;
  if (is_imag_)
    matrix /= ietl::two_norm(matrix);
};

template<class Matrix, class SymmGroup>
template<class MatrixType>
void RKEvolver<Matrix, SymmGroup>::FinalRescaling(MatrixType& input_mat, bool is_forward) const {
  typedef typename MatrixType::value_type ArgType;
  if (!is_imag_) {
    if (!is_forward)
      inverse_final_rescale<ArgType, MatrixType>(input_mat);
    else
      final_rescale<ArgType, MatrixType>(input_mat);
  } else {
    input_mat *= -1.;
  }
}

