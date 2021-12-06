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

template<class Matrix, class SymmGroup>
template<class SiteProblem, class MatrixType>
void RKEvolver<Matrix, SymmGroup>::evolve_kernel(SiteProblem const& site_problem, MatrixType& matrix, bool is_forward,
                                                 time_type time_current, time_type time_step) const
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

