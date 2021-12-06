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

#ifndef TIME_EVOLVER_H
#define TIME_EVOLVER_H

#ifdef DMRG_TD

#include <sstream>
#include <algorithm>
#include <numeric>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "TimeEvolutionAlgorithm.h"
#include "LanczosEvolver.h"
#include "RKEvolver.h"

/** @brief Wrapper class for time-evolution methods */
template<class Matrix, class SymmGroup, class ParameterType>
class TimeEvolver {
 public:

  /* Types definition */
  using dim_type = std::size_t;
  using scalar_type = typename MPSTensor<Matrix, SymmGroup>::scalar_type;
  using time_evolution_algorithm = TimeEvolutionAlgorithm<Matrix, SymmGroup>;
  using time_type = double;
  using LanczosTI = LanczosEvolver<Matrix, SymmGroup, TimeStepDistributor::Standard>;
  using LanczosEMR = LanczosEvolver<Matrix, SymmGroup, TimeStepDistributor::ExponentialMidpoint>;
  using LanczosFourthOrder = LanczosEvolver<Matrix, SymmGroup, TimeStepDistributor::FourthOrderMagnus>;

  /* Class constructors */
  //TimeEvolver();
  explicit TimeEvolver(ParameterType& parms);
  /* Getters */
  bool isImag() const { return is_imag_; };
  time_type get_time() const;
  time_type get_time_current();
  void add_to_current_time(time_type time_step);

  /* Time-evolution methods */
  template<class SiteProblem, class MatrixType>
  void evolve(SiteProblem const& site_problem, MatrixType& matrix, bool is_forward) const;

 private:
  /* Private attributes */
  bool is_imag_, has_td_part_;
  dim_type max_iterations_;
  time_type accuracy_, time_step_, time_current_;
  std::unique_ptr<time_evolution_algorithm> time_evolution_algorithm_;
};

#include "timeevolver.cpp"

#endif // DMRG_TD

#endif
