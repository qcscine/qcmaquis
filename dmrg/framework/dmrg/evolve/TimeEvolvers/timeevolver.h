/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
  void evolve(SiteProblem const& site_problem, MatrixType& matrix, bool is_forward, bool isTerminal) const;

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
