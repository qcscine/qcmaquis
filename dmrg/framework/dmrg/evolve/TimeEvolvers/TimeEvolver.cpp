/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#include <memory>
#include "dmrg/evolve/TimeEvolvers/RKEvolver.h"
#include "utils/io.hpp"
#include <ostream>
#include <stdexcept>
#include <string>

#include "dmrg/evolve/TimeEvolvers/TimeEvolver.h"

// -- Empty constructor --
/*
template<class Matrix, class SymmGroup, class ParameterType>
TimeEvolver<Matrix, SymmGroup, ParameterType>::TimeEvolver() :
accuracy_(1.0E-6), is_imag_(false), max_iterations_(100), time_step_(1.)
{
    time_step_ /= (2.418884326505E-2*219474.63E0) ;
    time_step_ /= 2. ;
};
*/

// -- Constructor from parameter object --
template <class Matrix, class SymmGroup, class ParameterType>
TimeEvolver<Matrix, SymmGroup, ParameterType>::TimeEvolver(ParameterType& parms)
    : accuracy_(parms["propagator_accuracy"]),
      is_imag_(parms["imaginary_time"] == "yes"),
      has_td_part_(false),
      max_iterations_(parms["propagator_maxiter"]),
      time_step_(parms["time_step"]),
      time_current_(0.) {
  // Prints header
  // This 2. factor comes from the second-order Trotter formula
  time_step_ /= 2.;
  // Set the time step. Check also all relevant units conversion.
  if (parms["hamiltonian_units"] == "Hartree") {
    if (parms["time_units"] == "fs") {
      time_step_ *= 41.341374575751;
    } else if (parms["time_units"] == "as") {
      time_step_ *= 0.041341374575751;
    } else {
      throw std::runtime_error("Units for the time variable not yet supported");
    }
  } else if (parms["hamiltonian_units"] == "cm-1") {
    if (parms["time_units"] == "fs") {
      time_step_ *= 41.341374575751 / 219474.63;
    } else if (parms["time_units"] == "as") {
      time_step_ *= 0.041341374575751 / 219474.63;
    } else {
      throw std::runtime_error("Units for the time variable not yet supported");
    }
  } else {
    throw std::runtime_error("Units for the Hamiltonian not yet supported");
  }

  maquis::cout << std::endl;
  maquis::cout << "+------------------------------------------+" << std::endl;
  maquis::cout << " Time Evolution" << std::endl;
  maquis::cout << " - Propagation type: ";
  maquis::cout << ((is_imag_) ? "Imaginary-time" : "Real-time");
  maquis::cout << std::endl;
  maquis::cout << " - Time step: " << time_step_ << std::endl;
  maquis::cout << " - Propagation algorithm: ";
  bool verbose = parms["verbose"];
  // Checks if it has a TD part
  if (parms.is_set("TD_perturbation")) {
    has_td_part_ = true;
    std::string intAlgo = parms["TD_integration_algorithm"];
    if (intAlgo == "RungeKutta") {
      time_evolution_algorithm_ =
          std::make_unique<RKEvolver<Matrix, SymmGroup> >(
              time_step_, has_td_part_, is_imag_
          );
      maquis::cout << "Fourth-order Runge Kutta\n";
    } else if (intAlgo == "EMR2") {
      time_evolution_algorithm_ = std::make_unique<LanczosEMR>(
          time_step_, has_td_part_, is_imag_, accuracy_, max_iterations_,
          verbose
      );
      maquis::cout << "Exponential midpoint rule-based Lanczos\n";
    } else if (intAlgo == "CF4") {
      time_evolution_algorithm_ = std::make_unique<LanczosFourthOrder>(
          time_step_, has_td_part_, is_imag_, accuracy_, max_iterations_,
          verbose
      );
      maquis::cout << "Fouth-order Lanczos\n";
    } else {
      throw std::runtime_error("TD integration algorithm not recognized");
    }
  } else {
    if (parms["transcorrelated_hamiltonian"] == true) {
      time_evolution_algorithm_ = std::make_unique<ArnoldiEvolverType>(
          time_step_, has_td_part_, is_imag_, accuracy_, max_iterations_,
          verbose
      );
      maquis::cout << "Arnoldi approximation of the propagator\n";
    } else {
      time_evolution_algorithm_ = std::make_unique<LanczosTI>(
          time_step_, has_td_part_, is_imag_, accuracy_, max_iterations_,
          verbose
      );
      maquis::cout << "Lanczos approximation of the propagator\n";
    }
  }
  maquis::cout << "+------------------------------------------+" << std::endl;
};

template <class Matrix, class SymmGroup, class ParameterType>
typename TimeEvolver<Matrix, SymmGroup, ParameterType>::time_type
TimeEvolver<Matrix, SymmGroup, ParameterType>::get_time() const {
  return this->time_step_;
}

template <class Matrix, class SymmGroup, class ParameterType>
typename TimeEvolver<Matrix, SymmGroup, ParameterType>::time_type
TimeEvolver<Matrix, SymmGroup, ParameterType>::get_time_current() {
  return this->time_current_;
}

template <class Matrix, class SymmGroup, class ParameterType>
void TimeEvolver<Matrix, SymmGroup, ParameterType>::add_to_current_time(
    time_type time_step
) {
  time_current_ += time_step;
}

// +-------+
//  METHODS
// +-------+

template <class Matrix, class SymmGroup, class ParameterType>
template <class SiteProblem, class MatrixType>
void TimeEvolver<Matrix, SymmGroup, ParameterType>::evolve(
    const SiteProblem& siteProblem, MatrixType& matrix, bool isForward,
    bool isTerminal
) const {
  auto actualTimeStep = (isTerminal) ? 2 * time_step_ : time_step_;
  time_evolution_algorithm_->evolve(
      siteProblem, matrix, isForward, time_current_, actualTimeStep
  );
}
