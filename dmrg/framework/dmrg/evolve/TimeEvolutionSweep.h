/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef TIMEEVOLUTIONSWEEP_H
#define TIMEEVOLUTIONSWEEP_H

#ifdef DMRG_TD

#if not defined(WIN32) && not defined(WIN64)
#include <sys/time.h>
#define HAVE_GETTIMEOFDAY
#endif

#include <boost/algorithm/string.hpp>
#include "utils/sizeof.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/results_collector.h"
#include "dmrg/utils/storage.h"
#include "dmrg/utils/time_limit_exception.h"
#include "dmrg/utils/checks.h"
#include "dmrg/evolve/TimeEvolvers/timeevolver.h"
#include "dmrg/evolve/siteshifter.h"
#include "dmrg/evolve/perturber.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/mp_tensors/zerositeproblem.h"
#include "dmrg/mp_tensors/boundary.h"

#include <memory>

/**
 * @brief TimeEvolutionSweep class
 *
 * Class representing an algorithm that performs the time-evolution of an MPS in a sweep=based fashion.
 * From this base virtual class, two classes are derived, namely the single-site and the two-sites
 * evolution.
 */

template<class Matrix, class SymmGroup, class Storage>
class TimeEvolutionSweep
{
  using contr = contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>;
  using boundary = Boundary<typename storage::constrained<Matrix>::type, SymmGroup>;
  using boundaries_type = std::vector< boundary >;
  using PerturberType = Perturber<Matrix, SymmGroup>;
  using scalar_type = typename MPSTensor<Matrix, SymmGroup>::scalar_type;
  using magnitude_type = typename MPSTensor<Matrix, SymmGroup>::magnitude_type;
  using TimeEvolverType = TimeEvolver<Matrix, SymmGroup, BaseParameters>;
  using SiteShifterType = SiteShifter<Matrix, SymmGroup, TimeEvolverType, PerturberType >;
public:
  //! Class constructor
  TimeEvolutionSweep(MPS<Matrix, SymmGroup>& mps , MPO<Matrix, SymmGroup> const & mpo,
                     BaseParameters & parms, boost::function<bool ()> stop_callback_, int site=0)
    : mps_(mps), mpo_(mpo), parms_(parms), stop_callback(stop_callback_),
      do_backpropagation_(true), time_step_(parms["time_step"]),
      isHermitian_(true), initial_site(site)
{
    // Sets overall sizes
    L_ = mps_.length();
    left_.resize(L_+1);
    right_.resize(L_+1);
    // Generate classes needed for propagation
    time_evolver_ = std::make_shared<TimeEvolverType>(parms_);
    time_step_ = time_evolver_->get_time();
    if (parms_["TD_backpropagation"] == "no")
      do_backpropagation_ = false;
    mps_.canonize(0);
    init_left_right(mpo, site);
    maquis::cout << "Done init_left_right" << std::endl;
    // Objects taking care of the movement to the next site of the optimization
    perturber_ = std::make_shared< PerturberType >(left_, right_, mpo_, parms_);
    site_shifter_ = std::make_shared< SiteShifterType >(mps_, time_evolver_, perturber_, site, do_backpropagation_);
  }

  /** Virtual destructor */
  virtual ~TimeEvolutionSweep() = default;

  /**
   * @brief Virtual method that must be overloaded by different implementations of the [TimeEvolutionSweep] algorithm.
   * @param sweep: index of the sweep that is currently done.
   * @param d: direction of the sweep (enum, can be LeftOnly, RightOnly or Both).
   */
  virtual void evolve_sweep(int sweep) = 0;

  /**
   * @brief Method called at the end of each (forward/backward) sweep.
   */
  void performFinalOperations() {
    initial_site = -1;
  }

  /** Getter for the struct containing the results. */
  results_collector const& iteration_results() const {
    return iteration_results_;
  }

  /**
   * @brief Get the energy at the current iteration.
   * @return magnitude_type energy of the MPS
   */
  magnitude_type get_energy() {
    return energy;
  }

protected:
  /**
   * @brief Method that builds the boundaries for the initial MPS.
   * Calls the same methods that are called when the boundaries are moved one step forward/backward.
   * @param mpo: Matrix Product Operator representation of the Hamiltonian.
   * @param site: initial size of the sweep optimization.
   */
  void init_left_right(MPO<Matrix, SymmGroup> const & mpo, int site)
  {
    Storage::drop(left_[0]);
    left_[0] = mps_.left_boundary();
    for (size_t i = 0; i < site; ++i) {
        Storage::drop(left_[i+1]);
        boundary_left_step(mpo, i);
        Storage::StoreToFile(left_[i]);
        Storage::sync();
    }
    Storage::StoreToFile(left_[site]);
    maquis::cout << "Partial initialization of boundaries completed.\n";
    Storage::drop(right_[L_]);
    right_[L_] = mps_.right_boundary();
    for (int i = L_-1; i >= site; --i) {
      Storage::drop(right_[i]);
      boundary_right_step(mpo, i);
      Storage::StoreToFile(right_[i+1]);
      Storage::sync();
    }
    Storage::StoreToFile(right_[site]);
    maquis::cout << "Full initialization of boundaries completed.\n";
  }

  /**
   * @brief Given the MPS and the left boundaries at site i, calculates the new boundaries at site i+1.
   * @param mpo: MPO representation of the operator associated to the boundaries.
   * @param site: site for which the boundaries are shifted.
   */
  inline void boundary_left_step(MPO<Matrix, SymmGroup> const & mpo, int site)
  {
    left_[site+1] = contr::overlap_mpo_left_step(mps_[site], mps_[site], left_[site], mpo_[site]);
  }

  /**
   * @brief Given the MPS and the tight boundaries at site i+1, calculates the new boundaries at site i.
   * @param mpo: MPO representation of the operator associated to the boundaries.
   * @param site: site for which the boundaries are shifted.
   */
  inline void boundary_right_step(MPO<Matrix, SymmGroup> const & mpo, int site)
  {
    right_[site] = contr::overlap_mpo_right_step(mps_[site], mps_[site], right_[site+1], mpo_[site]);
  }

  /**
   * @brief Logarithmic interpolation routine
   * @param y0: initial value for the interpolation.
   * @param y1: final value for the interpolation.
   * @param N: overall number of points.
   * @param i: point for which the interpolated curve is evaluated.
   * @return interpolated value.
   */
  inline double log_interpolate(double y0, double y1, int N, int i) const
  {
    if (N < 2)
      return y1;
    if (y0 == 0)
      return 0;
    double x = log(y1/y0)/(N-1);
    return y0*exp(x*i);
  }

  /**
   * @brief Getter for the cutoff parameter.
   *
   * This method retrieves, for a given sweep, the value of the threshold to be used in the truncation of the MPS.
   * Note that, for the first ngrowsweeps, the threshold is calculated according to a logarithmic decay between
   * the parameters truncation_initial and truncation_final. After that, the threshold is set to truncation_final
   * and kept constant.
   *
   * @param sweep: index of the sweep.
   * @return threshold for the truncation of the MPS.
   */
  double get_cutoff(int sweep) const
  {
    double cutoff;
    if (sweep >= parms_.template get<int>("ngrowsweeps"))
      cutoff = parms_.template get<double>("truncation_final");
    else
      cutoff = this->log_interpolate(parms_.template get<double>("truncation_initial"),
                                     parms_.template get<double>("truncation_final"),
                                     parms_.template get<int>("ngrowsweeps"), sweep);
    std::cout << "Cutoff - " << cutoff << std::endl ;
    return cutoff;
  }

  /**
   * @brief Retrieves the maximum value of the bond dimension that is employed in the MPS truncation.
   * @param sweep: index of the current sweep.
   * @return: the maximum value of the bond dimension.
   */
  std::size_t get_Mmax(int sweep) const
  {
    std::size_t Mmax;
    if (parms_.is_set("sweep_bond_dimensions")) {
      std::vector<std::size_t> ssizes = parms_.template get<std::vector<std::size_t> >("sweep_bond_dimensions");
      if (sweep >= ssizes.size())
          Mmax = *ssizes.rbegin();
      else
          Mmax = ssizes[sweep];
    } else
      Mmax = parms_.template get<std::size_t>("max_bond_dimension");
    return Mmax;
  }

  /* PRIVATE ATTRIBUTES */
  /* Size of the lattice and initial site of the optimization */
  int L_, initial_site;
  /* Struct collecting the result of the sweep */
  results_collector iteration_results_;
  /* Parameter container */
  BaseParameters& parms_;
  magnitude_type energy;
  boost::function<bool ()> stop_callback;
  boundaries_type left_, right_;
  std::shared_ptr< TimeEvolver< Matrix, SymmGroup, BaseParameters > > time_evolver_;
  std::shared_ptr< SiteShifter< Matrix, SymmGroup, TimeEvolver<Matrix, SymmGroup, BaseParameters>, PerturberType > > site_shifter_;
  std::shared_ptr<PerturberType> perturber_;
  // MPSs and MPOs
  MPS<Matrix, SymmGroup>& mps_;
  MPO<Matrix, SymmGroup> const& mpo_;
  // Quantum dynamics simulations
  bool do_backpropagation_, isHermitian_;
  double time_step_;
};

#include "SingleSiteTimeEvolution.hpp"
#include "TwoSiteTimeEvolution.hpp"

#endif // BUILD_DMRG_EVOLVE

#endif
