/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SINGLESITETIMEEVOLUTION_H
#define SINGLESITETIMEEVOLUTION_H

#include <chrono>

#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/evolve/TimeEvolutionSweep.h"
#include "dmrg/utils/DmrgParameters.h"


/**
 * @brief Class implementing the single-site time-evolution algorithm with a sweep-based Trotter decomposition.
 * @tparam Matrix: type of the matrix underlying the MPSTensor class.
 * @tparam SymmGroup: class representing the group of the molecule.
 * 2tparam Storage: class for the storage.
 */

template<class Matrix, class SymmGroup, class Storage>
class SingleSiteTimeEvolution : public TimeEvolutionSweep<Matrix, SymmGroup, Storage>
{
public:
  // Types definition
  typedef typename Matrix::value_type value_type;
  typedef TimeEvolutionSweep<Matrix, SymmGroup, Storage> base;
  // Members and methods inherited from the base class.
  using base::energy;
  using base::initial_site;
  using base::isHermitian_;
  using base::iteration_results_;
  using base::L_;
  using base::left_;
  using base::mpo_;
  using base::mps_;
  using base::parms_;
  using base::performFinalOperations;
  using base::perturber_;
  using base::right_;
  using base::site_shifter_;
  using base::stop_callback;
  using base::time_evolver_;
  using base::time_step_;

  /**
   * @brief Constructor of the SingleSiteTimeEvolution class.
   * @param mps: MPS to be propagated.
   * @param mpo: MPO representation of the Hamiltonian
   * @param tdmpo: (optional) MPO representation of the time-dependent part of the MPO.
   * @param has_td_perturb: if true, a time-dependent perturbation is present.
   * @param parms: parameter container.
   * @param stop_callback: stopping criteria.
   * @param initial_site_: site in which the optimization is started.
   */
  SingleSiteTimeEvolution(MPS<Matrix, SymmGroup> & mps, MPO<Matrix, SymmGroup> const & mpo,
                          BaseParameters & parms, boost::function<bool ()> stop_callback, int initial_site_ = 0)
        : base(mps, mpo, parms, stop_callback, to_site(mps.length(), initial_site_)) { };

  /**
   * @brief Convert a number in [0,2L] to the actual site in [0,L] by mapping [L+1,2L] to [0,L] according to a backward sweep.
   * @param L: overall size of the chain.
   * @param i: index to be converted.
   * @return the converted index.
   */
  inline int to_site(const int L, const int i) const
  {
      if (i < 0) return 0;
      return (i < L) ? i : 2*L - 1 - i;
  }

  /**
   * @brief Routine to perform a sweep (that can be forward, backward or both of them.
   * @param sweep: index of the current sweep.
   * @param d: enum indicating if the sweep is forward, backward or both of them.
   */
  void evolve_sweep(int sweep) override {
    // Initialization
    std::chrono::high_resolution_clock::time_point sweep_now = std::chrono::high_resolution_clock::now();
    // Clears the vector with the results of the previous iteration
    iteration_results_.clear();
    typename MPSTensor<Matrix, SymmGroup>::scalar_type dipole;
    std::size_t L = mps_.length();
    int _site = 0, site;
    if (initial_site != -1)
      _site = initial_site;
    Storage::prefetch(left_[_site]);
    Storage::prefetch(right_[_site+1]);
    // Main loop
    for (; _site < 2*L_; ++_site) {
      std::chrono::high_resolution_clock::time_point now, then;
      site = to_site(L_, _site);
      int lr = (_site < L) ? +1 : -1;
      print_header(sweep, site, lr);
      if (lr == -1 && site == L - 1) {
        maquis::cout << "Syncing storage" << std::endl;
        Storage::sync();
      }
      // Storage management.
      Storage::fetch(left_[site]);
      Storage::fetch(right_[site+1]);
      if (lr == +1 && site+2 <= L)
        Storage::prefetch(right_[site+2]);
      if (lr == -1 && site > 0)
        Storage::prefetch(left_[site-1]);
      assert( left_[site].reasonable());
      assert( right_[site+1].reasonable());
      std::pair<typename maquis::traits::real_type<value_type>::type, MPSTensor<Matrix, SymmGroup> > res;
      SiteProblem<Matrix, SymmGroup> sp(left_[site], right_[site+1], mpo_[site]);
      // +-----------------+
      //  MAIN PART: sweep.
      // +-----------------+
      std::cout << std::endl;
      std::cout << "Forward propagating the " << site << "-th MPSTensor " << std::endl;
      time_evolver_->evolve(sp, mps_[site], true, false);
      if (site == 0) {
          energy = ietl::get_energy(sp, mps_[site]);
          res.first = energy;
          maquis::cout << " Energy after = " << res.first << std::endl;
          maquis::cout << " Energy with core = " << res.first + mpo_.getCoreEnergy() << std::endl;
      }
      res = std::make_pair(energy, mps_[site]);
      // +---------------------+
      //  Collection of results
      // +---------------------+
      auto prec = maquis::cout.precision();
      maquis::cout.precision(15);
      maquis::cout.precision(prec);
      iteration_results_["Energy"] << res.first + mpo_.getCoreEnergy();
      // Loads the noise parameter
      double alpha, time_step_effective;
      int ngs = parms_.template get<int>("ngrowsweeps"), nms = parms_.template get<int>("nmainsweeps");
      if (sweep < ngs)
        alpha = parms_.template get<double>("alpha_initial");
      else if (sweep < ngs + nms)
        alpha = parms_.template get<double>("alpha_main");
      else
        alpha = parms_.template get<double>("alpha_final");
      //
      if (sweep > ngs && parms_.is_set("time_step_larger"))
        time_step_effective = parms_["time_step_larger"];
      else
        time_step_effective = time_step_;
      perturber_->update_alpha(alpha) ;
      // Update of the MPS
      double cutoff = this->get_cutoff(sweep);
      int Mmax = this->get_Mmax(sweep);
      truncation_results trunc;
      // +---------------------+
      //  TRUNCATION OF THE MPS
      // +---------------------+
      // -- FORWARD SWEEP --
      if (lr == +1) {
        if (site < L - 1) {
          maquis::cout << " Alpha = " << alpha << std::endl;
          trunc = site_shifter_->shiftSiteTD(mpo_, left_, right_, alpha, cutoff, Mmax, Mmax, perturber_->is_perturb_dm());
          Storage::drop(right_[site+1]);
          Storage::StoreToFile(left_[site]);
          site_shifter_->forward_shift();
        } else {
          time_evolver_->add_to_current_time(time_step_effective);
          site_shifter_->revert();
          this->boundary_left_step(mpo_, site);
        }
      // -- BACKWARD SWEEP --
      } else if (lr == -1) {
        if (site > 0) {
          maquis::cout << " Alpha = " << alpha << std::endl;
          trunc = site_shifter_->shiftSiteTD(mpo_, left_, right_, alpha, cutoff, Mmax, Mmax, perturber_->is_perturb_dm());
          Storage::drop(left_[site]);
          Storage::StoreToFile(right_[site+1]);
          site_shifter_->backward_shift();
        } else {
          time_evolver_->add_to_current_time(time_step_effective);
          site_shifter_->revert();
          this->boundary_right_step(mpo_, site);
        }
      }
      // +----------------+
      //  Final operations
      // +----------------+
      iteration_results_["BondDimension"] << trunc.bond_dimension;
      iteration_results_["TruncatedWeight"] << trunc.truncated_weight;
      iteration_results_["SmallestEV"] << trunc.smallest_ev;
      std::chrono::high_resolution_clock::time_point sweep_then = std::chrono::high_resolution_clock::now();
      double elapsed = std::chrono::duration<double>(sweep_then - sweep_now).count();
      maquis::cout << " Sweep has been running for " << elapsed << " seconds. \n" << std::endl;
      parallel::meminfo();
      if (stop_callback())
        throw dmrg::time_limit(sweep, _site + 1);
    }
    performFinalOperations();
  };

private:
  /**
   * @brief Prints the header in the log file for each sweep.
   * @param sweep: index of the sweep.
   * @param site: site of the optimization.
   * @param lr: direction of the optimization
   */
  void print_header(int& sweep, int& site, int& lr) {
    char buffer[40];
    int n, a;
    if (lr == 1) {
        a = 2*sweep+1;
        n = sprintf(buffer, "  Sweep number %3d - site number %3d", a, site);
    } else {
        a = 2*sweep+2;
        n = sprintf(buffer, "  Sweep number %3d - site number %3d", a, site);
    }
    std::cout << " +-----------------------------------+" << std::endl;
    std::cout << buffer << std::endl;
    std::cout << " +-----------------------------------+" << std::endl;
  }
};

#endif
