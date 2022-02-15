/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 by Alberto Baiardi <abaiardi@ethz.ch>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
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

#ifndef TWOSITETIMEEVOLUTION_H
#define TWOSITETIMEEVOLUTION_H

#include "dmrg/mp_tensors/twositetensor.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include <boost/tuple/tuple.hpp>
#include "dmrg/evolve/TimeEvolutionSweep.h"
#include <boost/tuple/tuple.hpp> // Needed for std::tie

/**
 * @brief Class implementing the two-site time-evolution algorithm with a sweep-based Trotter decomposition.
 * @tparam Matrix: type of the matrix underlying the MPSTensor class.
 * @tparam SymmGroup: class representing the group of the molecule.
 * @tparam Storage: class for the storage.
 */

template<class Matrix, class SymmGroup, class Storage>
class TwoSiteTimeEvolution : public TimeEvolutionSweep<Matrix, SymmGroup, Storage>
{
public:
  // Types definition
  typedef typename Matrix::value_type value_type;
  typedef TimeEvolutionSweep<Matrix, SymmGroup, Storage> base;
  using base::do_backpropagation_;
  using base::energy;
  using base::initial_site;
  using base::iteration_results_;
  using base::L_;
  using base::left_;
  using base::mpo_;
  using base::mps_;
  using base::parms_;
  using base::performFinalOperations;
  using base::right_;
  using base::stop_callback;
  using base::time_evolver_;
  using base::time_step_;

  /**
   * @brief Constructor of the TwoSiteTimeEvolution class.
   *
   * @param mps: MPS to be propagated.
   * @param mpo: MPO representation of the Hamiltonian
   * @param tdmpo: (optional) MPO representation of the time-dependent part of the MPO.
   * @param has_td_perturb: if true, a time-dependent perturbation is present.
   * @param parms: parameter container.
   * @param stop_callback: stopping criteria.
   * @param initial_site_: site in which the optimization is started.
   */
  TwoSiteTimeEvolution(MPS<Matrix, SymmGroup>& mps, MPO<Matrix, SymmGroup> const & mpo, BaseParameters & parms_,
                       boost::function<bool ()> stop_callback_, int initial_site_ = 0)
    : base(mps, mpo, parms_, stop_callback_, to_site(mps.length(), initial_site_))
  {
    parallel::guard::serial guard;
    make_ts_cache_mpo(mpo, ts_cache_mpo, mps);
  }

  /**
   * @brief Convert a number in [0,2L] to the actual site in [0,L] by mapping [L+1,2L] to [0,L] according to a backward sweep.
   *
   * @param L: overall size of the chain.
   * @param i: index to be converted.
   * @return the converted index.
   */
  inline int to_site(const int L, const int i) const
  {
    if (i < 0)
      return 0;
    else
      return (i < L-1) ? i : 2*L - 2 - i;
  }

  /**
   * @brief Routine to perform a sweep (that can be forward, backward or both of them.
   *
   * @param sweep: index of the current sweep.
   * @param d: enum indicating if the sweep is forward, backward or both of them.
   */
  void evolve_sweep(int sweep)
  {
    // Initialization
    typename MPSTensor<Matrix, SymmGroup>::scalar_type dipole;
    boost::chrono::high_resolution_clock::time_point sweep_now = boost::chrono::high_resolution_clock::now();
    iteration_results_.clear();
    // Definition of the initial site
    int _site = 0, site;
    if (initial_site != -1) {
      _site = initial_site;
      site = to_site(L_, _site);
    }
    if (_site < L_-1) {
      Storage::prefetch(left_[site]);
      Storage::prefetch(right_[site+2]);
    } 
    else {
      Storage::prefetch(left_[site-1]);
      Storage::prefetch(right_[site+1]);
    }
    // ==  MAIN LOOP - SWEEP OPTIMIZATION ==
    for (; _site < 2*L_-2; ++_site) {
      int lr, site1, site2;
      site = to_site(L_, _site);
      if (_site < L_-1) {
        lr = 1;
        site1 = site;
        site2 = site+1;
        ts_cache_mpo[site1].placement_l = mpo_[site1].placement_l;
        ts_cache_mpo[site1].placement_r = parallel::get_right_placement(ts_cache_mpo[site1],
                                                                        mpo_[site1].placement_l,
                                                                        mpo_[site2].placement_r);
      } else {
        lr = -1;
        site1 = site-1;
        site2 = site;
        ts_cache_mpo[site1].placement_l = parallel::get_left_placement(ts_cache_mpo[site1],
                                                                       mpo_[site1].placement_l,
                                                                       mpo_[site2].placement_r);
        ts_cache_mpo[site1].placement_r = mpo_[site2].placement_r;
      }
      // Printing
      print_header(sweep, site1, site2, lr);
      // Memory management
      if (_site != L_-1) {
        Storage::fetch(left_[site1]);
        Storage::fetch(right_[site2+1]);
      }
      if (lr == +1) {
        if (site2+2 < right_.size())
          Storage::prefetch(right_[site2+2]);
      } else {
        if (site1 > 0)
          Storage::prefetch(left_[site1-1]);
      }
     	// Create TwoSite objects.
      MPSTensor<Matrix, SymmGroup> twin_mps;
      TwoSiteTensor<Matrix, SymmGroup> two_vec(mps_[site1], mps_[site2]);
      twin_mps = two_vec.make_mps();
      // Creates the SiteProblem
      SiteProblem<Matrix, SymmGroup> sp(left_[site1], right_[site2+1], ts_cache_mpo[site1]);
      // == MAIN PART: performs the sweep ==
      std::pair< typename MPSTensor<Matrix, SymmGroup>::magnitude_type, MPSTensor<Matrix, SymmGroup> > res;
      time_evolver_->evolve(sp, twin_mps, false);
      if (site1 == 0)
        energy = ietl::get_energy(sp, twin_mps);
      res = std::make_pair(energy, twin_mps);
      two_vec << twin_mps;
      // +---------------------+
      //  Collection of results
      // +---------------------+
      twin_mps.clear();
      {
        maquis::cout.precision(15);
        if (site1 == 0)
          maquis::cout << " Energy = " << res.first + mpo_.getCoreEnergy() << std::endl;
        iteration_results_["Energy"] << res.first + mpo_.getCoreEnergy();
      }
      auto prec = maquis::cout.precision();
      maquis::cout.precision(prec);
      // +------------------------------------+
      //  Setting up parameters for truncation
      // +------------------------------------+
      double alpha;
      int ngs = parms_["ngrowsweeps"], nms = parms_["nmainsweeps"];
      if (sweep < ngs)
        alpha = parms_["alpha_initial"];
      else if (sweep < ngs + nms)
        alpha = parms_["alpha_main"];
      else
        alpha = parms_["alpha_final"];
      auto cutoff = this->get_cutoff(sweep);
      auto Mmax = this->get_Mmax(sweep);
      truncation_results trunc;
      // +--------------------------------------+
      //  Truncation of the TwoSiteTensor object
      // +--------------------------------------+
      // -- Forward sweep --
      if (lr == +1) {
        if (parms_["twosite_truncation"] == "svd")
          boost::tie(mps_[site1], mps_[site2], trunc) = two_vec.split_mps_l2r(Mmax, cutoff);
        else
          boost::tie(mps_[site1], mps_[site2], trunc) = two_vec.predict_split_l2r(Mmax, cutoff, alpha, left_[site1], mpo_[site1]);
        mps_[site2] /= ietl::two_norm(mps_[site2]);
        two_vec.clear();
        this->boundary_left_step(mpo_, site1);
        if (site2 != L_-1) {
          if (do_backpropagation_) {
            SiteProblem<Matrix, SymmGroup> sp2(left_[site2], right_[site2+1], mpo_[site2]);
            time_evolver_->evolve(sp2, mps_[site2], true);
          }
        } else {
          time_evolver_->add_to_current_time(time_step_);
        }
        if (site1 != L_-2)
          Storage::drop(right_[site2+1]);
      // -- Backward sweep --
      } else if (lr == -1) {
        if (parms_["twosite_truncation"] == "svd")
          boost::tie(mps_[site1], mps_[site2], trunc) = two_vec.split_mps_r2l(Mmax, cutoff);
        else
          boost::tie(mps_[site1], mps_[site2], trunc) = two_vec.predict_split_r2l(Mmax, cutoff, alpha, right_[site2+1], mpo_[site2]);
        two_vec.clear();
        mps_[site1] /= ietl::two_norm(mps_[site1]);
        this->boundary_right_step(mpo_, site2);
        if (site1 != 0) {
          if (do_backpropagation_) {
            SiteProblem<Matrix, SymmGroup> sp2(left_[site1], right_[site1+1], mpo_[site1]);
            time_evolver_->evolve(sp2, mps_[site1], true);
          }
        } else {
          time_evolver_->add_to_current_time(time_step_);
        }
        if(site1 != 0)
          Storage::drop(left_[site1]);
      }
      iteration_results_["BondDimension"]     << trunc.bond_dimension;
      iteration_results_["TruncatedWeight"]   << trunc.truncated_weight;
      iteration_results_["TruncatedFraction"] << trunc.truncated_fraction;
      iteration_results_["SmallestEV"]        << trunc.smallest_ev;
      iteration_results_["Dipole"]            << dipole;
      parallel::meminfo();
      boost::chrono::high_resolution_clock::time_point sweep_then = boost::chrono::high_resolution_clock::now();
      double elapsed = boost::chrono::duration<double>(sweep_then - sweep_now).count();
      maquis::cout << "Sweep has been running for " << elapsed << " seconds." << std::endl;
      if (stop_callback())
        throw dmrg::time_limit(sweep, _site+1);
    }
    performFinalOperations();
  }

private:

    /** Object where the MPO centered on two consecutive sites is stored */
    MPO<Matrix, SymmGroup> ts_cache_mpo, ts_cache_td_mpo;

    /**
     * @brief Printing of the header for a sweep in a two-sites optimization.
     * @param sweep: number of the sweep.
     * @param site1: left site of the two-sites optimization.
     * @param site2: right site of the two-sites optimization.
     * @param lr: direction of the sweep
     */
    void print_header(int& sweep, int& site1, int& site2, int& lr){
        char buffer[50] ;
        int a = (lr == 1) ? 2*sweep+1 : 2*sweep+2;
        std::cout << " +--------------------------------------------+" << std::endl ;
        std::cout << "  Sweep number " << a << " - site numbers " << site1 << " and " << site2 << std::endl;
        std::cout << " +--------------------------------------------+" << std::endl;
    }
};

#endif