/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_DMRG_SIM_H
#define APP_DMRG_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/shared_ptr.hpp>

#include "dmrg/sim/sim.h"
#include "dmrg/optimize/optimize.h"


template <class Matrix, class SymmGroup>
class dmrg_sim : public sim<Matrix, SymmGroup> {
    // Types definition
    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    typedef typename base::status_type status_type;
    typedef typename base::measurements_type measurements_type;
    // Inheritance of the attributes
    using base::mps_guess;
    using base::mps_sa;
    using base::mps_partial_overlap;
    using base::mpo;
    using base::mpo_squared;
    using base::parms;
    using base::all_measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::init_site;
    using base::n_states;
    using base::rfile_sa;

public:

    dmrg_sim (DmrgParameters & parms_) :

    base(parms_, parms_["n_states_sa"], fill_checkpoint_result_name(parms_, parms_["n_states_sa"]))

    { }

    void run()
    {
        int meas_each = parms["measure_each"];
        int chkp_each = parms["chkp_each"];
        // MPO creation
        MPO<Matrix, SymmGroup> mpoc = mpo;
        if (parms["use_compressed"])
            mpoc.compress(1e-12);
        //
        // Optimizer initialization
        // ------------------------
        boost::shared_ptr<opt_base_t> optimizer;
        if (parms["optimization"] == "singlesite") {
            optimizer.reset( new ss_optimize<Matrix, SymmGroup, storage::disk>
                            (mps_sa, mpoc, parms, stop_callback, init_site, mps_guess, mps_partial_overlap, mpo_squared) );
        } else if(parms["optimization"] == "twosite") {
            optimizer.reset( new ts_optimize<Matrix, SymmGroup, storage::disk>
                            (mps_sa, mpoc, parms, stop_callback, init_site, mps_guess, mps_partial_overlap, mpo_squared) );
        } else {
            throw std::runtime_error("Don't know this optimizer");
        }
        measurements_type always_measurements = this->iteration_measurements(init_sweep);
        //
        // DMRG Sweep optimization
        // -----------------------
        try {
            for (int sweep=init_sweep; sweep < parms["nsweeps"]; ++sweep) {
                // Do the sweep
                optimizer->sweep(sweep, Both);
                storage::disk::sync();
                // Check convergence and see if he has to write something
                bool converged = false;
                typedef typename maquis::traits::real_type<Matrix>::type real_type;
                // Energies for the convergence check, array of the size [n_states][sweeps]
                std::vector<boost::optional<real_type > > ediff_for_convergence_check(n_states, boost::none);

                for (std::size_t state = 0; state < n_states; state++) {
                    if ((sweep+1) % meas_each == 0 || (sweep+1) == parms["nsweeps"]) {
                        {
                            storage::archive ar(rfile_sa[state], "w");
                            ar[results_archive_path(sweep) + "/parameters"] << parms;
                            ar[results_archive_path(sweep) + "/results"] << optimizer->iteration_results()[state];

                            // stop simulation if an energy threshold has been specified
                            // FIXME: this does not work for complex numbers - stknecht feb 2016
                            int prev_sweep = sweep - meas_each;
                            if (prev_sweep >= 0 && parms["conv_thresh"] > 0.)
                            {
                                std::vector<real_type> energies;

                                ar[results_archive_path(sweep) + "/results/Energy/mean/value"] >> energies;
                                real_type emin = *std::min_element(energies.begin(), energies.end());
                                ar[results_archive_path(prev_sweep) + "/results/Energy/mean/value"] >> energies;
                                real_type emin_prev = *std::min_element(energies.begin(), energies.end());

                                ediff_for_convergence_check[state] = std::abs(emin - emin_prev);


                            }
                        }
                        // measure observables specified in 'always_measure'
                        if (always_measurements.size() > 0)
                            this->measure(this->results_archive_path(sweep) + "/results/", always_measurements);
                    }
                }
                // the sweep converged only if the largest energy difference of all states is below the convergence threshold

                bool first_sweep = std::all_of(ediff_for_convergence_check.begin(), ediff_for_convergence_check.end(), [](boost::optional<real_type> i){ return i == boost::none; });
                if (!first_sweep)
                {
                    real_type e_diff_all = (*std::max_element(ediff_for_convergence_check.begin(), ediff_for_convergence_check.end())).get();
                    if (e_diff_all < parms["conv_thresh"] && sweep > init_sweep)
                        converged = true;
                }
                // write checkpoint
                bool stopped = stop_callback() || converged;
                if (stopped || (sweep+1) % chkp_each == 0 || (sweep+1) == parms["nsweeps"])
                    checkpoint_simulation(mps_sa, sweep, -1);
                // Exit condition
                if (stopped) break;
            }
        } catch (dmrg::time_limit const& e) {
            maquis::cout << e.what() << " checkpointing partial result." << std::endl;
            checkpoint_simulation(mps_sa, e.sweep(), e.site());
            for (std::size_t state = 0; state < n_states; state++)
            {
                storage::archive ar(rfile_sa[state], "w");
                ar[results_archive_path(e.sweep()) + "/parameters"] << parms;
                ar[results_archive_path(e.sweep()) + "/results"] << optimizer->iteration_results()[state];
            }
        }
    }

    ~dmrg_sim()
    {
        storage::disk::sync();
    }

private:
    std::string results_archive_path(int sweep) const
    {
        status_type status;
        status["sweep"] = sweep;
        return base::results_archive_path(status);
    }

    void checkpoint_simulation(std::vector< class MPS<Matrix, SymmGroup> > const& state_vec,
                               int sweep,
                               int site)
    {
        status_type status;
        status["sweep"] = sweep;
        status["site"]  = site;
        return base::checkpoint_simulation(state_vec, status);
    }

    std::pair<std::vector<std::string>, std::vector<std::string> > fill_checkpoint_result_name(BaseParameters& p, int n_states_)
    {

//     Initialise checkpoint names for an SA calculation
//     Construct file names for checkpoints for each state
//     First, check if we have a state name in the provided chkpfile
//     (i.e. if checkpoint and/or resultfile match the pattern "something.X.something.h5" where X is the state #
//     if yes, replace it with the correct state number
//     if not, append an underscore and the state number to the state.
      typename std::stringstream ss;
      std::vector<std::string> chkpfile_sa_, rfile_sa_;

      std::string chkpfile = boost::trim_right_copy_if(p["chkpfile"].str(), boost::is_any_of("/ "));
      std::string rfile = p["resultfile"].str();

      boost::regex chkpfile_expr("(.+?)\\.[0-9]+(.*\\.h5)$");
      boost::smatch what;

      bool match = boost::regex_match(chkpfile,what,chkpfile_expr);

      for (std::size_t i = 0; i < n_states_; i++) {
          ss.str("") ;
          ss << i ;
          std::string replacement = "$1."+ss.str()+"$2";
          std::string str = match ? boost::regex_replace(chkpfile, chkpfile_expr, replacement) : chkpfile + '_' + ss.str();
          chkpfile_sa_.push_back(str) ;
      }

      ss.clear();

      // Initialise result file names for an SA calculation, similarly to the checkpoint name initialisation
      match = boost::regex_match(rfile,what,chkpfile_expr);

      for (std::size_t i = 0; i < n_states_; i++) {
          ss.str("") ;
          ss << i ;
          std::string replacement = "$1."+ss.str()+".h5";
          std::string str = match ? boost::regex_replace(rfile, chkpfile_expr, replacement) : rfile + '_' + ss.str();
          rfile_sa_.push_back(str);
      }

      return std::make_pair(chkpfile_sa_, rfile_sa_);
    }

};

#endif
