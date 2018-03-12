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
    
    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    typedef typename base::status_type status_type;
    typedef typename base::measurements_type measurements_type;
    
    using base::mps;
    using base::mps_sa;
    using base::mpo;
    using base::parms;
    using base::all_measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::init_site;
    using base::rfile;
    
public:
    
    dmrg_sim (DmrgParameters & parms_)
    : base(parms_)
    { }
    
    void run()
    {
        int meas_each = parms["measure_each"];
        int chkp_each = parms["chkp_each"];
        // MPO creation
        if (parms["MODEL"] == std::string("quantum_chemistry"))
            throw std::runtime_error("quantum chemistry simulations disabled");
        MPO<Matrix, SymmGroup> mpoc = mpo;
        if (parms["use_compressed"])
            mpoc.compress(1e-12);
        //
        // Optimizer initialization
        // ------------------------
        boost::shared_ptr<opt_base_t> optimizer;
        if (parms["optimization"] == "singlesite") {
            optimizer.reset( new ss_optimize<Matrix, SymmGroup, storage::disk>
                            (mps, mps_sa, mpoc, parms, stop_callback, init_site) );
        } else if(parms["optimization"] == "twosite") {
            optimizer.reset( new ts_optimize<Matrix, SymmGroup, storage::disk>
                            (mps, mps_sa, mpoc, parms, stop_callback, init_site) );
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
                if ((sweep+1) % meas_each == 0 || (sweep+1) == parms["nsweeps"]) {
                    {
                        storage::archive ar(rfile, "w");
                        ar[results_archive_path(sweep) + "/parameters"] << parms;
                        ar[results_archive_path(sweep) + "/results"] << optimizer->iteration_results();
                        // ar[results_archive_path(sweep) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);

                        // stop simulation if an energy threshold has been specified
                        int prev_sweep = sweep - meas_each;
                        if (prev_sweep >= 0 && parms["conv_thresh"] > 0.)
                        {
                            typedef typename maquis::traits::real_type<Matrix>::type real_type;
                            std::vector<real_type> energies;

                            ar[results_archive_path(sweep) + "/results/Energy/mean/value"] >> energies;
                            real_type emin = *std::min_element(energies.begin(), energies.end());
                            ar[results_archive_path(prev_sweep) + "/results/Energy/mean/value"] >> energies;
                            real_type emin_prev = *std::min_element(energies.begin(), energies.end());
                            real_type e_diff = std::abs(emin - emin_prev);

                            if (e_diff < parms["conv_thresh"])
                                converged = true;
                        }
                    }
                    // measure observables specified in 'always_measure'
                    if (always_measurements.size() > 0)
                        this->measure(this->results_archive_path(sweep) + "/results/", always_measurements);
                }
                // write checkpoint
                bool stopped = stop_callback() || converged;
                if (stopped || (sweep+1) % chkp_each == 0 || (sweep+1) == parms["nsweeps"])
                    checkpoint_simulation(mps, mps_sa, sweep, -1);
                // Exit condition
                if (stopped) break;
            }
        } catch (dmrg::time_limit const& e) {
            maquis::cout << e.what() << " checkpointing partial result." << std::endl;
            checkpoint_simulation(mps, mps_sa, e.sweep(), e.site());
            
            {
                storage::archive ar(rfile, "w");
                ar[results_archive_path(e.sweep()) + "/parameters"] << parms;
                ar[results_archive_path(e.sweep()) + "/results"] << optimizer->iteration_results();
                // ar[results_archive_path(e.sweep()) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
            }
        }
    }
    
    ~dmrg_sim() { storage::disk::sync() ; }

private:
    std::string results_archive_path(int sweep) const {
        status_type status;
        status["sweep"] = sweep;
        return base::results_archive_path(status);
    }
    
    void checkpoint_simulation(MPS<Matrix, SymmGroup> const& state,
                               std::vector< class MPS<Matrix, SymmGroup> > const& state_vec,
                               int sweep,
                               int site)
    {
        status_type status;
        status["sweep"] = sweep;
        status["site"]  = site;
        return base::checkpoint_simulation(state, state_vec, status);
    }

};

#endif
