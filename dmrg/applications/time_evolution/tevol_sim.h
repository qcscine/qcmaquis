/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_TEVOL_SIM_H
#define APP_DMRG_TEVOL_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include "dmrg/sim/sim.h"
#include "dmrg/sim/te_utils.hpp"

template <class Matrix, class SymmGroup, class TimeEvolver>
class tevol_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    typedef typename base::status_type status_type;
    
    using base::mps;
    using base::mpo;
    using base::lat;
    using base::H;
    using base::parms;
    using base::model;
    using base::measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::rfile;
    
public:
    tevol_sim(DmrgParameters const & parms_, ModelParameters const  & model_)
    : base(parms_, model_)
    { }
    
    void run()
    {
        int meas_each    = parms["measure_each"];
        int chkp_each    = parms["chkp_each"];
        int update_each  = parms["update_each"];
        int nsweeps     = parms["nsweeps"];
        int nsweeps_img = parms["nsweeps_img"];
        
        parms = parms.get_at_index("t", init_sweep);
        model = model.get_at_index("t", init_sweep);

        this->model_init();
        this->mps_init();
        
        /// compute nsteps as the min of the three above
        int nsteps = parms["nsweeps"];
        if (meas_each > -1)
            nsteps = std::min(nsteps, meas_each);
        if (chkp_each > -1)
            nsteps = std::min(nsteps, chkp_each);
        if (update_each > -1)
            nsteps = std::min(nsteps, update_each);
        
        #define CHECK_MULTIPLICITY(var)                                               \
        if (var > 0 && var % nsteps != 0)                                             \
            throw std::runtime_error("var must be a multiple of 'nsteps'.");  
        CHECK_MULTIPLICITY(nsweeps)
        CHECK_MULTIPLICITY(nsweeps_img)
        CHECK_MULTIPLICITY(meas_each)
        CHECK_MULTIPLICITY(chkp_each)
        CHECK_MULTIPLICITY(update_each)
        #undef CHECK_MULTIPLICITY
        
        TimeEvolver evolver(&parms, &mps, lat, &H, init_sweep);
        
        int n = nsweeps / nsteps;
        for (int i=0; i < n; ++i) {
            // TODO: introduce some timings
            
            int sweep = i*nsteps;
            if (update_each > -1 && (sweep % update_each) == 0)
            {
                int pc = 0, mc = 0;
                parms = parms.get_at_index("t", sweep, &pc);
                model = model.get_at_index("t", sweep, &mc);
                if (mc > 0 || pc > 0) {
                    this->model_init(sweep);
                    evolver = TimeEvolver(&parms, &mps, lat, &H, sweep);
                }
            } else if (sweep == nsweeps_img) {
                    // since this is just a change in the time step, there is
                    // no need to split the hamiltonian in non-overlapping terms.
                    evolver.prepare_te_terms();
            }
            
            /// time evolution
            evolver(nsteps);
            sweep = evolver.sweep();
            
            /// measurements
            if ((sweep+1) % meas_each == 0 || (sweep+1) == parms["nsweeps"])
            {
                /// measure energy
                double energy = maquis::real(expval(mps, mpo));
                maquis::cout << "Energy " << energy << std::endl;
                
                /// measure observables specified in 'always_measure'
                if (!parms["always_measure"].empty())
                    this->measure(this->results_archive_path(sweep) + "/results/", measurements.sublist(parms["always_measure"]));

                /// write iteration results
                {
                    storage::archive ar(rfile, "w");
                    ar[this->results_archive_path(sweep) + "/parameters"] << parms;
                    ar[this->results_archive_path(sweep) + "/parameters"] << model;
                    ar[this->results_archive_path(sweep) + "/results"] << evolver.iteration_results();
                    ar[this->results_archive_path(sweep) + "/results/Energy/mean/value"] << std::vector<double>(1, energy);
                    // ar[this->results_archive_path(sweep) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
                }
            }
            
            /// write checkpoint
            bool stopped = stop_callback();
            if (stopped || (sweep+1) % chkp_each == 0 || (sweep+1) == parms["nsweeps"])
                checkpoint_state(mps, sweep);
            
            if (stopped) break;
        }
    }
    
private:
    std::string results_archive_path(int sweep) const
    {
        status_type status;
        status["sweep"] = sweep;
        return base::results_archive_path(status);
    }
    
    void checkpoint_state(MPS<Matrix, SymmGroup> const& state, int sweep)
    {
        status_type status;
        status["sweep"] = sweep;
        status["site"]  = -1;
        return base::checkpoint_state(state, status);
    }
};

#endif
