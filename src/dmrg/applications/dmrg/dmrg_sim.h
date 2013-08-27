/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
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
#include "dmrg/mp_tensors/optimize.h"


template <class Matrix, class SymmGroup>
class dmrg_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    typedef typename base::status_type status_type;

    using base::mps;
    using base::mpo;
    using base::mpoc;
    using base::parms;
    using base::model;
    using base::measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::init_site;
    using base::rfile;
    
public:
    
    dmrg_sim (DmrgParameters & parms_, ModelParameters & model_)
    : base(parms_, model_)
    { }
    
    void run()
    {
        int meas_each = parms["measure_each"];
        int chkp_each = parms["chkp_each"];
        
        this->model_init();
        this->mps_init();
        
        maquis::cout("MPS initialization has finished..."); // MPS restored now

        boost::shared_ptr<opt_base_t> optimizer;
        if (parms["optimization"] == "singlesite")
        {
            optimizer.reset( new ss_optimize<Matrix, SymmGroup, storage::disk>
                            (mps, mpoc, parms, stop_callback, init_sweep, init_site) );
        }
        
        else if(parms["optimization"] == "twosite")
        {
            optimizer.reset( new ts_optimize<Matrix, SymmGroup, storage::disk>
                            (mps, mpoc, parms, stop_callback, init_sweep, init_site) );
        }
        else {
            throw std::runtime_error("Don't know this optimizer");
        }
        
        try {
            for (int sweep=init_sweep; sweep < parms["nsweeps"]; ++sweep) {
                // TODO: introduce some timings
                
                optimizer->sweep(sweep, Both);
                storage::disk::sync();
                
                if ((sweep+1) % meas_each == 0 || (sweep+1) == parms["nsweeps"])
                {
                    /// write iteration results
                    {
                        storage::archive ar(rfile, "w");
                        ar[results_archive_path(sweep) + "/parameters"] << parms;
                        ar[results_archive_path(sweep) + "/parameters"] << model;
                        ar[results_archive_path(sweep) + "/results"] << optimizer->iteration_results();
                        // ar[results_archive_path(sweep) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
                    }
                    
                    /// measure observables specified in 'always_measure'
                    if (!parms["always_measure"].empty())
                        this->measure(this->results_archive_path(sweep) + "/results/", measurements.sublist(parms["always_measure"]));
                }
                
                /// write checkpoint
                bool stopped = stop_callback();
                if (stopped || (sweep+1) % chkp_each == 0 || (sweep+1) == parms["nsweeps"])
                    checkpoint_state(mps, sweep, -1);
                
                if (stopped) break;
            }
        } catch (dmrg::time_limit const& e) {
            maquis::cout << e.what() << " checkpointing partial result." << std::endl;
            checkpoint_state(mps, e.sweep(), e.site());
            
            {
                storage::archive ar(rfile, "w");
                ar[results_archive_path(e.sweep()) + "/parameters"] << parms;
                ar[results_archive_path(e.sweep()) + "/parameters"] << model;
                ar[results_archive_path(e.sweep()) + "/results"] << optimizer->iteration_results();
                // ar[results_archive_path(e.sweep()) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
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
    
    void checkpoint_state(MPS<Matrix, SymmGroup> const& state, int sweep, int site)
    {
        status_type status;
        status["sweep"] = sweep;
        status["site"]  = site;
        return base::checkpoint_state(state, status);
    }

};

#endif
