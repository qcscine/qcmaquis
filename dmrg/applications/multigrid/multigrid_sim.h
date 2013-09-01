/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_MULTIGRID_SIM_H
#define APP_MULTIGRID_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/shared_ptr.hpp>

#include "dmrg/sim/sim.h"
#include "dmrg/models/continuum/factory.h"

#include "dmrg/mp_tensors/optimize.h"
#include "dmrg/mp_tensors/multigrid.h"

template <class Matrix, class SymmGroup>
class multigrid_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    typedef typename base::status_type status_type;

    enum measure_t {sweep_measure, mg_measure};
    
    using base::mps;
    using base::mpo;
    using base::lat;
    using base::mpoc;
    using base::parms;
    using base::model;
    using base::measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::init_site;
    using base::rfile;
    
public:
    multigrid_sim(DmrgParameters & parms_, ModelParameters & model_)
    : base(parms_, model_)
    , initial_graining(0)
    {
        assert(parms["lattice_library"] == "continuum");
        assert(parms["model_library"] == "continuum");
        
        if (this->restore)
        {
            storage::archive ar(this->chkpfile+"/props.h5");
            ar["/status/graining"] >> initial_graining;
        }
    }
    
    void run()
    {
        /// Set current status in parms
        parms = parms.get_at_index("graining", initial_graining);
        model = model.get_at_index("graining", initial_graining);
        /// Build current model and load/build MPS
        this->model_init();
        this->mps_init();
        
        for (int graining=initial_graining; graining < parms["ngrainings"]; ++graining)
        {
            /// usual optimization
            if (init_sweep < parms["nsweeps"])
                dmrg_run(graining);
            
            if ( stop_callback() ) {
                maquis::cout << "Time limit reached." << std::endl;
                break;
            }

            /// fine graining
            if (graining < parms["ngrainings"]-1) {
                maquis::cout << "*** Starting grainings ***" << std::endl;
                
                parms = parms.get_at_index("graining", graining+1);
                model = model.get_at_index("graining", graining+1);
                this->model_init();
                
                boost::shared_ptr<mps_initializer<Matrix, SymmGroup> > initializer = boost::shared_ptr<mps_initializer<Matrix, SymmGroup> > (new empty_mps_init<Matrix, SymmGroup>());
                MPS<Matrix, SymmGroup> new_mps = MPS<Matrix, SymmGroup>(lat->size(), 1, this->phys, this->initc, *initializer);
                
                int curL = mps.length();
                BaseParameters oldmodel = model.get_at_index("graining", graining);
                
                std::vector<MPO<Matrix, SymmGroup> > mpo_mix(curL+1, MPO<Matrix, SymmGroup>(0));
                double r = lat->size() / curL;
                for (int i=0; i<=curL; ++i)
                    mpo_mix[i] = mixed_mpo(base::model, r*i, oldmodel, curL-i);
                
                results_collector graining_results;
                if (curL < new_mps.length())
                    graining_results = multigrid::extension_optim(base::parms, mps, new_mps, mpo_mix);
                else if (this->mps.length() > new_mps.length())
                    throw std::runtime_error("Restriction operation not really implemented.");
                // graining_results = multigrid::restriction(this->mps, initial_mps);
                
                /// swap mps
                swap(mps, new_mps);
                
                
                /// write iteration results
                {
                    storage::archive ar(rfile, "w");
                    ar[results_archive_path(0, graining, mg_measure) + "/parameters"] << parms;
                    ar[results_archive_path(0, graining, mg_measure) + "/parameters"] << model;
                    ar[results_archive_path(0, graining, mg_measure) + "/results"] << graining_results;
                }
                
                /// measure observables specified in 'always_measure'
                if (!parms["always_measure"].empty())
                    this->measure(results_archive_path(0, graining, mg_measure) + "/results/", measurements.sublist(parms["always_measure"]));
                
                /// checkpoint new mps
                this->checkpoint_simulation(mps, 0, -1, graining+1);
            }
            
            if ( stop_callback() ) {
                maquis::cout << "Time limit reached." << std::endl;
                break;
            }
        }
    }
    
    ~multigrid_sim()
    {
        storage::disk::sync();
    }
    
private:
    
    std::string results_archive_path(status_type const& status) const
    {
        throw std::runtime_error("do not use in multigrid.");
    }
    
    void checkpoint_simulation(MPS<Matrix, SymmGroup> const& state, int sweep, int site, int graining)
    {
        status_type status;
        status["sweep"]    = sweep;
        status["site"]     = site;
        status["graining"] = graining;
        return base::checkpoint_simulation(state, status);
    }

    std::string results_archive_path(int sweep, int graining, measure_t m_type) const
    {
        std::ostringstream oss;
        oss.str("");
        switch (m_type) {
            case sweep_measure:
                oss << "/simulation/iteration/graining/" << graining << "/sweep" << sweep;
                break;
            case mg_measure:
                oss << "/simulation/iteration/graining/" << graining;
                break;
        }
        return oss.str();
    }
    
    void dmrg_run(int graining)
    {
        int meas_each = parms["measure_each"];
        int chkp_each = parms["chkp_each"];
        
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
                        ar[results_archive_path(sweep, graining, sweep_measure) + "/parameters"] << parms;
                        ar[results_archive_path(sweep, graining, sweep_measure) + "/parameters"] << model;
                        ar[results_archive_path(sweep, graining, sweep_measure) + "/results"] << optimizer->iteration_results();
                        // ar[results_archive_path(sweep, graining, sweep_measure) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
                    }
                    
                    /// measure observables specified in 'always_measure'
                    if (!parms["always_measure"].empty())
                        this->measure(results_archive_path(sweep, graining, sweep_measure) + "/results/", measurements.sublist(parms["always_measure"]));
                }
                
                /// write checkpoint
                bool stopped = stop_callback();
                if (stopped || (sweep+1) % chkp_each == 0 || (sweep+1) == parms["nsweeps"])
                    this->checkpoint_simulation(mps, sweep, -1, graining);
                
                if (stopped) break;
            }
        } catch (dmrg::time_limit const& e) {
            maquis::cout << e.what() << " checkpointing partial result." << std::endl;
            this->checkpoint_simulation(mps, e.sweep(), e.site(), graining);
            
            {
                storage::archive ar(rfile, "w");
                ar[results_archive_path(e.sweep(), graining, sweep_measure) + "/parameters"] << parms;
                ar[results_archive_path(e.sweep(), graining, sweep_measure) + "/parameters"] << model;
                ar[results_archive_path(e.sweep(), graining, sweep_measure) + "/results"] << optimizer->iteration_results();
                // ar[results_archive_path(e.sweep(), graining, sweep_measure) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
            }
        }
        
        /// for the next graining level
        init_sweep = 0;
        init_site  = -1;
    }

    
    MPO<Matrix, SymmGroup> mixed_mpo(BaseParameters & parms1, int L1, BaseParameters & parms2, int L2)
    {
        assert( parms1["LATTICE"] == parms2["LATTICE"] );
        
        Lattice_ptr lat;
        if (parms1["LATTICE"] == "continuous_chain"
            || parms1["LATTICE"] == std::string("continuous_left_chain"))
            lat = Lattice_ptr(new MixedContChain(parms1, L1, parms2, L2));
        else if (parms2["LATTICE"] == std::string("continuous_center_chain"))
            lat = Lattice_ptr(new MixedContChain_c(parms1, L1, parms2, L2));
        else
            throw std::runtime_error("Don't know this lattice!");
        
#ifndef NDEBUG
        // debugging output, to be removed soon!
//        maquis::cout << "MIXED LATTICE ( " << L1 << ", " <<  L2 << " )" << std::endl;
//        for (int p=0; p<lat->size(); ++p) {
//            maquis::cout << lat->get_prop<std::string>("label", p) << ": " << lat->get_prop<double>("dx", p) << std::endl;
//            maquis::cout << lat->get_prop<std::string>("label", p, p+1) << ": " << lat->get_prop<double>("dx", p, p+1) << std::endl;
//            maquis::cout << lat->get_prop<std::string>("label", p, p-1) << ": " << lat->get_prop<double>("dx", p, p-1) << std::endl;
//        }
#endif
         
        typename model_traits<Matrix, SymmGroup>::model_ptr tmpmodel;
        tmpmodel = cont_model_factory<Matrix, SymmGroup>::parse(*lat, parms1);
        Hamiltonian<Matrix, SymmGroup> H = tmpmodel->H();
        MPO<Matrix, SymmGroup> mpo = make_mpo(lat->size(), H);
        
        return mpo;
    }
    
private:
    int initial_graining;
};

#endif
