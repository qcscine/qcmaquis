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
    
    typedef std::vector<MPOTensor<Matrix, SymmGroup> > mpo_t;
    typedef Boundary<Matrix, SymmGroup> boundary_t;
    
    enum measure_t {sweep_measure, mg_measure};
    
public:
    
    multigrid_sim (DmrgParameters & parms_, ModelParameters & model_)
    : base(parms_, model_, false)
    , parms_orig(parms_)
    , model_orig(model_)
    , graining(0)
    , m_type(sweep_measure)
    {
        assert(parms_orig["lattice_library"] == "continuum");
        assert(parms_orig["model_library"] == "continuum");
        
        if (this->restore)
        {
            storage::archive ar(this->chkpfile);
            ar["/status/graining"] >> graining;
        }
        
        base::parms = parms_orig.get_at_index("graining", graining);
        base::model = model_orig.get_at_index("graining", graining);
        
        base::model_init();
        base::mps_init();
        cur_graining = graining;
    }
    
    std::string sweep_archive_path ()
    {
        std::ostringstream oss;
        oss.str("");
        switch (m_type) {
            case sweep_measure:
                oss << "/simulation/iteration/graining/" << graining << "/sweep" << this->sweep;
                break;
            case mg_measure:
                oss << "/simulation/iteration/graining/" << graining;
                break;
        }
        return oss.str();
    }
    
    int do_sweep (double time_limit = -1)
    {
        int exit_site = optimizer->sweep(base::sweep, Both,
                                         base::site, time_limit);
        storage::disk::sync();
        
        base::mps = optimizer->get_current_mps();
        
        return exit_site;
    }
    
    bool run ()
    {
        bool early_exit = false;
        do {
            MPS<Matrix, SymmGroup> initial_mps;
            
            
            m_type = mg_measure;

            if (cur_graining != graining)
            {
                maquis::cout << "*** Starting grainings ***" << std::endl;
                
                base::parms = parms_orig.get_at_index("graining", graining);
                base::model = model_orig.get_at_index("graining", graining);
                base::model_init();
                
                boost::shared_ptr<mps_initializer<Matrix, SymmGroup> > initializer = boost::shared_ptr<mps_initializer<Matrix, SymmGroup> > (new empty_mps_init<Matrix, SymmGroup>());
                initial_mps = MPS<Matrix, SymmGroup>(this->lat->size(), 1, this->phys, this->initc, *initializer);
                
                int curL = this->mps.length();
                BaseParameters oldmodel = model_orig.get_at_index("graining", graining-1);
                
                std::vector<MPO<Matrix, SymmGroup> > mpo_mix(curL+1, MPO<Matrix, SymmGroup>(0));
                double r = this->lat->size() / curL;
                for (int i=0; i<=curL; ++i)
                    mpo_mix[i] = mixed_mpo(base::model, r*i, oldmodel, curL-i);
                
                //            maquis::cout << "Old MPS:" << std::endl << initial_mps.description() << std::endl;
                if (curL < initial_mps.length())
                    multigrid::extension_optim(base::parms,
                                               this->mps, initial_mps, mpo_mix);
                else if (this->mps.length() > initial_mps.length())
                    throw std::runtime_error("Restriction operation not really implemented.");
//                    multigrid::restriction(this->mps, initial_mps);
                //            maquis::cout << "New MPS:" << std::endl << initial_mps.description();
                
                this->mps = initial_mps;
                cur_graining = graining;
                
                this->do_sweep_measure();
                { // TODO: port this to a function in the base class!
                    storage::archive ar(this->rfile, "w");
                    
                    ar[this->sweep_archive_path() + "/parameters"] << this->parms;
                    ar[this->sweep_archive_path() + "/parameters"] << this->model;
                    
                    ar[this->sweep_archive_path() + "/results"] << storage::log;
                }
                if (!this->dns)
                {
                    storage::archive ar(this->chkpfile, "w");
                    
                    ar["/state"] << this->mps;
                    ar["/status/sweep"] << this->sweep;
                    ar["/status/graining"] << this->cur_graining;
                    ar["/status/site"] << -1;
                }
            } else {
                { // TODO: port this to a function in the base class!
                    storage::archive ar(this->rfile, "w");
                    
                    ar[this->sweep_archive_path() + "/parameters"] << this->parms;
                    ar[this->sweep_archive_path() + "/parameters"] << this->model;
                }
            }
                
            
            // OPTIM
            m_type = sweep_measure;
            init_optimizer();
            early_exit = base::run();
            if (!this->dns)
            {
                storage::archive ar(this->chkpfile, "w");
                ar["/status/graining"] << this->cur_graining;
            }
            if (early_exit)
                break;
            
            ++graining;
            this->sweep = 0;
        } while (graining < parms_orig.get<int>("ngrainings"));
    
        return early_exit;
    }
    
    
    ~multigrid_sim()
    {
        storage::disk::sync();
    }    
private:
    
    void init_optimizer()
    {
        if (base::parms["optimization"] == "singlesite")
        {
            optimizer = 
            boost::shared_ptr<opt_base_t> ( new ss_optimize<Matrix, SymmGroup, storage::disk>
                                           (base::mps, base::mpoc,
                                            base::parms) );
        } 
        
        else if (base::parms["optimization"] == "twosite")
        {
            optimizer = 
            boost::shared_ptr<opt_base_t> ( new ts_optimize<Matrix, SymmGroup, storage::disk>
                                           (base::mps, base::mpoc, base::ts_cache_mpo,
                                            base::parms) );
        }
        
        else
        {
            throw std::runtime_error("Don't know this optimization!");
        }
    }
    
    MPO<Matrix, SymmGroup> mixed_mpo (BaseParameters & parms1, int L1, BaseParameters & parms2, int L2)
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
        maquis::cout << "MIXED LATTICE ( " << L1 << ", " <<  L2 << " )" << std::endl;
        for (int p=0; p<lat->size(); ++p) {
            maquis::cout << lat->get_prop<std::string>("label", p) << ": " << lat->get_prop<double>("dx", p) << std::endl;
            maquis::cout << lat->get_prop<std::string>("label", p, p+1) << ": " << lat->get_prop<double>("dx", p, p+1) << std::endl;
            maquis::cout << lat->get_prop<std::string>("label", p, p-1) << ": " << lat->get_prop<double>("dx", p, p-1) << std::endl;
        }
#endif
         
        typename model_traits<Matrix, SymmGroup>::model_ptr tmpmodel;
        tmpmodel = cont_model_factory<Matrix, SymmGroup>::parse(*lat, parms1);
        Hamiltonian<Matrix, SymmGroup> H = tmpmodel->H();
        MPO<Matrix, SymmGroup> mpo = make_mpo(lat->size(), H);
        
        return mpo;
    }
    
    BaseParameters parms_orig;
    BaseParameters model_orig;
    int graining, cur_graining;
    measure_t m_type;
    boost::shared_ptr<opt_base_t> optimizer;
};

#endif
