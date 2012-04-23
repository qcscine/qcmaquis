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

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/sim/sim.h"
#include "dmrg/models/continuum/factory.h"

#include "dmrg/mp_tensors/optimize.h"
#include "dmrg/mp_tensors/multigrid.h"

using namespace app;
template <class Matrix, class SymmGroup>
class multigrid_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, StreamStorageMaster> opt_base_t;
    
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
        assert(parms_orig.get<std::string>("lattice_library") == "continuum");
        assert(parms_orig.get<std::string>("model_library") == "continuum");
        
        if (this->restore)
        {
            alps::hdf5::archive h5ar_in(this->chkpfile);
            h5ar_in >> alps::make_pvp("/status/graining", graining);
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
    
    int do_sweep (Logger& iteration_log, double time_limit = -1)
    {
        int exit_site = optimizer->sweep(base::sweep, iteration_log, Both,
                                         base::site, time_limit);
        base::ssm.sync();
        
        base::mps = optimizer->get_current_mps();
        
        return exit_site;
    }
    
    void run ()
    {
        bool early_exit = false;
        do {
            MPS<Matrix, SymmGroup> initial_mps;
            
            
            static Timer multigrid_t("Multigrid");
            m_type = mg_measure;

            if (cur_graining != graining)
            {
                multigrid_t.begin();
                cout << "*** Starting grainings ***" << std::endl;
                Logger iteration_log;
                
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
                
                //            cout << "Old MPS:" << std::endl << initial_mps.description() << std::endl;
                if (curL < initial_mps.length())
                    multigrid::extension_optim(base::parms, iteration_log,
                                               this->mps, initial_mps, mpo_mix);
                else if (this->mps.length() > initial_mps.length())
                    multigrid::restriction(this->mps, initial_mps);
                //            cout << "New MPS:" << std::endl << initial_mps.description();
                multigrid_t.end();
                
                this->mps = initial_mps;
                cur_graining = graining;
                
                this->do_sweep_measure(iteration_log);
                { // TODO: port this to a function in the base class!
                    alps::hdf5::archive h5ar(this->rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                    
                    h5ar << alps::make_pvp(this->sweep_archive_path() + "/parameters", this->parms);
                    h5ar << alps::make_pvp(this->sweep_archive_path() + "/parameters", this->model);
                    
                    h5ar << alps::make_pvp(this->sweep_archive_path() + "/results", iteration_log);
                }
                if (!this->dns)
                {
                    alps::hdf5::archive h5ar(this->chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                    
                    h5ar << alps::make_pvp("/state", this->mps);
                    h5ar << alps::make_pvp("/status/sweep", this->sweep);
                    h5ar << alps::make_pvp("/status/graining", this->cur_graining);
                    h5ar << alps::make_pvp("/status/site", -1);
                }
            } else {
                { // TODO: port this to a function in the base class!
                    alps::hdf5::archive h5ar(this->rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                    
                    h5ar << alps::make_pvp(this->sweep_archive_path() + "/parameters", this->parms);
                    h5ar << alps::make_pvp(this->sweep_archive_path() + "/parameters", this->model);
                }
            }
                
            
            // OPTIM
            m_type = sweep_measure;
            init_optimizer();
            early_exit = base::exec_sweeps();
            if (!this->dns)
            {
                alps::hdf5::archive h5ar(this->chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                h5ar << alps::make_pvp("/status/graining", this->cur_graining);
            }
            if (early_exit)
                break;
            
            ++graining;
            this->sweep = 0;
        } while (graining < parms_orig.get<int>("ngrainings"));
        
    }
    
    
    ~multigrid_sim()
    {
        base::ssm.sync();
    }    
private:
    
    void init_optimizer()
    {
        if (base::parms.template get<std::string>("optimization") == "singlesite")
        {
            optimizer = 
            boost::shared_ptr<opt_base_t> ( new ss_optimize<Matrix, SymmGroup, StreamStorageMaster>
                                           (base::mps,
                                            base::parms.template get<int>("use_compressed") == 0 ? base::mpo : base::mpoc,
                                            base::parms, base::ssm) );
        } 
        
        else if (base::parms.template get<std::string>("optimization") == "twosite")
        {
            optimizer = 
            boost::shared_ptr<opt_base_t> ( new ts_optimize<Matrix, SymmGroup, StreamStorageMaster>
                                           (base::mps,
                                            base::parms.template get<int>("use_compressed") == 0 ? base::mpo : base::mpoc,
                                            base::parms, base::ssm) );
        }
        
        else
        {
            throw std::runtime_error("Don't know this optimization!");
        }
    }
    
    MPO<Matrix, SymmGroup> mixed_mpo (BaseParameters & parms1, int L1, BaseParameters & parms2, int L2)
    {
        assert( parms1.get<std::string>("LATTICE") == parms2.get<std::string>("LATTICE") );
        
        Lattice_ptr lat;
        if (parms1.get<std::string>("LATTICE") == "continuous_chain"
            || parms1.get<std::string>("LATTICE") == std::string("continuous_left_chain"))
            lat = Lattice_ptr(new MixedContChain(parms1, L1, parms2, L2));
        else if (parms2.get<std::string>("LATTICE") == std::string("continuous_center_chain"))
            lat = Lattice_ptr(new MixedContChain_c(parms1, L1, parms2, L2));
        else
            throw std::runtime_error("Don't know this lattice!");
        
#ifndef NDEBUG
        // debugging output, to be removed soon!
        cout << "MIXED LATTICE ( " << L1 << ", " <<  L2 << " )" << std::endl;
        for (int p=0; p<lat->size(); ++p) {
            cout << lat->get_prop<std::string>("label", p) << ": " << lat->get_prop<double>("dx", p) << std::endl;
            cout << lat->get_prop<std::string>("label", p, p+1) << ": " << lat->get_prop<double>("dx", p, p+1) << std::endl;
            cout << lat->get_prop<std::string>("label", p, p-1) << ": " << lat->get_prop<double>("dx", p, p-1) << std::endl;
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
