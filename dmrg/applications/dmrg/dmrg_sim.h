/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_SIM_H
#define APP_DMRG_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include <boost/filesystem.hpp>

#include "utils/data_collector.hpp"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"

#ifdef USE_MTM
#include "types/mt_matrix/algorithms.hpp"
#endif

#ifdef USE_GPU
#include <cublas.h>
#endif

#include <alps/hdf5.hpp>

#include "dmrg/utils/DmrgParameters2.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "dmrg/utils/stream_storage.h"
#include "dmrg/utils/logger.h"
#include "dmrg/utils/random.hpp"
#include "utils/timings.h"

#include "dmrg/mp_tensors/ss_optimize.h"

#include "dmrg/models/factory.h"

using namespace app;


//template <class Matrix, class SymmGroup>
//struct do_fix_density {
//    
//    static void apply (MPS<Matrix, SymmGroup> & mps, std::vector<block_matrix<Matrix, SymmGroup> > const & dens_ops,
//                       typename SymmGroup::charge initc)
//    { }
//};
//template <class Matrix>
//struct do_fix_density<Matrix, TwoU1> {
//    
//    static void apply (MPS<Matrix, TwoU1> & mps, std::vector<block_matrix<Matrix, TwoU1> > const & dens_ops,
//                       TwoU1::charge initc)
//    {
//        std::vector<std::vector<double> > dens(2, std::vector<double>(mps.length(), 0.));
//        std::fill(dens[0].begin(), dens[0].end(), double(initc[0])/(mps.length()/2));
//        std::fill(dens[1].begin(), dens[1].end(), double(initc[1])/(mps.length()/2));
//        
//        fix_density(mps, dens_ops, dens);
//    }
//};

template <class Matrix, class SymmGroup>
class dmrg_sim {
        
    typedef std::vector<MPOTensor<Matrix, SymmGroup> > mpo_t;
    typedef Boundary<Matrix, SymmGroup> boundary_t;
    
public:
    
    dmrg_sim (DmrgParameters & parms_, ModelParameters & model_)
    : parms(parms_)
    , model(model_)
    , sweep(0)
    , chkpfile(parms.get<std::string>("chkpfile"))
    , rfile(parms.get<std::string>("resultfile"))
    , ssm(parms.get<std::string>("storagedir"))
    , dns( (parms.get<int>("donotsave") != 0) )
    { 
        DCOLLECTOR_GROUP(gemm_collector, "init")
        DCOLLECTOR_GROUP(svd_collector, "init")
#ifdef USE_GPU
        cublasInit();
#endif
        bool restore = false;
        {
			boost::filesystem::path p(chkpfile);
			if (boost::filesystem::exists(p) && boost::filesystem::is_regular_file(p))
            /* struct stat tmp;
            if (stat(chkpfile.c_str(), &tmp) == 0 && S_ISREG(tmp.st_mode)) */
            {
                cout << "Restoring state." << endl;
                restore = true;
            }
        }
        
#ifndef WIN32
		// drand48 does not seem to be used anymore anyway
        srand48(parms.get<int>("seed"));
#endif
        dmrg_random::engine.seed(parms.get<int>("seed"));
        
        typename SymmGroup::charge initc;
        model_parser(parms.get<std::string>("lattice_library"), parms.get<std::string>("model_library"), model, lat, H, initc, measurements);
        phys = H.get_phys();
        
        /*
         std::cout << "initc: " << initc << std::endl;
         std::cout << "phys:" << std::endl << phys << std::endl;
         std::cout << measurements << std::endl;
         std::cout << "Hamiltonian:" << std::endl << H << std::endl;
         */
        
        if (!parms.get<std::string>("always_measure").empty()) {
            meas_always.set_identity(measurements.get_identity());
            std::vector<std::string> meas_list = parms.get<std::vector<std::string> >("always_measure");
            for (int i=0; i<meas_list.size(); ++i)
                meas_always.add_term(measurements.get(meas_list[i]));
        }
        
        mpo = make_mpo(lat->size(), H);
        mpoc = mpo;
        if (parms.get<int>("use_compressed") > 0)
            mpoc.compress(1e-12);
        
        mps = MPS<Matrix, SymmGroup>(lat->size(),
                               parms.get<std::size_t>("init_bond_dimension"),
                               phys, initc,
                               *initializer_factory(parms));
        
        if (restore) {
            alps::hdf5::archive h5ar_in(chkpfile);
            h5ar_in >> alps::make_pvp("/state", mps);
            h5ar_in >> alps::make_pvp("/status/sweep", sweep);
            ++sweep;
        } else if (parms.get<std::string>("initfile").size() > 0) {
            alps::hdf5::archive h5ar_in(parms.get<std::string>("initfile"));
            h5ar_in >> alps::make_pvp("/state", mps);
        }
        
        gettimeofday(&now, NULL);
        
    }
    
    
    void optimize ()
    {
        bool early_exit = false;
        
        ss_optimize<Matrix, SymmGroup, StreamStorageMaster> optimizer(mps,
                                                                parms.get<int>("use_compressed") == 0 ? mpo : mpoc,
                                                                parms, ssm);
        
        
        for ( ; sweep < parms.get<int>("nsweeps"); ++sweep) {
            DCOLLECTOR_GROUP(gemm_collector, "sweep"+boost::lexical_cast<std::string>(sweep))
            DCOLLECTOR_GROUP(svd_collector, "sweep"+boost::lexical_cast<std::string>(sweep))
            gettimeofday(&snow, NULL);
            
            Logger iteration_log;
            
            optimizer.sweep(sweep, iteration_log);
            ssm.sync();
            
            MPS<Matrix, SymmGroup> cur_mps = optimizer.get_current_mps();
            
            std::vector<double> entropies = calculate_bond_entropies(cur_mps);
            
            gettimeofday(&sthen, NULL);
            double elapsed = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
            
            {
                alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                
                h5ar << alps::make_pvp("/parameters", parms);
                h5ar << alps::make_pvp("/parameters", model);
                
                std::ostringstream oss;
                
                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results";
                h5ar << alps::make_pvp(oss.str().c_str(), iteration_log);
                
                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results/Iteration Entropies/mean/value";
                h5ar << alps::make_pvp(oss.str().c_str(), entropies);
                
                cout << "Sweep done after " << elapsed << " seconds." << endl;
                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results/Runtime/mean/value";
                h5ar << alps::make_pvp(oss.str().c_str(), std::vector<double>(1, elapsed));                
            }
            {
                std::ostringstream oss;
                oss << "/simulation/sweep" << sweep << "/results/";
                if (meas_always.n_terms() > 0)
                    measure_on_mps(cur_mps, *lat, meas_always, rfile, oss.str());
            }
            
            if (!dns)
            {
                alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                
                h5ar << alps::make_pvp("/parameters", parms);
                h5ar << alps::make_pvp("/parameters", model);
                h5ar << alps::make_pvp("/state", cur_mps);
                h5ar << alps::make_pvp("/status/sweep", sweep);
            }
            
            mps = cur_mps;
            
            gettimeofday(&then, NULL);
            elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);            
            int rs = parms.get<int>("run_seconds");
            if (rs > 0 && elapsed > rs) {
                early_exit = true;
                break;
            }
        }
        
        ssm.sync();
    }
    
    
    void measure ()
    {
        
        {
            alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);

            h5ar << alps::make_pvp("/parameters", parms);
            h5ar << alps::make_pvp("/parameters", model);
        }

        cout << "Measurements." << endl;
        measure_on_mps(mps, *lat, measurements, rfile);
        
        Timer tvn("vN entropy"), tr2("Renyi n=2");
        cout << "Calculating vN entropy." << endl;
        tvn.begin(); std::vector<double> entropies = calculate_bond_entropies(mps); tvn.end();
        cout << "Calculating n=2 Renyi entropy." << endl;
        tr2.begin(); std::vector<double> renyi2 = calculate_bond_renyi_entropies(mps, 2); tr2.end();
        
        {
            alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            if (entropies.size() > 0)
                h5ar << alps::make_pvp("/spectrum/results/Entropy/mean/value", entropies);
            if (renyi2.size() > 0)
                h5ar << alps::make_pvp("/spectrum/results/Renyi2/mean/value", renyi2);
        }
        
        double energy = expval(mps, mpoc);
        cout << "Energy before: " << expval(mps, mpo) << endl;
        cout << "Energy after: " << expval(mps, mpoc) << endl;
        {
            alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            h5ar << alps::make_pvp("/spectrum/results/Energy/mean/value", std::vector<double>(1, energy));
        }
        
        if (parms.get<int>("calc_h2") > 0) {
            Timer tt1("square"), tt2("compress");
            tt1.begin(); MPO<Matrix, SymmGroup> mpo2 = square_mpo(mpo); tt1.end();
            tt2.begin(); mpo2.compress(1e-12); tt2.end();
            
            Timer t3("expval mpo2"), t4("expval mpo2c");
            
            t4.begin();
            double energy2 = expval(mps, mpo2, true);
            t4.end();
            
            cout << "Energy^2: " << energy2 << endl;
            cout << "Variance: " << energy2 - energy*energy << endl;
            
            {
                alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                h5ar << alps::make_pvp("/spectrum/results/Energy^2/mean/value", std::vector<double>(1, energy2));
                h5ar << alps::make_pvp("/spectrum/results/EnergyVariance/mean/value",
                                       std::vector<double>(1, energy2 - energy*energy));
            }
        }
    }
    
    
    
    ~dmrg_sim()
    {
#ifdef USE_GPU
        cublasShutdown();
#endif
    }
    
    
private:
    
    mps_initializer<Matrix, SymmGroup> * initializer_factory(BaseParameters & params)
    {
        if (params.get<std::string>("init_state") == "default")
            return new default_mps_init<Matrix, SymmGroup>();
        else if (params.get<std::string>("init_state") == "const")
            return new const_mps_init<Matrix, SymmGroup>();
        else if (params.get<std::string>("init_state") == "thin")
            return new thin_mps_init<Matrix, SymmGroup>();
        else if (params.get<std::string>("init_state") == "thin_const")
            return new thin_const_mps_init<Matrix, SymmGroup>();
        else {
            throw std::runtime_error("Don't know this initial state.");
            return NULL;
        }
    }
    
    
    DmrgParameters & parms;
    ModelParameters & model;
    
    timeval now, then, snow, sthen;

	bool dns;
    int sweep;
    std::string chkpfile;
    std::string rfile;
    
    Lattice * lat;
    Hamiltonian<Matrix, SymmGroup> H;
    Index<SymmGroup> phys;
    MPS<Matrix, SymmGroup> mps;
    MPO<Matrix, SymmGroup> mpo, mpoc;
    Measurements<Matrix, SymmGroup> measurements;
    Measurements<Matrix, SymmGroup> meas_always;

    StreamStorageMaster ssm;
    
};



#endif
