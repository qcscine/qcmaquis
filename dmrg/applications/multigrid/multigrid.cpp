//
//  multigrid.cpp
//  
//
//  Created by Michele Dolfi on 01.07.11.
//  Copyright 2011 ETHZ. All rights reserved.
//

#include "dmrg_version.h"

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

#include <boost/shared_ptr.hpp>

using std::cerr;
using std::cout;
using std::endl;

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
//#include "types/dense_matrix/dense_matrix_algorithms.h"
#include "types/dense_matrix/algorithms/algorithms.hpp"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"

#ifdef USE_GPU
#include <cublas.h>
#endif

#ifdef USE_MTM
#define USE_MTM_MAIN
#include "types/mt_matrix/mt_matrix.h"
#include "types/mt_matrix/algorithms.hpp"
typedef maquis::types::mt_matrix<double> Matrix;
#else
typedef maquis::types::dense_matrix<double> Matrix;
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
#include "dmrg/mp_tensors/multigrid.h"

#include "dmrg/utils/stream_storage.h"
#include "dmrg/utils/logger.h"
#include "dmrg/utils/random.hpp"

#include "dmrg/mp_tensors/optimize.h"

#include "dmrg/models/factory.h"
#include "dmrg/models/continuum/factory.h"
#include "dmrg/models/continuum/lattice.hpp"
//#ifdef UseTwoU1
//#include "dmrg/models/continuum/models_2u1.hpp"
//#else
//#ifdef UseNULL
//#include "dmrg/models/continuum/models_none.hpp"
//#else
//#include "dmrg/models/continuum/models_u1.hpp"
//#endif
//#endif

using namespace app;

#ifdef UseTwoU1
typedef TwoU1 grp;
#else
#ifdef UseNULL
typedef TrivialGroup grp;
#else
typedef U1 grp;
#endif
#endif

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

template<class Matrix>
mps_initializer<Matrix, grp> * initializer_factory(BaseParameters & params)
{
    if (params.get<std::string>("init_state") == "default")
        return new default_mps_init<Matrix, grp>();
    else if (params.get<std::string>("init_state") == "const")
        return new const_mps_init<Matrix, grp>();
    else if (params.get<std::string>("init_state") == "thin")
        return new thin_mps_init<Matrix, grp>();
    else if (params.get<std::string>("init_state") == "thin_const")
        return new thin_const_mps_init<Matrix, grp>();
    else {
        throw std::runtime_error("Don't know this initial state.");
        return NULL;
    }
}


MPO<Matrix, grp> mixed_mpo (BaseParameters & parms1, int L1, BaseParameters & parms2, int L2)
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
    
//    std::cout << "MIXED LATTICE ( " << L1 << ", " <<  L2 << " )" << std::endl;
//    for (int p=0; p<lat->size(); ++p) {
//        std::cout << lat->get_prop<std::string>("label", p, p+1) << ": " << lat->get_prop<double>("dx", p, p+1) << std::endl;
//        std::cout << lat->get_prop<std::string>("label", p, p-1) << ": " << lat->get_prop<double>("dx", p, p-1) << std::endl;
//    }
    
    model_traits<Matrix, grp>::model_ptr model = cont_model_factory<Matrix, grp>::parse(*lat, parms1);
    Hamiltonian<Matrix, grp> H = model->H();
    MPO<Matrix, grp> mpo = make_mpo(lat->size(), H);
        
    return mpo;
}

int main(int argc, char ** argv)
{
    cout << DMRG_VERSION_STRING << endl;
    if (argc != 3)
    {
        cout << "Usage: <parms> <model_parms>" << endl;
        exit(1);
    }
    
    cout.precision(10);
    
#ifdef USE_GPU
	cublasInit();
#endif
    
    std::ifstream param_file(argv[1]);
    if (!param_file) {
        cerr << "Could not open parameter file. (" << argv[1] << ")" << endl;
        exit(1);
    }
    DmrgParameters raw_parms(param_file);
    
    std::ifstream model_file(argv[2]);
    if (!model_file) {
        cerr << "Could not open model file. (" << argv[2] << ")" << endl;
        exit(1);
    }
    ModelParameters raw_model(model_file);
    
    srand48(raw_parms.get<int>("seed"));
    dmrg_random::engine.seed(raw_parms.get<int>("seed"));
    
    std::string chkpfile = raw_parms.get<std::string>("chkpfile");
    std::string rfile = raw_parms.get<std::string>("resultfile");
    bool dns = (raw_parms.get<int>("donotsave") != 0);
    
    int graining = 0;
    int sweep = 0;
    bool restore = false;
    {
        struct stat tmp;
        if (stat(chkpfile.c_str(), &tmp) == 0 && S_ISREG(tmp.st_mode))
        {
            cout << "Restoring state." << endl;
            restore = true;
            
            alps::hdf5::archive h5ar_in(chkpfile);
            h5ar_in >> alps::make_pvp("/status/sweep", sweep);
            ++sweep;
            h5ar_in >> alps::make_pvp("/status/graining", graining);
            
        }
    }
    
    
    {
        alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        h5ar << alps::make_pvp("/parameters", raw_parms);
        h5ar << alps::make_pvp("/parameters", raw_model);
        h5ar << alps::make_pvp("/version", DMRG_VERSION_STRING);
    }
    if (!dns) {
        alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        h5ar << alps::make_pvp("/parameters", raw_parms);
        h5ar << alps::make_pvp("/parameters", raw_model);
        h5ar << alps::make_pvp("/version", DMRG_VERSION_STRING);
    }
    
    Lattice_ptr lat;
    model_traits<Matrix, grp>::model_ptr phys_model;
    grp::charge initc;
    Hamiltonian<Matrix, grp> H;
    Measurements<Matrix, grp> measurements;
    Index<grp> phys;
    
    MPO<Matrix, grp> mpo(0), mpoc(0);
    MPS<Matrix, grp> cur_mps;
    
    BaseParameters old_model;
    {
        int t_graning = (graining > 0) ? graining-1 : 0;
                
        BaseParameters parms = raw_parms.get_at_index("graining", t_graning);
        BaseParameters model = raw_model.get_at_index("graining", t_graning);
        
        old_model = model;
        
        model_parser<Matrix, grp>("continuum", "continuum", model, lat, phys_model);
        H = phys_model->H();
        measurements = phys_model->measurements();
        initc = phys_model->initc(model);
        phys = H.get_phys();
        std::cout << "initc: " << initc << std::endl;
        
        mpo = make_mpo(lat->size(), H);
        mpoc = mpo;
        if (parms.get<int>("use_compressed") > 0)
            mpoc.compress(1e-12);
    }
    
    if (restore) {
        alps::hdf5::archive h5ar_in(chkpfile);
        h5ar_in >> alps::make_pvp("/state", cur_mps);
    } else if (raw_parms.get<std::string>("initfile").size() > 0) {
        alps::hdf5::archive h5ar_in(raw_parms.get<std::string>("initfile"));
        h5ar_in >> alps::make_pvp("/state", cur_mps);
    }
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    do {
        
        BaseParameters parms = raw_parms.get_at_index("graining", graining);
        BaseParameters model = raw_model.get_at_index("graining", graining);
        
        model_parser<Matrix, grp>("continuum", "continuum", model, lat, phys_model);
        H = phys_model->H();
        measurements = phys_model->measurements();
        phys = H.get_phys();

#ifndef NDEBUG
//        std::cout << parms << std::endl;
//        std::cout << model << std::endl;
//        std::cout << measurements << std::endl;
//        std::cout << H << std::endl;
//
//        std::cout << "LATTICE:" << std::endl;
//        for (int i=0; i<lat->size(); ++i)
//            std::cout << i << ": " << lat->get_prop<double>("x", i) << std::endl;
#endif
        
        Measurements<Matrix, grp> meas_always;
        if (!parms.get<std::string>("always_measure").empty()) {
            meas_always.set_identity(measurements.get_identity());
            std::vector<std::string> meas_list = parms.get<std::vector<std::string> >("always_measure");
            for (int i=0; i<meas_list.size(); ++i)
                meas_always.add_term(measurements.get(meas_list[i]));
        }
        
        MPO<Matrix, grp> t_mpo = make_mpo(lat->size(), H);
        MPO<Matrix, grp> t_mpoc = t_mpo;
        if (parms.get<int>("use_compressed") > 0)
            t_mpoc.compress(1e-12);
        
        
        MPS<Matrix, grp> initial_mps;

        
        static Timer multigrid_t("Multigrid");
        if (cur_mps.length() > 0 && cur_mps.length() != lat->size())
        {
            multigrid_t.begin();
            std::cout << "*** Starting grainings ***" << std::endl;
            Logger iteration_log;
            
            boost::shared_ptr<mps_initializer<Matrix, grp> > initializer = boost::shared_ptr<mps_initializer<Matrix, grp> > (new empty_mps_init<Matrix, grp>());
            initial_mps = MPS<Matrix, grp>(lat->size(), 1, phys, initc, *initializer);
                        
            int oldL = old_model.get<double>("Ndiscr") * old_model.get<double>("L");
            std::vector<MPO<Matrix, grp> > mpo_mix(oldL+1, MPO<Matrix, grp>(0));
            double r = model.get<double>("Ndiscr")/old_model.get<double>("Ndiscr");
            for (int i=0; i<=oldL; ++i)
                mpo_mix[i] = mixed_mpo(model, r*i, old_model, oldL-i);
            
//            std::cout << "Old MPS:" << std::endl << initial_mps.description() << std::endl;
            if (cur_mps.length() < initial_mps.length())
                multigrid::extension_optim(parms, iteration_log,
                                           cur_mps, initial_mps, mpo_mix);
            else if (cur_mps.length() > initial_mps.length())
                multigrid::restriction(cur_mps, initial_mps);
//            std::cout << "New MPS:" << std::endl << initial_mps.description();
            multigrid_t.end();

            std::vector<double> energies, entropies;            
            entropies = calculate_bond_entropies(initial_mps);
                        
            {
                alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                
                std::ostringstream oss;
                
                oss.str("");
                oss << "/simulation/iteration/graining/" << graining << "/parameters";
                h5ar << alps::make_pvp(oss.str(), parms);
                h5ar << alps::make_pvp(oss.str(), model);
                
                oss.str("");
                oss << "/simulation/iteration/graining/" << graining << "/results";
                
                h5ar << alps::make_pvp(oss.str(), iteration_log);
                h5ar << alps::make_pvp(oss.str()+"/Iteration Entropies/mean/value", entropies);
            }
            {
                std::ostringstream oss;
                oss << "/simulation/iteration/graining/" << graining << "/results/";
                if (meas_always.n_terms() > 0)
                    measure_on_mps(initial_mps, *lat, meas_always, rfile, oss.str());
            }            
        } else if (cur_mps.length() > 0) {
            initial_mps = cur_mps;
        } else {
            initial_mps = MPS<Matrix, grp>(lat->size(),
                                         parms.get<std::size_t>("init_bond_dimension"),
                                         phys, initc,
                                         *initializer_factory<Matrix>(parms));
        }
        
        cur_mps = initial_mps;
        mpo = t_mpo;
        mpoc = t_mpoc;
        old_model = model;
        
#ifndef MEASURE_ONLY
        
        //    BaseStorageMaster * bsm = bsm_factory(parms);
        StreamStorageMaster ssm(parms.get<std::string>("storagedir"));
        
        
        bool early_exit = false;
        {   
            std::cout << "*** Starting optimization ***" << std::endl;
            ss_optimize<Matrix, grp, StreamStorageMaster> optimizer(initial_mps,
                                                                    parms.get<int>("use_compressed") == 0 ? mpo : mpoc,
                                                                    parms, ssm);
            
            for ( ; sweep < parms.get<int>("nsweeps"); ++sweep) {
                gettimeofday(&snow, NULL);
                
                Logger iteration_log;
                
                optimizer.sweep(sweep, iteration_log);
                ssm.sync();
                
                cur_mps = optimizer.get_current_mps();
                
                std::vector<double> energies, entropies, renyi2;
                std::vector<std::size_t> truncations;

                entropies = calculate_bond_entropies(cur_mps);
                
                gettimeofday(&sthen, NULL);
                double elapsed = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
                cout << "Sweep done after " << elapsed << " seconds." << endl;
                
                {
                    alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                    
                    std::ostringstream oss;
                    
                    oss.str("");
                    oss << "/simulation/iteration/graining/" << graining << "/parameters";
                    h5ar << alps::make_pvp(oss.str(), parms);
                    h5ar << alps::make_pvp(oss.str(), model);
                    
                    oss.str("");
                    oss << "/simulation/iteration/graining/" << graining << "/sweep" << sweep << "/results";
                    
                    h5ar << alps::make_pvp(oss.str(), iteration_log);                    
                    h5ar << alps::make_pvp(oss.str()+"/Iteration Entropies/mean/value", entropies);
                    h5ar << alps::make_pvp(oss.str()+"/Runtime/mean/value", std::vector<double>(1, elapsed));                
                }
                {
                    std::ostringstream oss;
                    oss << "/simulation/iteration/graining/" << graining << "/sweep" << sweep << "/results/";
                    if (meas_always.n_terms() > 0)
                        measure_on_mps(cur_mps, *lat, meas_always, rfile, oss.str());
                }
                
                if (!dns)
                {
                    alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                    
                    h5ar << alps::make_pvp("/state", cur_mps);
                    h5ar << alps::make_pvp("/status/sweep", sweep);
                    h5ar << alps::make_pvp("/status/graining", graining);
                }
                
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
        sweep = 0;
#endif

        ++graining;
    } while (graining < raw_parms.get<int>("ngrainings"));
    
#ifdef MEASURE_ONLY
    {
        MPS<Matrix, grp> mps = cur_mps;
        
        std::vector<double> energies, entropies, renyi2;
        std::vector<std::size_t> truncations;
        
        cout << "Measurements." << endl;
        measure_on_mps(mps, *lat, measurements, rfile);
        
        Timer tvn("vN entropy"), tr2("Renyi n=2");
        cout << "Calculating vN entropy." << endl;
        tvn.begin(); entropies = calculate_bond_entropies(mps); tvn.end();
        cout << "Calculating n=2 Renyi entropy." << endl;
        tr2.begin(); renyi2 = calculate_bond_renyi_entropies(mps, 2); tr2.end();
        
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
        
        if (raw_parms.get<int>("calc_h2") > 0) {
            Timer tt1("square"), tt2("compress");
            tt1.begin(); MPO<Matrix, grp> mpo2 = square_mpo(mpo); tt1.end();
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
#endif
    
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
    
    cout << "Task took " << elapsed << " seconds." << endl;
    
#ifdef USE_GPU
	cublasShutdown();
#endif
}
