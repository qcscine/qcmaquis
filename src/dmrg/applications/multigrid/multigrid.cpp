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

#include "dmrg/block_matrix/detail/alps.hpp"

typedef alps::numeric::matrix<double> Matrix;

#include "dmrg/utils/DmrgParameters2.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/multigrid.h"

#include "dmrg/utils/storage.h"
#include "dmrg/utils/logger.h"
#include "dmrg/utils/random.hpp"

#include "dmrg/mp_tensors/optimize.h"

#include "dmrg/models/factory.h"
#include "dmrg/models/continuum/factory.h"
#include "dmrg/models/continuum/lattice.hpp"

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
    if (params["init_state"] == "default")
        return new default_mps_init<Matrix, grp>();
    else if (params["init_state"] == "const")
        return new const_mps_init<Matrix, grp>();
    else if (params["init_state"] == "thin")
        return new thin_mps_init<Matrix, grp>();
    else if (params["init_state"] == "thin_const")
        return new thin_const_mps_init<Matrix, grp>();
    else {
        throw std::runtime_error("Don't know this initial state.");
        return NULL;
    }
}


MPO<Matrix, grp> mixed_mpo (BaseParameters & parms1, int L1, BaseParameters & parms2, int L2)
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
    
//    maquis::cout << "MIXED LATTICE ( " << L1 << ", " <<  L2 << " )" << std::endl;
//    for (int p=0; p<lat->size(); ++p) {
//        maquis::cout << lat->get_prop<std::string>("label", p, p+1) << ": " << lat->get_prop<double>("dx", p, p+1) << std::endl;
//        maquis::cout << lat->get_prop<std::string>("label", p, p-1) << ": " << lat->get_prop<double>("dx", p, p-1) << std::endl;
//    }
    
    model_traits<Matrix, grp>::model_ptr model = cont_model_factory<Matrix, grp>::parse(*lat, parms1);
    Hamiltonian<Matrix, grp> H = model->H();
    MPO<Matrix, grp> mpo = make_mpo(lat->size(), H);
        
    return mpo;
}

int main(int argc, char ** argv)
{
    maquis::cout << DMRG_VERSION_STRING << std::endl;
#ifdef MEASURE_ONLY
    maquis::cout << "Only measuring." << std::endl;
#endif
#ifdef UseTwoU1
    maquis::cout << "TwoU1 symmetry" << std::endl;
#else
#ifdef UseNULL
    maquis::cout << "No symmetry" << std::endl;
#else
    maquis::cout << "U1 symmetry" << std::endl;
#endif
#endif

    if (argc != 3)
    {
        maquis::cout << "Usage: <parms> <model_parms>" << std::endl;
        exit(1);
    }
    
    maquis::cout.precision(10);
    
    std::ifstream param_file(argv[1]);
    if (!param_file) {
        maquis::cerr << "Could not open parameter file. (" << argv[1] << ")" << std::endl;
        exit(1);
    }
    DmrgParameters raw_parms(param_file);
    
    std::ifstream model_file(argv[2]);
    if (!model_file) {
        maquis::cerr << "Could not open model file. (" << argv[2] << ")" << std::endl;
        exit(1);
    }
    ModelParameters raw_model(model_file);
    
    srand48(raw_parms["seed"]);
    dmrg_random::engine.seed(raw_parms["seed"]);
    
    std::string chkpfile = raw_parms["chkpfile"];
    std::string rfile = raw_parms["resultfile"];
    bool dns = (raw_parms["donotsave"] != 0);
    
    int graining = 0;
    int sweep = 0;
    bool restore = false;
    {
        struct stat tmp;
        if (stat(chkpfile.c_str(), &tmp) == 0 && S_ISREG(tmp.st_mode))
        {
            maquis::cout << "Restoring state." << std::endl;
            restore = true;
            
            storage::archive ar(chkpfile);
            ar["/status/sweep"] >> sweep; ++sweep;
            ar["/status/graining"] >> graining;
        }
    }
    
    
    {
        storage::archive ar(rfile, "w");
        ar["/parameters"] << raw_parms;
        ar["/parameters"] << raw_model;
        ar["/version"] << DMRG_VERSION_STRING;
    }
    if (!dns) {
        storage::archive ar(chkpfile, "w");
        ar["/parameters"] << raw_parms;
        ar["/parameters"] << raw_model;
        ar["/version"] << DMRG_VERSION_STRING;
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
        maquis::cout << "initc: " << initc << std::endl;
//        maquis::cout << "Hamiltonian: " << H << std::endl;
        
        mpo = make_mpo(lat->size(), H);
        mpoc = mpo;
        if (parms["use_compressed"] > 0)
            mpoc.compress(1e-12);
    }
    
    if (restore) {
        storage::archive ar(chkpfile);
        ar["/state"] >> cur_mps;
    } else if (raw_parms["initfile"].size() > 0) {
        storage::archive ar(raw_parms["initfile"].str());
        ar["/state"] >> cur_mps;
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
//        maquis::cout << parms << std::endl;
//        maquis::cout << model << std::endl;
//        maquis::cout << measurements << std::endl;
//        maquis::cout << H << std::endl;
//
//        maquis::cout << "LATTICE:" << std::endl;
//        for (int i=0; i<lat->size(); ++i)
//            maquis::cout << i << ": " << lat->get_prop<double>("x", i) << std::endl;
#endif
        
        Measurements<Matrix, grp> meas_always;
        if (!parms["always_measure"].empty()) {
            meas_always.set_identity(measurements.get_identity());
            std::vector<std::string> meas_list = parms.get<std::vector<std::string> >("always_measure");
            for (int i=0; i<meas_list.size(); ++i)
                meas_always.add_term(measurements.get(meas_list[i]));
        }
        
        MPO<Matrix, grp> t_mpo = make_mpo(lat->size(), H);
        MPO<Matrix, grp> t_mpoc = t_mpo;
        if (parms["use_compressed"] > 0)
            t_mpoc.compress(1e-12);
        
        
        MPS<Matrix, grp> initial_mps;

        
        if (cur_mps.length() > 0 && cur_mps.length() != lat->size())
        {
            maquis::cout << "*** Starting grainings ***" << std::endl;
            
            boost::shared_ptr<mps_initializer<Matrix, grp> > initializer = boost::shared_ptr<mps_initializer<Matrix, grp> > (new empty_mps_init<Matrix, grp>());
            initial_mps = MPS<Matrix, grp>(lat->size(), 1, phys, initc, *initializer);
                        
            int oldL = old_model["Ndiscr"] * old_model["L"];
            std::vector<MPO<Matrix, grp> > mpo_mix(oldL+1, MPO<Matrix, grp>(0));
            double r = model["Ndiscr"]/old_model["Ndiscr"];
            for (int i=0; i<=oldL; ++i)
                mpo_mix[i] = mixed_mpo(model, r*i, old_model, oldL-i);
            
//            maquis::cout << "Old MPS:" << std::endl << initial_mps.description() << std::endl;
            if (cur_mps.length() < initial_mps.length())
                multigrid::extension_optim(parms, cur_mps, initial_mps, mpo_mix);
            else if (cur_mps.length() > initial_mps.length())
                multigrid::restriction(cur_mps, initial_mps);
//            maquis::cout << "New MPS:" << std::endl << initial_mps.description();

            std::vector<double> energies, entropies;            
            entropies = calculate_bond_entropies(initial_mps);
                        
            {
                storage::archive ar(rfile, "w");
                
                std::ostringstream oss;
                
                oss.str("");
                oss << "/simulation/iteration/graining/" << graining << "/parameters";
                ar[oss.str()] << parms;
                ar[oss.str()] << model;
                
                oss.str("");
                oss << "/simulation/iteration/graining/" << graining << "/results";
                
                ar[oss.str()] << storage::log;
                ar[oss.str()+"/Iteration Entropies/mean/value"] << entropies;
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
                                         parms["init_bond_dimension"],
                                         phys, initc,
                                         *initializer_factory<Matrix>(parms));
        }
        
        cur_mps = initial_mps;
        mpo = t_mpo;
        mpoc = t_mpoc;
        old_model = model;
        
#ifndef MEASURE_ONLY
        
        storage::setup(parms);
        
        bool early_exit = false;
        {   
            maquis::cout << "*** Starting optimization ***" << std::endl;
            ss_optimize<Matrix, grp, storage::disk> optimizer(initial_mps,
                                                                      parms["use_compressed"] == 0 ? mpo : mpoc,
                                                                      parms);
            
            for ( ; sweep < parms["nsweeps"]; ++sweep) {
                gettimeofday(&snow, NULL);
                
                optimizer.sweep(sweep);
                storage::disk::sync();
                
                cur_mps = optimizer.get_current_mps();
                
                std::vector<double> energies, entropies, renyi2;
                std::vector<std::size_t> truncations;

                entropies = calculate_bond_entropies(cur_mps);
                
                gettimeofday(&sthen, NULL);
                double elapsed = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
                maquis::cout << "Sweep done after " << elapsed << " seconds." << std::endl;
                
                {
                    storage::archive ar(rfile, "w");
                    
                    std::ostringstream oss;
                    
                    oss.str("");
                    oss << "/simulation/iteration/graining/" << graining << "/parameters";
                    ar[oss.str()] << parms;
                    ar[oss.str()] << model;
                    
                    oss.str("");
                    oss << "/simulation/iteration/graining/" << graining << "/sweep" << sweep << "/results";
                    
                    ar[oss.str()] << storage::log;                    
                    ar[oss.str()+"/Iteration Entropies/mean/value"] << entropies;
                    ar[oss.str()+"/Runtime/mean/value"] << std::vector<double>(1, elapsed);
                }
                {
                    std::ostringstream oss;
                    oss << "/simulation/iteration/graining/" << graining << "/sweep" << sweep << "/results/";
                    if (meas_always.n_terms() > 0)
                        measure_on_mps(cur_mps, *lat, meas_always, rfile, oss.str());
                }
                
                if (!dns)
                {
                    storage::archive ar(chkpfile, "w");
                    
                    ar["/state"] << cur_mps;
                    ar["/status/sweep"] << sweep;
                    ar["/status/graining"] << graining;
                }
                
                gettimeofday(&then, NULL);
                elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);            
                int rs = parms["run_seconds"];
                if (rs > 0 && elapsed > rs) {
                    early_exit = true;
                    break;
                }
            }
            storage::disk::sync();
        }
        sweep = 0;
#endif

        ++graining;
    } while (graining < raw_parms["ngrainings"]);
    
#ifdef MEASURE_ONLY
    {
        MPS<Matrix, grp> mps = cur_mps;
        
        std::vector<double> energies, entropies, renyi2;
        std::vector<std::size_t> truncations;
        
        maquis::cout << "Measurements." << std::endl;
        measure_on_mps(mps, *lat, measurements, rfile);
        
        maquis::cout << "Calculating vN entropy." << std::endl;
        entropies = calculate_bond_entropies(mps);
        maquis::cout << "Calculating n=2 Renyi entropy." << std::endl;
        renyi2 = calculate_bond_renyi_entropies(mps, 2);
        
        {
            storage::archive ar(rfile, "w");
            if (entropies.size() > 0)
                ar["/spectrum/results/Entropy/mean/value"] << entropies;
            if (renyi2.size() > 0)
                ar["/spectrum/results/Renyi2/mean/value"] << renyi2;
        }
        
        double energy = maquis::real(expval(mps, mpoc));
        maquis::cout << "Energy before: " << maquis::real(expval(mps, mpo)) << std::endl;
        maquis::cout << "Energy after: " << maquis::real(expval(mps, mpoc)) << std::endl;
        {
            storage::archive ar(rfile, "w");
            ar["/spectrum/results/Energy/mean/value"] << std::vector<double>(1, energy);
        }
        
        if (raw_parms["calc_h2"] > 0) {
            MPO<Matrix, grp> mpo2 = square_mpo(mpo);
            mpo2.compress(1e-12);
            
            double energy2 = maquis::real(expval(mps, mpo2, true));
            
            maquis::cout << "Energy^2: " << energy2 << std::endl;
            maquis::cout << "Variance: " << energy2 - energy*energy << std::endl;
            
            {
                storage::archive ar(rfile, "w");
                ar["/spectrum/results/Energy^2/mean/value"] << std::vector<double>(1, energy2);
                ar["/spectrum/results/EnergyVariance/mean/value"] << std::vector<double>(1, energy2 - energy*energy);
            }
        }
    }
#endif
    
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
    
    maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
    
}
