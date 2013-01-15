/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

namespace detail {
    
    template <class Matrix, class SymmGroup>
    struct call_linear_init {
        static mps_initializer<Matrix, SymmGroup> * call()
        {
            throw std::runtime_error("Linear MPS init is available only for U1 symmetry group.");
            return new default_mps_init<Matrix, SymmGroup>();
        }
    };
    template <class Matrix>
    struct call_linear_init<Matrix, U1> {
        static mps_initializer<Matrix, U1> * call()
        {
            return new linear_mps_init<Matrix>();
        }
    };

    template <class Matrix, class SymmGroup>
    struct call_hf_init {
        static mps_initializer<Matrix, SymmGroup> * call(BaseParameters & params)
        {
            throw std::runtime_error("HF MPS init is available only for TwoU1 symmetry group.");
            return new default_mps_init<Matrix, SymmGroup>();
        }
    };
    template <class Matrix>
    struct call_hf_init<Matrix, TwoU1> {
        static mps_initializer<Matrix, TwoU1> * call(BaseParameters & params)
        {
            return new hf_mps_init<Matrix>(params);
        }
    }; 
    
}

template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::sim(DmrgParameters const & parms_, ModelParameters const & model_, bool fullinit)
: parms(parms_)
, model(model_)
, sweep(0)
, site(-1)
, chkpfile(parms.get<std::string>("chkpfile"))
, rfile(parms.get<std::string>("resultfile"))
, ssm(parms.get<std::string>("storagedir"))
, dns( (parms.get<int>("donotsave") != 0) )
{ 
    maquis::cout << DMRG_VERSION_STRING << std::endl;
    
    DCOLLECTOR_GROUP(gemm_collector, "init")
    DCOLLECTOR_GROUP(svd_collector, "init")
    
    gettimeofday(&now, NULL);

#ifdef USE_GPU
    cublasInit();
#endif
    restore = false;
    {
		boost::filesystem::path p(chkpfile);
		if (boost::filesystem::exists(p) && boost::filesystem::is_regular_file(p))
        {
            alps::hdf5::archive h5ar_in(chkpfile);
            if (h5ar_in.is_group("/state") && h5ar_in.is_scalar("/status/sweep"))
            {
                h5ar_in >> alps::make_pvp("/status/sweep", sweep);
                
                if (h5ar_in.is_data("/status/site") && h5ar_in.is_scalar("/status/site"))
                    h5ar_in >> alps::make_pvp("/status/site", site);
                
                if (site == -1)
                    ++sweep;
                    
                maquis::cout << "Restoring state." << std::endl;
                maquis::cout << "Will start again at site " << site << " in sweep " << sweep << std::endl;
                restore = true;
            } else {
                maquis::cout << "Invalid checkpoint, overwriting." << std::endl; 
            }
        }
    }
    
#ifndef WIN32
    // drand48 does not seem to be used anymore anyway
    srand48(parms.get<int>("seed"));
#endif
    dmrg_random::engine.seed(parms.get<int>("seed"));
    
    if (fullinit) {
        model_init();
        mps_init();
    }

    {
        alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
        h5ar << alps::make_pvp("/version", DMRG_VERSION_STRING);
    }
    if (!dns)
    {
        alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
        h5ar << alps::make_pvp("/version", DMRG_VERSION_STRING);
    }
    
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::model_init()
{
    model_parser<Matrix, SymmGroup>(parms.get<std::string>("lattice_library"), parms.get<std::string>("model_library"), model, lat, phys_model);
    initc = phys_model->initc(model);
    H = phys_model->H();
    measurements = phys_model->measurements();
    parse_overlaps(model, sweep, measurements);
    phys = H.get_phys();
    
    /*
     maquis::cout << "initc: " << initc << std::endl;
     maquis::cout << "phys:" << std::endl << phys << std::endl;
     maquis::cout << measurements << std::endl;
     maquis::cout << "Hamiltonian:" << std::endl << H << std::endl;
     */
    
    mpo = make_mpo(lat->size(), H);
    mpoc = mpo;
    if (parms.get<int>("use_compressed") > 0)
        mpoc.compress(1e-12);

    if (parms.get<std::string>("optimization") == "twosite") {
        Timer t("TS_MPO"); t.begin();
        ts_cache_mpo = make_ts_cache_mpo(mpoc, phys);
        t.end();
    }
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::mps_init()
{
    assert(lat.get() != NULL);
    
    if (restore) {
        mps = MPS<Matrix, SymmGroup>(lat->size());
        alps::hdf5::archive h5ar_in(chkpfile);
        h5ar_in >> alps::make_pvp("/state", mps);
    } else if (parms.get<std::string>("initfile").size() > 0) {
        maquis::cout << "Loading init state from " << parms.get<std::string>("initfile") << std::endl;
        mps = MPS<Matrix, SymmGroup>(lat->size());
        alps::hdf5::archive h5ar_in(parms.get<std::string>("initfile"));
        h5ar_in >> alps::make_pvp("/state", mps);
    } else {
        mps = MPS<Matrix, SymmGroup>(lat->size(),
                                     parms.get<std::size_t>("init_bond_dimension"),
                                     phys, initc,
                                     *initializer_factory(parms));
    }
}


template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::~sim()
{
    ssm.sync();
#ifdef USE_GPU
    cublasShutdown();
#endif
}

template <class Matrix, class SymmGroup>
int sim<Matrix, SymmGroup>::advance (Logger& iteration_log, int nsteps, double time_limit)
{
    int site = -1;
    int ns = sweep + nsteps;
    for (; sweep < ns && site < 0; ++sweep) // do_sweep requires the correct sweep number!
    {
        parms.set("sweep", sweep);
        site = do_sweep(iteration_log, time_limit);
    }
    // sweep = ns !
    --sweep;
    return site;
}

template <class Matrix, class SymmGroup>
bool sim<Matrix, SymmGroup>::run ()
{
    bool early_exit = false;
    
    int nsteps = parms.get<int>("nsweeps");
    if (parms.get<int>("measure_each") > -1)
        nsteps = std::min(nsteps, parms.get<int>("measure_each"));
    if (parms.get<int>("measure_each") > -1)
        nsteps = std::min(nsteps, parms.get<int>("chkp_each"));
    if (parms.get<int>("update_each") > -1)
        nsteps = std::min(nsteps, parms.get<int>("update_each"));
    
    while (sweep < parms.get<int>("nsweeps")) {
        DCOLLECTOR_GROUP(gemm_collector, "sweep"+boost::lexical_cast<std::string>(sweep))
        DCOLLECTOR_GROUP(svd_collector, "sweep"+boost::lexical_cast<std::string>(sweep))
        gettimeofday(&snow, NULL);
        
        Logger iteration_log;
        int rs = parms.get<int>("run_seconds");
        
        gettimeofday(&then, NULL);
        double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);            
        
        int sweep_after = sweep + nsteps - 1;
        site = advance(iteration_log, std::min(parms.get<int>("nsweeps")-sweep, nsteps), rs > 0 ? rs-elapsed : -1);
        early_exit = (site >= 0);
        assert(sweep == sweep_after);
        
        gettimeofday(&sthen, NULL);
        double elapsed_sweep = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
        maquis::cout << "Sweep " << sweep << " done " << nsteps << " steps after " << elapsed_sweep << " seconds." << std::endl;
        
        
        gettimeofday(&then, NULL);
        elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);            
        if (rs > 0 && elapsed > rs)
            early_exit = true;
        
        
        if (!early_exit &&
            (   (sweep+1) % parms.get<int>("measure_each") == 0 
             || (sweep+1) == parms.get<int>("nsweeps") ))
        {
            
            do_sweep_measure(iteration_log);
            
            gettimeofday(&sthen, NULL);
            double elapsed_measure = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
            maquis::cout << "Sweep measure done after " << elapsed_measure << " seconds." << std::endl;
            
            gettimeofday(&then, NULL);
            elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);            
            if (rs > 0 && elapsed > rs)
                early_exit = true;
            
            
            // Write results
            {
                alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                
                h5ar << alps::make_pvp(sweep_archive_path() + "/parameters",
                                       parms);
                h5ar << alps::make_pvp(sweep_archive_path() + "/parameters",
                                       model);
                
                
                h5ar << alps::make_pvp(sweep_archive_path() + "/results",
                                       iteration_log);
                
                h5ar << alps::make_pvp(sweep_archive_path() + "/results/Runtime/mean/value",
                                       std::vector<double>(1, elapsed_sweep + elapsed_measure));                
            }
        }
        
        
        
        // Write checkpoint
        if (!dns && (early_exit ||
                     (sweep+1) % parms.get<int>("chkp_each") == 0 ||
                     (sweep+1) == parms.get<int>("nsweeps")))
        {
            alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            
            h5ar << alps::make_pvp("/state", mps);
            h5ar << alps::make_pvp("/status/sweep", sweep);
            h5ar << alps::make_pvp("/status/site", site);
        }
        
        
        if (early_exit)
            return true;
        
        ++sweep;
    }
    return false;
}

template <class Matrix, class SymmGroup>
std::string sim<Matrix, SymmGroup>::sweep_archive_path ()
{
    std::ostringstream oss;
    oss.str("");
    oss << "/simulation/sweep" << sweep;
    return oss.str();
}


template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::do_sweep_measure (Logger&)
{
    if (!parms.get<std::string>("always_measure").empty()) {
        meas_always.clear();
        meas_always.set_identity(measurements.get_identity());
        std::vector<std::string> meas_list = parms.get<std::vector<std::string> >("always_measure");
        for (int i=0; i<meas_list.size(); ++i) {
            try {
                meas_always.add_term(measurements.get(meas_list[i]));
            } catch (std::runtime_error & e) {
                maquis::cout << "WARNING: " << e.what() << std::endl;
            }
        }
    }
    
    maquis::cout << "Calculating vN entropy." << std::endl;
    std::vector<double> entropies = calculate_bond_entropies(mps);
    
    {
        alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);

        h5ar << alps::make_pvp(sweep_archive_path() + "/results/Iteration Entropies/mean/value",
                               entropies);
    }
    
    {
        if (meas_always.n_terms() > 0)
            measure_on_mps(mps, *lat, meas_always, rfile, sweep_archive_path() + "/results/");
        
    }
}


template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::measure ()
{
    {
        alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
        h5ar << alps::make_pvp("/version", DMRG_VERSION_STRING);
    }
    
    maquis::cout << "Measurements." << std::endl;
    measure_on_mps(mps, *lat, measurements, rfile);
    
    maquis::cout << "Calculating vN entropy." << std::endl;
    std::vector<double> entropies = calculate_bond_entropies(mps);
    maquis::cout << "Calculating n=2 Renyi entropy." << std::endl;
    std::vector<double> renyi2 = calculate_bond_renyi_entropies(mps, 2);
    
    {
        alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        if (entropies.size() > 0)
            h5ar << alps::make_pvp("/spectrum/results/Entropy/mean/value", entropies);
        if (renyi2.size() > 0)
            h5ar << alps::make_pvp("/spectrum/results/Renyi2/mean/value", renyi2);
    }
    
    double energy = maquis::real(expval(mps, mpoc));
    maquis::cout << "Energy before: " << maquis::real(expval(mps, mpo)) << std::endl;
    maquis::cout << "Energy after: " << maquis::real(expval(mps, mpoc)) << std::endl;
    {
        alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        h5ar << alps::make_pvp("/spectrum/results/Energy/mean/value", std::vector<double>(1, energy));
    }
    
    if (parms.get<int>("calc_h2") > 0) {
        MPO<Matrix, SymmGroup> mpo2 = square_mpo(mpo);
        mpo2.compress(1e-12);
        
        double energy2 = maquis::real(expval(mps, mpo2, true));
        
        maquis::cout << "Energy^2: " << energy2 << std::endl;
        maquis::cout << "Variance: " << energy2 - energy*energy << std::endl;
        
        {
            alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            h5ar << alps::make_pvp("/spectrum/results/Energy^2/mean/value", std::vector<double>(1, energy2));
            h5ar << alps::make_pvp("/spectrum/results/EnergyVariance/mean/value",
                                   std::vector<double>(1, energy2 - energy*energy));
        }
    }
}



template <class Matrix, class SymmGroup>
mps_initializer<Matrix, SymmGroup> * sim<Matrix, SymmGroup>::initializer_factory(BaseParameters & params)
{
    if (params.get<std::string>("init_state") == "default")
        return new default_mps_init<Matrix, SymmGroup>();
    else if (params.get<std::string>("init_state") == "linear")
        return detail::call_linear_init<Matrix, SymmGroup>::call();
    else if (params.get<std::string>("init_state") == "const")
        return new const_mps_init<Matrix, SymmGroup>();
    else if (params.get<std::string>("init_state") == "thin")
        return new thin_mps_init<Matrix, SymmGroup>();
    else if (params.get<std::string>("init_state") == "thin_const")
        return new thin_const_mps_init<Matrix, SymmGroup>();
    else if (params.get<std::string>("init_state") == "hf")
        return detail::call_hf_init<Matrix, SymmGroup>::call(params);
    else {
        throw std::runtime_error("Don't know this initial state.");
        return NULL;
    }
}

