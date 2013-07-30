/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::sim(DmrgParameters const & parms_, ModelParameters const & model_, bool fullinit)
: parms(parms_)
, model(model_)
, sweep(0)
, site(-1)
, chkpfile(parms["chkpfile"].str())
, rfile(parms["resultfile"].str())
, dns( (parms["donotsave"] != 0) )
{ 
    maquis::cout << DMRG_VERSION_STRING << std::endl;
    storage::setup(parms);
    
    DCOLLECTOR_GROUP(gemm_collector, "init")
    DCOLLECTOR_GROUP(svd_collector, "init")
    
    gettimeofday(&now, NULL);

    restore = false;
    {
		boost::filesystem::path p(chkpfile);
		if (boost::filesystem::exists(p) && boost::filesystem::is_regular_file(p))
        {
            storage::archive ar_in(chkpfile);
            if (ar_in.is_group("/state") && ar_in.is_scalar("/status/sweep"))
            {
                ar_in["/status/sweep"] >> sweep;
                
                if (ar_in.is_data("/status/site") && ar_in.is_scalar("/status/site"))
                    ar_in["/status/site"] >> site;
                
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
    srand48(parms["seed"]);
#endif
    dmrg_random::engine.seed(parms["seed"]);
    
    if (fullinit) {
        model_init();
        mps_init();
    }

    {
        storage::archive ar(rfile, "w");
        
        ar["/parameters"] << parms;
        ar["/parameters"] << model;
        ar["/version"] << DMRG_VERSION_STRING;
    }
    if (!dns)
    {
        storage::archive ar(chkpfile, "w");
        
        ar["/parameters"] << parms;
        ar["/parameters"] << model;
        ar["/version"] << DMRG_VERSION_STRING;
    }
    
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::model_init()
{
    model_parser<Matrix, SymmGroup>(parms["lattice_library"], parms["model_library"], model, lat, phys_model);
    initc = phys_model->initc(model);
    measurements = phys_model->measurements();
    parse_overlaps(model, sweep, measurements);

    /*
     maquis::cout << "initc: " << initc << std::endl;
     maquis::cout << "phys:" << std::endl << phys << std::endl;
     maquis::cout << measurements << std::endl;
     maquis::cout << "Hamiltonian:" << std::endl << H << std::endl;
     */

    if (model["MODEL"] == std::string("quantum_chemistry"))
    {  
        typedef typename alps::numeric::associated_one_matrix<Matrix>::type MPOMatrix;
        MPO<MPOMatrix, SymmGroup> scratch_mpo;

        Hamiltonian<MPOMatrix, SymmGroup> Hloc = phys_model->H_chem();
        phys = Hloc.get_phys();

        make_compressed_mpo(scratch_mpo, 1e-12, lat->size(), Hloc);
        Timer t("DENSE_MPO conversion"); t.begin();
        compressor<MPOMatrix, SymmGroup>::convert_to_dense_matrix(scratch_mpo, mpoc);
        t.end();

        if (parms["optimization"] == "twosite") {
            Timer t("TS_MPO"); t.begin();
            make_ts_cache_mpo(scratch_mpo, ts_cache_mpo, phys);
            t.end();
        }
    }
    else
    {  
        H = phys_model->H();
        phys = H.get_phys();

        mpo = make_mpo(lat->size(), H); 
        mpoc = mpo;

        if (parms["use_compressed"] > 0)
            mpoc.compress(1e-12);

        if (parms["optimization"] == "twosite")
            make_ts_cache_mpo(mpoc, ts_cache_mpo, phys);
    }
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::mps_init()
{
    assert(lat.get() != NULL);
    
    if (restore) {
        mps = MPS<Matrix, SymmGroup>(lat->size());
        storage::archive ar_in(chkpfile);
        ar_in["/state"] >> mps;
    } else if (!parms["initfile"].empty()) {
        maquis::cout << "Loading init state from " << parms["initfile"] << std::endl;
        mps = MPS<Matrix, SymmGroup>(lat->size());
        storage::archive ar_in(parms["initfile"].str());
        ar_in["/state"] >> mps;
    } else {
        mps = MPS<Matrix, SymmGroup>(lat->size(),
                                     parms["init_bond_dimension"],
                                     phys, initc,
                                     *(phys_model->initializer(parms)));
    }
}


template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::~sim()
{
}

template <class Matrix, class SymmGroup>
int sim<Matrix, SymmGroup>::advance (int nsteps, double time_limit)
{
    int site = -1;
    int ns = sweep + nsteps;
    for (; sweep < ns && site < 0; ++sweep) // do_sweep requires the correct sweep number!
    {
        parms.set("sweep", sweep);
        site = do_sweep(time_limit);
    }
    // sweep = ns !
    --sweep;
    return site;
}

template <class Matrix, class SymmGroup>
bool sim<Matrix, SymmGroup>::run ()
{
    bool early_exit = false;
    
    int measure_each = parms["measure_each"];
    int chkp_each    = parms["chkp_each"];
    int update_each  = parms["update_each"];
    
    int nsteps = parms["nsweeps"];
    if (measure_each > -1)
        nsteps = std::min(nsteps, measure_each);
    if (chkp_each > -1)
        nsteps = std::min(nsteps, chkp_each);
    if (update_each > -1)
        nsteps = std::min(nsteps, update_each);
    
    while (sweep < parms["nsweeps"]) {
        DCOLLECTOR_GROUP(gemm_collector, "sweep"+boost::lexical_cast<std::string>(sweep))
        DCOLLECTOR_GROUP(svd_collector, "sweep"+boost::lexical_cast<std::string>(sweep))
        gettimeofday(&snow, NULL);
        
        int rs = parms["run_seconds"];
        
        gettimeofday(&then, NULL);
        double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);            
        
        int sweep_after = sweep + nsteps - 1;
        site = advance(std::min(parms["nsweeps"]-sweep, nsteps), rs > 0 ? rs-elapsed : -1);
        early_exit = (site >= 0);
        assert(sweep == sweep_after);
        
        gettimeofday(&sthen, NULL);
        double elapsed_sweep = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
        maquis::cout << "Sweep " << sweep << " done " << nsteps << " steps after " << elapsed_sweep << " seconds." << std::endl;
        
        
        gettimeofday(&then, NULL);
        elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);            
        if (rs > 0 && elapsed > rs)
            early_exit = true;
        
        
        if (early_exit ||
            (   (sweep+1) % measure_each == 0 
             || (sweep+1) == parms["nsweeps"] ))
        {
            double elapsed_measure = 0.;
            if (!early_exit) {
                do_sweep_measure();
                
                gettimeofday(&sthen, NULL);
                elapsed_measure = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
                maquis::cout << "Sweep measure done after " << elapsed_measure << " seconds." << std::endl;
                
                gettimeofday(&then, NULL);
                elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
                if (rs > 0 && elapsed > rs)
                    early_exit = true;
            }
            
            // Write results
            {
                storage::archive ar(rfile, "w");
                ar[sweep_archive_path() + "/parameters"] << parms;
                ar[sweep_archive_path() + "/parameters"] << model;
                ar[sweep_archive_path() + "/results"] << storage::log;
                ar[sweep_archive_path() + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
            }
        }
        
        // Write checkpoint
        if (!dns && (early_exit ||
                     (sweep+1) % chkp_each == 0 ||
                     (sweep+1) == parms["nsweeps"]))
        {
            storage::archive ar(chkpfile, "w");
            ar["/state"] << mps;
            ar["/status/sweep"] << sweep;
            ar["/status/site"] << site;
        }
        
        if (early_exit) return true;
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
void sim<Matrix, SymmGroup>::do_sweep_measure()
{
    if (!parms["always_measure"].empty())
        meas_always = measurements.sublist( parms.get<std::vector<std::string> >("always_measure") );
    
    std::vector< std::vector<double> > * spectra;
    if (parms["entanglement_spectra"])
        spectra = new std::vector< std::vector<double> >();
    else
        spectra = NULL;
    
    maquis::cout << "Calculating vN entropy." << std::endl;
    std::vector<double> entropies = calculate_bond_entropies(mps, spectra);
    
    {
        storage::archive ar(rfile, "w");

        ar[sweep_archive_path() + "/results/Iteration Entropies/mean/value"] << entropies;
        if (spectra != NULL)
            ar[sweep_archive_path() + "/results/Iteration Entanglement Spectra/mean/value"] << *spectra;
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
        storage::archive ar(rfile, "w");
        
        ar["/parameters"] << parms;
        ar["/parameters"] << model;
        ar["/version"] << DMRG_VERSION_STRING;
    }
    
    maquis::cout << "Measurements." << std::endl;
    measure_on_mps(mps, *lat, measurements, rfile);
    
    maquis::cout << "Calculating vN entropy." << std::endl;
    std::vector<double> entropies = calculate_bond_entropies(mps);
    maquis::cout << "Calculating n=2 Renyi entropy." << std::endl;
    std::vector<double> renyi2 = calculate_bond_renyi_entropies(mps, 2);
    
    {
        storage::archive ar(rfile, "w");
        if (entropies.size() > 0)
            ar["/spectrum/results/Entropy/mean/value"] << entropies;
        if (renyi2.size() > 0)
            ar["/spectrum/results/Renyi2/mean/value"] << renyi2;
    }
    
    double energy = maquis::real(expval(mps, mpoc));
    // MD: removed redundant energy calculation
    // maquis::cout << "Energy before: " << maquis::real(expval(mps, mpo)) << std::endl;
    maquis::cout << "Energy: " << maquis::real(expval(mps, mpoc)) << std::endl;
    {
        storage::archive ar(rfile, "w");
        ar["/spectrum/results/Energy/mean/value"] << std::vector<double>(1, energy);
    }
    
    if (parms["calc_h2"] > 0) {
        MPO<Matrix, SymmGroup> mpo2 = square_mpo(mpoc);
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

