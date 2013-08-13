/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::sim(DmrgParameters const & parms_, ModelParameters const & model_)
: parms(parms_)
, model(model_)
, init_sweep(0)
, init_site(-1)
, restore(false)
, chkpfile(parms["chkpfile"].str())
, rfile(parms["resultfile"].str())
, dns( (parms["donotsave"] != 0) )
, stop_callback(static_cast<double>(parms["run_seconds"]))
{ 
    maquis::cout << DMRG_VERSION_STRING << std::endl;
    storage::setup(parms);
    
    {
		boost::filesystem::path p(chkpfile);
		if (boost::filesystem::exists(p) && boost::filesystem::is_regular_file(p))
        {
            storage::archive ar_in(chkpfile);
            if (ar_in.is_group("/state") && ar_in.is_scalar("/status/sweep"))
            {
                ar_in["/status/sweep"] >> init_sweep;
                
                if (ar_in.is_data("/status/site") && ar_in.is_scalar("/status/site"))
                    ar_in["/status/site"] >> init_site;
                
                if (init_site == -1)
                    ++init_sweep;
                
                maquis::cout << "Restoring state." << std::endl;
                maquis::cout << "Will start again at site " << init_site << " in sweep " << init_sweep << std::endl;
                restore = true;
            } else {
                maquis::cout << "A fresh simulation will start." << std::endl;
            }
        }
    }
    
    dmrg_random::engine.seed(parms["seed"]);
    
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
void sim<Matrix, SymmGroup>::model_init(boost::optional<int> opt_sweep)
{
    model_parser<Matrix, SymmGroup>(parms["lattice_library"], parms["model_library"], model, lat, phys_model);
    initc = phys_model->initc(model);
    measurements = phys_model->measurements();
    if (opt_sweep)
        parse_overlaps(model, opt_sweep.get(), measurements);
    else
        parse_overlaps(model, init_sweep, measurements);

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
        
        // TODO: move to ts_optim init
        if (parms["optimization"] == "twosite") {
            Timer t("TS_MPO"); t.begin();
            throw std::runtime_error("compression should be moved inside ts_optim constructor.");
//            make_ts_cache_mpo(scratch_mpo, ts_cache_mpo, phys);
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
    }
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::mps_init()
{
    #ifdef AMBIENT_TRACKING
    ambient::overseer::log::region("sim::mps_init");
    #endif
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
        #ifdef AMBIENT_TRACKING
        for(int i = 0; i < mps.length(); ++i) __ambient_track_array(mps, i);
        #endif
    }
}


template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::~sim()
{
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::checkpoint_state(MPS<Matrix, SymmGroup> const& state, status_type const& status)
{
    if (!dns) {
        storage::archive ar(chkpfile, "w");
        ar["/state"]  << state;
        ar["/status"] << status;
    }
}

template <class Matrix, class SymmGroup>
std::string sim<Matrix, SymmGroup>::results_archive_path(status_type const& status) const
{
    std::ostringstream oss;
    oss.str("");
    oss << "/simulation/sweep" << status.at("sweep");
    return oss.str();
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::measure(std::string archive_path, Measurements<Matrix, SymmGroup> const& meas)
{
    maquis::cout << "Measurements." << std::endl;
    measure_on_mps(mps, *lat, meas, rfile, archive_path);
    
    // TODO: move into special measurement
    std::vector<double> entropies, renyi2;
    if (model["ENABLE_MEASURE[Entropy]"]) {
        maquis::cout << "Calculating vN entropy." << std::endl;
        entropies = calculate_bond_entropies(mps);
    }
    if (model["ENABLE_MEASURE[Renyi2]"]) {
        maquis::cout << "Calculating n=2 Renyi entropy." << std::endl;
        renyi2 = calculate_bond_renyi_entropies(mps, 2);
    }
    std::vector< std::vector<double> > * spectra;
    if (parms["entanglement_spectra"])
        spectra = new std::vector< std::vector<double> >();
    else
        spectra = NULL;

    {
        storage::archive ar(rfile, "w");
        if (entropies.size() > 0)
            ar[archive_path + "/Entropy/mean/value"] << entropies;
        if (renyi2.size() > 0)
            ar[archive_path + "/Renyi2/mean/value"] << renyi2;
        if (spectra != NULL)
            ar[archive_path + "/Entanglement Spectra/mean/value"] << *spectra;
    }
}

