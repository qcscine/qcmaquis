/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include "dmrg/models/measure_on_mps.h"

#include <boost/algorithm/string.hpp>

template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::sim(DmrgParameters const & parms_, ModelParameters const & model_)
: parms(parms_)
, model(model_)
, init_sweep(0)
, init_site(-1)
, restore(false)
, chkpfile(boost::trim_right_copy_if(parms["chkpfile"].str(), boost::is_any_of("/ ")))
, rfile(parms["resultfile"].str())
, dns( (parms["donotsave"] != 0) )
, stop_callback(static_cast<double>(parms["run_seconds"]))
{ 
    maquis::cout << DMRG_VERSION_STRING << std::endl;
    storage::setup(parms);
    
    {
        boost::filesystem::path p(chkpfile);
        if (boost::filesystem::exists(p) && boost::filesystem::exists(p / "mps0.h5"))
        {
            storage::archive ar_in(chkpfile+"/props.h5");
            if (ar_in.is_scalar("/status/sweep"))
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
        if (!boost::filesystem::exists(chkpfile))
            boost::filesystem::create_directory(chkpfile);
        storage::archive ar(chkpfile+"/props.h5", "w");
        
        ar["/parameters"] << parms;
        ar["/parameters"] << model;
        ar["/version"] << DMRG_VERSION_STRING;
    }
    
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::model_init(boost::optional<int> opt_sweep)
{
    lat = Lattice(parms, model);
    phys_model = Model<Matrix, SymmGroup>(lat, parms, model);
    measurements = phys_model.measurements();
    if (opt_sweep)
        parse_overlaps(model, opt_sweep.get(), measurements);
    else
        parse_overlaps(model, init_sweep, measurements);

    if (model["MODEL"] == std::string("quantum_chemistry") && (parms.get<int>("use_compressed") > 0))
    {  
        throw std::runtime_error("chem compression has been disabled");
        /*
        typedef typename alps::numeric::associated_one_matrix<Matrix>::type MPOMatrix;
        MPO<MPOMatrix, SymmGroup> scratch_mpo;

        Hamiltonian<MPOMatrix, SymmGroup> Hloc = phys_model->H_chem();
        phys = Hloc.get_phys();

        make_compressed_mpo(scratch_mpo, 1e-12, lat->size(), Hloc);
        Timer t("DENSE_MPO conversion"); t.begin();
        compressor<MPOMatrix, SymmGroup>::convert_to_dense_matrix(scratch_mpo, mpoc);
        t.end();
        */
    }
    else
    {  
        mpo = make_mpo(lat, phys_model, model);
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
    
    if (restore) {
        load(chkpfile, mps);
    } else if (!parms["initfile"].empty()) {
        maquis::cout << "Loading init state from " << parms["initfile"] << std::endl;
        load(parms["initfile"].str(), mps);
    } else {
        mps = MPS<Matrix, SymmGroup>(lat.size(), *(phys_model.initializer(lat, parms, model)));
        #ifdef AMBIENT_TRACKING
        for(int i = 0; i < mps.length(); ++i) ambient_track_array(mps, i);
        #endif
    }
    assert(mps.length() == lat.size());
}


template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::~sim()
{
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::checkpoint_simulation(MPS<Matrix, SymmGroup> const& state, status_type const& status)
{
    if (!dns) {
        /// create chkp dir
        #ifdef USE_AMBIENT
        if (ambient::master() && !boost::filesystem::exists(chkpfile))
            boost::filesystem::create_directory(chkpfile);
        #else
        if (!boost::filesystem::exists(chkpfile))
            boost::filesystem::create_directory(chkpfile);
        #endif
        /// save state to chkp dir
        save(chkpfile, state);
        
        /// save status
        #ifdef USE_AMBIENT
        if(!ambient::master()) return;
        #endif
        storage::archive ar(chkpfile+"/props.h5", "w");
        ar["/status"] << status;
    }
}

template <class Matrix, class SymmGroup>
std::string sim<Matrix, SymmGroup>::results_archive_path(status_type const& status) const
{
    std::ostringstream oss;
    oss.str("");
#if defined(__xlC__)
    typename status_type::const_iterator match = status.find("sweep");
    oss << "/simulation/sweep" << match->second;
#else
    oss << "/simulation/sweep" << status.at("sweep");
#endif
    return oss.str();
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::measure(std::string archive_path, Measurements<Matrix, SymmGroup> const& meas)
{
    maquis::cout << "Measurements." << std::endl;
    measure_on_mps(mps, lat, meas, rfile, archive_path);
    
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

