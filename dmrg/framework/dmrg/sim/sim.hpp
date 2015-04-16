/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#include <boost/algorithm/string.hpp>
#include "dmrg/version.h"

template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::sim(DmrgParameters const & parms_)
: parms(parms_)
, init_sweep(0)
, init_site(-1)
, restore(false)
, dns( (parms["donotsave"] != 0) )
, chkpfile(boost::trim_right_copy_if(parms["chkpfile"].str(), boost::is_any_of("/ ")))
, rfile(parms["resultfile"].str())
, stop_callback(static_cast<double>(parms["run_seconds"]))
{ 
    maquis::cout << DMRG_VERSION_STRING << std::endl;
    storage::setup(parms);
    dmrg_random::engine.seed(parms["seed"]);
    
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
    
    
    {
        storage::archive ar(rfile, "w");
        
        ar["/parameters"] << parms;
        ar["/version"] << DMRG_VERSION_STRING;
    }
    if (!dns)
    {
        if (!boost::filesystem::exists(chkpfile))
            boost::filesystem::create_directory(chkpfile);
        storage::archive ar(chkpfile+"/props.h5", "w");
        
        ar["/parameters"] << parms;
        ar["/version"] << DMRG_VERSION_STRING;
    }
    
    
    /// Model initialization
    lat = Lattice(parms);
    model = Model<Matrix, SymmGroup>(lat, parms);
    mpo = make_mpo(lat, model);
    all_measurements = model.measurements();
    all_measurements << overlap_measurements<Matrix, SymmGroup>(parms);
    
    
    /// MPS initialization
    if (restore) {
        load(chkpfile, mps);
    } else if (!parms["initfile"].empty()) {
        maquis::cout << "Loading init state from " << parms["initfile"] << std::endl;
        load(parms["initfile"].str(), mps);
    } else {
        mps = MPS<Matrix, SymmGroup>(lat.size(), *(model.initializer(lat, parms)));
    }
    
    assert(mps.length() == lat.size());
    maquis::cout << "MPS initialization has finished...\n"; // MPS restored now
    
}

template <class Matrix, class SymmGroup>
typename sim<Matrix, SymmGroup>::measurements_type
sim<Matrix, SymmGroup>::iteration_measurements(int sweep)
{
    measurements_type mymeas(all_measurements);
    mymeas << overlap_measurements<Matrix, SymmGroup>(parms, sweep);
    
    measurements_type sweep_measurements;
    if (!parms["ALWAYS_MEASURE"].empty())
        sweep_measurements = meas_sublist(mymeas, parms["ALWAYS_MEASURE"]);
    
    return sweep_measurements;
}


template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::~sim()
{
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::checkpoint_simulation(MPS<Matrix, SymmGroup> const& state, status_type const& status)
{
    if (!dns) {
        /// save state to chkp dir
        save(chkpfile, state);
        
        /// save status
        if(!parallel::master()) return;
        storage::archive ar(chkpfile+"/props.h5", "w");
        ar["/status"] << status;
    }
}

template <class Matrix, class SymmGroup>
std::string sim<Matrix, SymmGroup>::results_archive_path(status_type const& status) const
{
    std::ostringstream oss;
    oss.str("");
#if defined(__xlC__) || defined(__FCC_VERSION)
    typename status_type::const_iterator match = status.find("sweep");
    oss << "/spectrum/iteration/" << match->second;
#else
    oss << "/spectrum/iteration/" << status.at("sweep");
#endif
    return oss.str();
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::measure(std::string archive_path, measurements_type & meas)
{
    std::for_each(meas.begin(), meas.end(), measure_and_save<Matrix, SymmGroup>(rfile, archive_path, mps));

    // TODO: move into special measurement
    std::vector<int> * measure_es_where = NULL;
    entanglement_spectrum_type * spectra = NULL;
    if (parms.defined("entanglement_spectra")) {
        spectra = new entanglement_spectrum_type();
        measure_es_where = new std::vector<int>();
        *measure_es_where = parms.template get<std::vector<int> >("entanglement_spectra");
    }
    std::vector<double> entropies, renyi2;
    if (parms["MEASURE[Entropy]"]) {
        maquis::cout << "Calculating vN entropy." << std::endl;
        entropies = calculate_bond_entropies(mps);
    }
    if (parms["MEASURE[Renyi2]"]) {
        maquis::cout << "Calculating n=2 Renyi entropy." << std::endl;
        renyi2 = calculate_bond_renyi_entropies(mps, 2, measure_es_where, spectra);
    }

    {
        storage::archive ar(rfile, "w");
        if (entropies.size() > 0)
            ar[archive_path + "Entropy/mean/value"] << entropies;
        if (renyi2.size() > 0)
            ar[archive_path + "Renyi2/mean/value"] << renyi2;
        if (spectra != NULL)
            ar[archive_path + "Entanglement Spectra/mean/value"] << *spectra;
    }
}

