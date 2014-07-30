/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef ALPS_MPS_OPTIM_RUN_EIGENSTATE_SIM_HPP
#define ALPS_MPS_OPTIM_RUN_EIGENSTATE_SIM_HPP

#include "dmrg_sim.hpp"
#include "dmrg/models/measurements.h"
#include "mps_eigenstate_sims/simulation.hpp"

#include <algorithm>
#include <iterator>
#include <sstream>

#include <alps/parser/xmlstream.h>

class save_to_xml {
public:
    save_to_xml(alps::oxstream& o) : out(o) { }
    
    template <class Matrix, class SymmGroup>
    void operator()(measurement<Matrix, SymmGroup> const& m)
    {
        m.write_xml(out);
    }
    
private:
    alps::oxstream& out;
};

template <class Matrix, class SymmGroup>
void run_eigenstate_sim(BaseParameters parms, bool write_xml, run_type rt)
{
    int neigen = parms["NUMBER_EIGENVALUES"];
    std::vector<std::string> checkpoints(neigen);
    for (int eig=0; eig<neigen; ++eig) {
        DmrgParameters myparms(parms);
        std::string resfile = parms["resultfile"].str();
        std::string ckpfile = parms["chkpfile"].str();
        if (neigen > 1) {
            std::string eig_suffix = std::string(".") + boost::lexical_cast<std::string>(eig);
            resfile = boost::replace_last_copy(parms["resultfile"].str(), ".h5",   eig_suffix+".h5");
            ckpfile = boost::replace_last_copy(parms["chkpfile"].str(),   ".chkp", eig_suffix+".chkp");
            myparms.set("resultfile", resfile);
            myparms.set("chkpfile",   ckpfile);
            
            std::stringstream ortho_states;
            std::copy(checkpoints.begin(), checkpoints.begin()+eig, std::ostream_iterator<std::string>(ortho_states, ","));
            myparms.set("n_ortho_states", eig);
            myparms.set("ortho_states", ortho_states.str());
        }
        
        /// run optimization
        if (rt == optim_only || rt == optim_and_measure) {
            dmrg_sim<Matrix, SymmGroup> sim(myparms);
            sim.run();
        }
        
        checkpoints[eig] = ckpfile;
    }
    
    
    /// perform measurements
    if (rt == measure_only || rt == optim_and_measure) {
        /// Build model
        Lattice lattice(parms);
        Model<Matrix, SymmGroup> model(lattice, parms);
        MPO<Matrix, SymmGroup> mpo = make_mpo(lattice, model, parms);

        /// Get measurements
        typedef typename Model<Matrix, SymmGroup>::measurements_type measurements_type;
        measurements_type measurements = model.measurements();
        { // overlap measurements
            measurements_type m = (parms.defined("nsweeps")) ?  overlap_measurements<Matrix, SymmGroup>(parms, parms["nsweeps"]-1) : overlap_measurements<Matrix, SymmGroup>(parms);
            measurements.insert(measurements.end(), m.begin(), m.end());
        }
    
        std::string rfile = parms["resultfile"].str();
        {
            storage::archive ar(rfile, "w");
            ar["/parameters"] << parms;
        }
    
        alps::oxstream out(boost::replace_last_copy(rfile, ".h5", ".xml"));
        out << alps::header("UTF-8") << alps::stylesheet(alps::xslt_path("ALPS.xsl"));
        out << alps::start_tag("SIMULATION") << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
            << alps::attribute("xsi:noNamespaceSchemaLocation","http://xml.comp-phys.org/2003/10/ALPS.xsd");
    
        out << parms;
    
        out << alps::start_tag("EIGENSTATES") << alps::attribute("number", neigen);
    
        for (int eig=0; eig<neigen; ++eig) {
            MPS<Matrix, SymmGroup> mps;
            load(checkpoints[eig], mps);
        
            maquis::cout << "Measurements." << std::endl;
            std::for_each(measurements.begin(), measurements.end(), measure_and_save<Matrix, SymmGroup>(rfile, "/spectrum/results/", mps, eig));
        
            std::vector<int> * measure_es_where = NULL;
            entanglement_spectrum_type * spectra = NULL;
            if (parms.defined("entanglement_spectra")) {
                spectra = new entanglement_spectrum_type();
                measure_es_where = new std::vector<int>();
                *measure_es_where = parms.template get<std::vector<int> >("entanglement_spectra");
            }
            std::vector<double> entropies, renyi2;
            if (parms["MEASURE[Entropy]"]) {
                std::cout << "Calculating vN entropy." << std::endl;
                entropies = calculate_bond_entropies(mps);
            }
            if (parms["MEASURE[Renyi2]"]) {
                std::cout << "Calculating n=2 Renyi entropy." << std::endl;
                renyi2 = calculate_bond_renyi_entropies(mps, 2, measure_es_where, spectra);
            }
        
            double energy = maquis::real(expval(mps, mpo));
            std::cout << "Energy: " << energy << std::endl;
            {
                storage::archive ar(rfile, "w");
                save_val_at_index(ar, "/spectrum/results/Energy/mean/value", energy, eig);
                if (entropies.size() > 0) save_val_at_index(ar, "/spectrum/results/Entropy/mean/value", entropies, eig);
                if (renyi2.size() > 0)    save_val_at_index(ar, "/spectrum/results/Renyi2/mean/value", renyi2, eig);
                if (spectra != NULL)      save_val_at_index(ar, "/spectrum/results/Entanglement Spectra/mean/value", *spectra, eig);
            }

            if (parms["MEASURE[EnergyVariance]"]) {
                MPO<Matrix, SymmGroup> mpo2 = square_mpo(mpo);
                mpo2.compress(1e-12);
            
                double energy2 = maquis::real(expval(mps, mpo2, true));
            
                maquis::cout << "Energy^2: " << energy2 << std::endl;
                maquis::cout << "Variance: " << energy2 - energy*energy << std::endl;
            
                {
                    storage::archive ar(rfile, "w");
                    save_val_at_index(ar, "/spectrum/results/Energy^2/mean/value", energy2, eig);
                    save_val_at_index(ar, "/spectrum/results/EnergyVariance/mean/value", energy2 - energy*energy, eig);
                }
            }
        
            /// Output into xml file
            out << alps::start_tag("EIGENSTATE") << alps::attribute("number", eig);
            if (write_xml)
                std::for_each(measurements.begin(), measurements.end(), save_to_xml(out));
        
            out << alps::start_tag("SCALAR_AVERAGE") <<  alps::attribute("name", "Energy") << alps::no_linebreak
                << alps::start_tag("MEAN") <<  alps::no_linebreak << energy << alps::end_tag("MEAN")
                << alps::end_tag("SCALAR_AVERAGE");
        
            out << alps::end_tag("EIGENSTATE");
        }
    
        out << alps::end_tag("EIGENSTATES");
        out << alps::end_tag("SIMULATION");
    }
}


#endif
