/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_DMRG_TEVOL_SIM_H
#define APP_DMRG_TEVOL_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include "dmrg/sim/sim.h"
#include "dmrg/evolve/te_utils.hpp"

template <class Matrix, class SymmGroup, class TimeEvolver>
class tevol_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    typedef typename base::status_type status_type;
    typedef typename base::measurements_type measurements_type;

    using base::mps;
    using base::mpo;
    using base::lat;
    using base::parms;
    using base::model;
    using base::all_measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::rfile;
    
public:
    tevol_sim(DmrgParameters const & parms_, bool write_xml_)
    : base(parms_)
    , write_xml(write_xml_)
    {
        alps::oxstream out(boost::replace_last_copy(rfile, ".h5", ".xml"));
        out << alps::header("UTF-8") << alps::stylesheet(alps::xslt_path("ALPS.xsl"));
        out << alps::start_tag("SIMULATION") << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
            << alps::attribute("xsi:noNamespaceSchemaLocation","http://xml.comp-phys.org/2003/10/ALPS.xsd");
        out << parms;
        out << alps::end_tag("SIMULATION");
    }
    
    void run()
    {
        int meas_each    = parms["measure_each"];
        int chkp_each    = parms["chkp_each"];
        int update_each  = parms["update_each"];
        int nsweeps     = parms["nsweeps"];
        int nsweeps_img = parms["nsweeps_img"];
        
        parms << parms.iteration_params("Time", init_sweep);
        
        /// compute nsteps as the min of the three above
        int nsteps = parms["nsweeps"];
        if (meas_each > -1)
            nsteps = std::min(nsteps, meas_each);
        if (chkp_each > -1)
            nsteps = std::min(nsteps, chkp_each);
        if (update_each > -1)
            nsteps = std::min(nsteps, update_each);
        
        #define CHECK_MULTIPLICITY(var)                                               \
        if (var > 0 && var % nsteps != 0)                                             \
            throw std::runtime_error("var must be a multiple of 'nsteps'.");  
        CHECK_MULTIPLICITY(nsweeps)
        CHECK_MULTIPLICITY(nsweeps_img)
        CHECK_MULTIPLICITY(meas_each)
        CHECK_MULTIPLICITY(chkp_each)
        CHECK_MULTIPLICITY(update_each)
        #undef CHECK_MULTIPLICITY
        
        TimeEvolver evolver(&parms, &mps, lat, model, init_sweep);
        
        int n = nsweeps / nsteps;
        for (int i=init_sweep/nsteps; i < n; ++i) {
            // TODO: introduce some timings
            
            int sweep = i*nsteps;
            BaseParameters iteration_params = parms.iteration_params("Time", sweep);
            if (update_each > -1 && (sweep % update_each) == 0)
            {
                if (iteration_params.size() > 1) { // Time will always be set
                    parms << iteration_params;
                    meas_each    = parms["measure_each"];
                    chkp_each    = parms["chkp_each"];
                    update_each  = parms["update_each"];
                    model.update(parms);
                    mpo = make_mpo(lat, model, parms);
                    evolver = TimeEvolver(&parms, &mps, lat, model, sweep);
                }
            } else if (sweep == nsweeps_img) {
                    // since this is just a change in the time step, there is
                    // no need to split the hamiltonian in non-overlapping terms.
                    evolver.prepare_te_terms(sweep);
            }
            
            /// time evolution
            evolver(sweep, nsteps);
            sweep = (i+1)*nsteps - 1;
            iteration_params.set("Time", sweep);
            
            /// measurements
            if ((sweep+1) % meas_each == 0 || (sweep+1) == parms["nsweeps"])
            {
                /// measure energy
                double energy = maquis::real(expval(mps, mpo));
                maquis::cout << "Energy " << energy << std::endl;
                
                /// measure observables specified in 'always_measure'
                measurements_type always_measure = this->iteration_measurements(sweep); // todo: make measure() using const&
                if (!parms["ALWAYS_MEASURE"].empty())
                    this->measure(this->results_archive_path(sweep) + "/results/", always_measure);

                /// write iteration results
                {
                    storage::archive ar(rfile, "w");
                    ar[this->results_archive_path(sweep) + "/parameters"] << iteration_params;
                    ar[this->results_archive_path(sweep) + "/results"] << evolver.iteration_results();
                    ar[this->results_archive_path(sweep) + "/results/Energy/mean/value"] << std::vector<double>(1, energy);
                    // ar[this->results_archive_path(sweep) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
                }
            }
            
            /// write checkpoint
            bool stopped = stop_callback();
            if (stopped || (sweep+1) % chkp_each == 0 || (sweep+1) == parms["nsweeps"])
                checkpoint_simulation(mps, sweep);
            
            if (stopped) break;
        }
    }
    
private:
    using base::results_archive_path; // The following function is an overload, not the virtual function
    std::string results_archive_path(int sweep) const
    {
        status_type status;
        status["sweep"] = sweep;
        return base::results_archive_path(status);
    }
    
    using base::checkpoint_simulation; // The following function is an overload, not the virtual function
    void checkpoint_simulation(MPS<Matrix, SymmGroup> const& state, int sweep)
    {
        status_type status;
        status["sweep"] = sweep;
        status["site"]  = -1;
        return base::checkpoint_simulation(state, status);
    }
    
    bool write_xml;
};

#endif
