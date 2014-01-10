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
    
    using base::mps;
    using base::mpo;
    using base::lat;
    using base::parms;
    using base::model;
    using base::all_measurements;
    using base::sweep_measurements;
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
        
        parms << parms.iteration_params("t", init_sweep);

        this->model_init();
        this->mps_init();
        
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
            if (update_each > -1 && (sweep % update_each) == 0)
            {
                BaseParameters iteration_params = parms.iteration_params("t", sweep);
                if (iteration_params.size() > 0) {
                    parms <<iteration_params;
                    this->model_init(sweep);
                    meas_each    = parms["measure_each"];
                    chkp_each    = parms["chkp_each"];
                    update_each  = parms["update_each"];
                    evolver = TimeEvolver(&parms, &mps, lat, model, sweep);
                }
            } else if (sweep == nsweeps_img) {
                    // since this is just a change in the time step, there is
                    // no need to split the hamiltonian in non-overlapping terms.
                    evolver.prepare_te_terms();
            }
            
            /// time evolution
            evolver(nsteps);
            sweep = evolver.sweep();
            
            /// measurements
            if ((sweep+1) % meas_each == 0 || (sweep+1) == parms["nsweeps"])
            {
                /// measure energy
                double energy = maquis::real(expval(mps, mpo));
                maquis::cout << "Energy " << energy << std::endl;
                
                /// measure observables specified in 'always_measure'
                if (sweep_measurements.size() > 0)
                    this->measure(this->results_archive_path(sweep) + "/results/", sweep_measurements);

                /// write iteration results
                {
                    storage::archive ar(rfile, "w");
                    ar[this->results_archive_path(sweep) + "/parameters"] << parms;
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
    std::string results_archive_path(int sweep) const
    {
        status_type status;
        status["sweep"] = sweep;
        return base::results_archive_path(status);
    }
    
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
