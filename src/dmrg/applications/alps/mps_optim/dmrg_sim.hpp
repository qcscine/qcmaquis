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

#ifndef ALPS_MPS_OPTIM_DMRG_SIM_HPP
#define ALPS_MPS_OPTIM_DMRG_SIM_HPP

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/shared_ptr.hpp>

#include "dmrg/sim/sim.h"
#include "dmrg/mp_tensors/optimize.h"
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
class dmrg_sim : public sim<Matrix, SymmGroup> {
    
    typedef sim<Matrix, SymmGroup> base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    typedef typename base::status_type status_type;
    
    using base::mps;
    using base::mpo;
    using base::mpoc;
    using base::parms;
    using base::all_measurements;
    using base::sweep_measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::init_site;
    using base::rfile;
    
public:
    
    dmrg_sim (DmrgParameters & parms_)
    : base(parms_)
    { }
    
    void run()
    {
        int meas_each = parms["measure_each"];
        int chkp_each = parms["chkp_each"];
        
        this->model_init();
        this->mps_init();
        
        maquis::cout << "MPS initialization has finished...\n"; // MPS restored now
        
        boost::shared_ptr<opt_base_t> optimizer;
        if (init_sweep < parms["nsweeps"]) {
            if (parms["optimization"] == "singlesite")
                optimizer.reset( new ss_optimize<Matrix, SymmGroup, storage::disk>
                                                (mps, mpoc, parms, stop_callback, init_site) );
            else if(parms["optimization"] == "twosite")
                optimizer.reset( new ts_optimize<Matrix, SymmGroup, storage::disk>
                                                (mps, mpoc, parms, stop_callback, init_site) );
            else
                throw std::runtime_error("Don't know this optimizer");
        }
        
        try {
            bool stopped = false;
            for (int sweep=init_sweep; sweep < parms["nsweeps"]; ++sweep) {
                // TODO: introduce some timings
                
                optimizer->sweep(sweep, Both);
                storage::disk::sync();
                
                if ((sweep+1) % meas_each == 0 || (sweep+1) == parms["nsweeps"])
                {
                    /// write iteration results
                    {
                        storage::archive ar(rfile, "w");
                        ar[results_archive_path(sweep) + "/parameters"] << parms;
                        ar[results_archive_path(sweep) + "/results"] << optimizer->iteration_results();
                        // ar[results_archive_path(sweep) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
                    }
                    
                    /// measure observables specified in 'always_measure'
                    if (sweep_measurements.size() > 0)
                        this->measure(this->results_archive_path(sweep) + "/results/", sweep_measurements);
                }
                
                /// write checkpoint
                stopped = stop_callback();
                if (stopped || (sweep+1) % chkp_each == 0 || (sweep+1) == parms["nsweeps"])
                    checkpoint_simulation(mps, sweep, -1);
                
                if (stopped) break;
            }
            
            /// Final measurements
            if (!stopped) {
                this->measure("/spectrum/results/", all_measurements);
                
                alps::oxstream out(boost::replace_last_copy(rfile, ".h5", ".xml"));
                out << alps::header("UTF-8") << alps::stylesheet(alps::xslt_path("ALPS.xsl"));
                out << alps::start_tag("SIMULATION") << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
                    << alps::attribute("xsi:noNamespaceSchemaLocation","http://xml.comp-phys.org/2002/10/ALPS.xsd");

                out << parms;
                
                out << alps::start_tag("EIGENSTATES") << alps::attribute("number", 1);
                out << alps::start_tag("EIGENSTATE") << alps::attribute("number", 0);
                
                std::for_each(all_measurements.begin(), all_measurements.end(), save_to_xml(out));

                double energy = maquis::real(expval(mps, mpoc));
                maquis::cout << "Energy: " << maquis::real(expval(mps, mpoc)) << std::endl;
                {
                    storage::archive ar(rfile, "w");
                    ar["/spectrum/results/Energy/mean/value"] << std::vector<double>(1, energy);
                }
                out << alps::start_tag("SCALAR_AVERAGE") <<  alps::attribute("name", "Energy") << alps::no_linebreak
                    << alps::start_tag("MEAN") <<  alps::no_linebreak << energy << alps::end_tag("MEAN")
                    << alps::end_tag("SCALAR_AVERAGE");

                out << alps::end_tag("EIGENSTATE");
                out << alps::end_tag("EIGENSTATES");


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
                out << alps::end_tag("SIMULATION");
            }
            
        } catch (dmrg::time_limit const& e) {
            maquis::cout << e.what() << " checkpointing partial result." << std::endl;
            checkpoint_simulation(mps, e.sweep(), e.site());
            
            {
                storage::archive ar(rfile, "w");
                ar[results_archive_path(e.sweep()) + "/parameters"] << parms;
                ar[results_archive_path(e.sweep()) + "/results"] << optimizer->iteration_results();
                // ar[results_archive_path(e.sweep()) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
            }
        }
    }
    
    ~dmrg_sim()
    {
        storage::disk::sync();
    }
    
private:
    std::string results_archive_path(int sweep) const
    {
        status_type status;
        status["sweep"] = sweep;
        return base::results_archive_path(status);
    }
    
    void checkpoint_simulation(MPS<Matrix, SymmGroup> const& state, int sweep, int site)
    {
        status_type status;
        status["sweep"] = sweep;
        status["site"]  = site;
        return base::checkpoint_simulation(state, status);
    }
    
};


#endif
