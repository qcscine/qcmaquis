/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018         Leon Freitag <lefreita@ethz.ch>
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
#ifndef INTERFACE_SIM_H
#define INTERFACE_SIM_H
#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include "dmrg/sim/sim.h"
#include "dmrg/optimize/optimize.h"
#include "dmrg/models/chem/measure_transform.hpp"


// The sim class for interface-based DMRG runs and measurements
template <class Matrix, class SymmGroup>
class interface_sim : public sim<Matrix, SymmGroup>, public abstract_interface_sim<Matrix> {

    typedef sim<Matrix, SymmGroup> base;
    typedef abstract_interface_sim<Matrix> interface_base;
    typedef optimizer_base<Matrix, SymmGroup, storage::disk> opt_base_t;
    typedef typename base::status_type status_type;
    typedef typename base::measurements_type measurements_type;
    typedef typename interface_base::meas_with_results_type meas_with_results_type;
    typedef typename interface_base::results_map_type results_map_type;

    using base::mps;
    using base::mpo;
    using base::parms;
    using base::all_measurements;
    using base::stop_callback;
    using base::init_sweep;
    using base::init_site;
    using base::rfile;

public:

    interface_sim (DmrgParameters & parms_)
    : base(parms_)
    { }

    void run()
    {
        optimize();
        // or whatever
        // throw std::runtime_error("run() shouldn't be called from interface_sim");
    }

    void optimize()
    {
        int meas_each = parms["measure_each"];
        int chkp_each = parms["chkp_each"];

        /// MPO creation
        if (parms["MODEL"] == std::string("quantum_chemistry") && parms["use_compressed"])
            throw std::runtime_error("chem compression has been disabled");
        MPO<Matrix, SymmGroup> mpoc = mpo;
        if (parms["use_compressed"])
            mpoc.compress(1e-12);

        /// Optimizer initialization
        std::shared_ptr<opt_base_t> optimizer;
        if (parms["optimization"] == "singlesite")
        {
            optimizer.reset( new ss_optimize<Matrix, SymmGroup, storage::disk>
                            (mps, mpoc, parms, stop_callback, init_site) );
        }
        else if(parms["optimization"] == "twosite")
        {
            optimizer.reset( new ts_optimize<Matrix, SymmGroup, storage::disk>
                            (mps, mpoc, parms, stop_callback, init_site) );
        }
        else {
            throw std::runtime_error("Don't know this optimizer");
        }

        measurements_type always_measurements = this->iteration_measurements(init_sweep);

        try {
            for (int sweep=init_sweep; sweep < parms["nsweeps"]; ++sweep) {
                // TODO: introduce some timings

                optimizer->sweep(sweep, Both);
                storage::disk::sync();

                bool converged = false;

                if (!rfile.empty())
                {
                    if ((sweep+1) % meas_each == 0 || (sweep+1) == parms["nsweeps"])
                    {
                        /// write iteration results
                        {
                            storage::archive ar(rfile, "w");
                            ar[results_archive_path(sweep) + "/parameters"] << parms;
                            ar[results_archive_path(sweep) + "/results"] << optimizer->iteration_results();
                            // ar[results_archive_path(sweep) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);

                            // stop simulation if an energy threshold has been specified
                            // FIXME: this does not work for complex numbers - stknecht feb 2016
                            int prev_sweep = sweep - meas_each;
                            if (prev_sweep >= 0 && parms["conv_thresh"] > 0.)
                            {
                                typedef typename maquis::traits::real_type<Matrix>::type real_type;
                                std::vector<real_type> energies;

                                ar[results_archive_path(sweep) + "/results/Energy/mean/value"] >> energies;
                                real_type emin = *std::min_element(energies.begin(), energies.end());
                                ar[results_archive_path(prev_sweep) + "/results/Energy/mean/value"] >> energies;
                                real_type emin_prev = *std::min_element(energies.begin(), energies.end());
                                real_type e_diff = std::abs(emin - emin_prev);

                                if (e_diff < parms["conv_thresh"])
                                    converged = true;
                            }
                        }

                        /// measure observables specified in 'always_measure'
                        if (always_measurements.size() > 0)
                            this->measure(this->results_archive_path(sweep) + "/results/", always_measurements);
                    }
                }

                /// write checkpoint
                bool stopped = stop_callback() || converged;
                if (stopped || (sweep+1) % chkp_each == 0 || (sweep+1) == parms["nsweeps"])
                    checkpoint_simulation(mps, sweep, -1);

                if (stopped) break;
            }
        } catch (dmrg::time_limit const& e) {
            maquis::cout << e.what() << " checkpointing partial result." << std::endl;
            checkpoint_simulation(mps, e.sweep(), e.site());

            {
                if (!rfile.empty())
                {
                    storage::archive ar(rfile, "w");
                    ar[results_archive_path(e.sweep()) + "/parameters"] << parms;
                    ar[results_archive_path(e.sweep()) + "/results"] << optimizer->iteration_results();
                    // ar[results_archive_path(e.sweep()) + "/results/Runtime/mean/value"] << std::vector<double>(1, elapsed_sweep + elapsed_measure);
                }
            }
        }
    }

    void run_measure()
    {
        this->measure("/spectrum/results/", all_measurements);

        /// MPO creation
        MPO<Matrix, SymmGroup> mpoc = mpo;
        if (parms["use_compressed"])
            mpoc.compress(1e-12);

        double energy;

        if (parms["MEASURE[Energy]"]) {
            energy = maquis::real(expval(mps, mpoc)) + maquis::real(mpoc.getCoreEnergy());
            maquis::cout << "Energy: " << energy << std::endl;

            if (!rfile.empty())
            {
                storage::archive ar(rfile, "w");
                ar["/spectrum/results/Energy/mean/value"] << std::vector<double>(1, energy);
            }
        }

        if (parms["MEASURE[EnergyVariance]"] > 0) {
            MPO<Matrix, SymmGroup> mpo2 = square_mpo(mpoc);
            mpo2.compress(1e-12);

            if (!parms["MEASURE[Energy]"]) energy = maquis::real(expval(mps, mpoc)) + maquis::real(mpoc.getCoreEnergy());
            double energy2 = maquis::real(expval(mps, mpo2, true));

            maquis::cout << "Energy^2: " << energy2 << std::endl;
            maquis::cout << "Variance: " << energy2 - energy*energy << std::endl;

            if (!rfile.empty())
            {
                storage::archive ar(rfile, "w");
                ar["/spectrum/results/Energy^2/mean/value"] << std::vector<double>(1, energy2);
                ar["/spectrum/results/EnergyVariance/mean/value"] << std::vector<double>(1, energy2 - energy*energy);
            }
        }

        #if defined(HAVE_TwoU1) || defined(HAVE_TwoU1PG)
        if (parms.is_set("MEASURE[ChemEntropy]"))
            if (!rfile.empty())
                measure_transform<Matrix, SymmGroup>()(rfile, "/spectrum/results", base::lat, mps);
            else
                throw std::runtime_error("Transformed measurements not implemented yet without checkpoints");
        #endif
    }

    results_map_type measure_out()
    {
        results_map_type ret;

        // Run all measurements and fill the result map
        for(auto&& meas: all_measurements)
            ret[meas.name()] = measure_and_save<Matrix,SymmGroup>(mps).meas_out(meas);

        #if defined(HAVE_TwoU1) || defined(HAVE_TwoU1PG)
        if (parms.is_set("MEASURE[ChemEntropy]"))
        {
            // Obtain a map with transformed measurements
            results_map_type transformed_meas = measure_transform<Matrix, SymmGroup>().meas_out(base::lat, mps);

            // Merge transformed measurements with the remaining results
            ret.insert(transformed_meas.begin(), transformed_meas.end());
        }
        #endif

        return ret;
    }

    typename Matrix::value_type get_energy()
    {
        return expval(mps, mpo) + mpo.getCoreEnergy();
    }

    ~interface_sim()
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