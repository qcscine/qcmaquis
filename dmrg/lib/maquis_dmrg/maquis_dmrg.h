/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2019 by Leon Freitag <lefreita@ethz.ch>
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

#ifndef MAQUIS_DMRG_H
#define MAQUIS_DMRG_H

#include "maquis_dmrg_detail.h"


namespace maquis
{

    template <class V> // real or complex
    class DMRGInterface
    {
        public:
            // Types for measurement results
            // meas_with_results_type: one measurement = pair of vectors with labels and results
            typedef std::pair<std::vector<std::vector<int> >, std::vector<V> > meas_with_results_type;
            // All measurements -- a map with measurement names as keys and results as values
            typedef std::map<std::string, meas_with_results_type> results_map_type;

            DMRGInterface(DmrgParameters& parms_);
            ~DMRGInterface();

            // Run DMRG optimization
            void optimize();
            // Get energy after the optimization
            V energy();

            // Get sweep statistics
            results_collector& get_iteration_results();
            int get_last_sweep();

            // Run dmrg_meas (measure all measurements and save them to a HDF5 file specified in parameters)
            // Do not store the measurements in a map
            void run_measure();

            // Run all measurements and save them as a map
            // Currently, all measurements must be added via parameters first, and then measure() may be executed
            // In the future, we should support single measurements directly from the interface
            void measure();

            const results_map_type & measurements();

            //update the integrals and re-initialize the model
            void update_integrals(const integral_map<V> & integrals);

            // TODO: This does not work for 2U1/2U1PG symmetry because "oneptdm" measurement is not recognised by the model!
            // Fix the model to recognise it!
            const meas_with_results_type & onerdm();
            const meas_with_results_type & twordm();

            // Load an MPS from a given checkpoint and measure the overlap with the current MPS
            V overlap(const std::string& aux_mps_name);
        private:
            DmrgParameters& parms;
            results_map_type measurements_;

            struct Impl;
            std::unique_ptr<Impl> impl_;
    };
}
#endif
