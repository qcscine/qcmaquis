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

            // Run DMRG optimization
            void optimize();

            // Run dmrg_meas (measure all measurements and save them to a HDF5 file specified in parameters)
            void run_measure();

            // Get energy after the optimization
            V energy();

            // Run all measurements and return them as a map
            // Currently, all measurements must be added via parameters first, and then measure() may be executed
            // In the future, we should support single measurements directly from the interface
            results_map_type measure();
        private:
            DmrgParameters& parms;
            typename simulation_traits<V>::shared_ptr sim;
    };
}
#endif
