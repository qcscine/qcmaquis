/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2019 by Leon Freitag <lefreita@ethz.ch>
 *               2020- by Robin Feldmann <robinfe@phys.chem.ethz.ch>
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

using chem::Hamiltonian;

namespace maquis {
    
// Types for measurement results
// meas_with_results_type: one measurement = pair of vectors with labels and results
template <class V>
using meas_with_results_type = std::pair<std::vector<std::vector<int> >, std::vector<V> >;

// All measurements -- a map with measurement names as keys and results as values
template <class V>
using results_map_type = std::map<std::string, meas_with_results_type<V> >;

template <class V, Hamiltonian HamiltonianType=Hamiltonian::Electronic> // real or complex
class DMRGInterface
{
public:
    typedef maquis::meas_with_results_type<V> meas_with_results_type;
    typedef maquis::results_map_type<V> results_map_type;

    /** @brief Class constructor */
    DMRGInterface(DmrgParameters& parms_);

    /** @brief Class destructor */
    ~DMRGInterface();

    /** @brief Run a DMRG optimization */
    void optimize();

    /** @brief Run a DMRG propagation */
    void evolve();

    /** @brief Gets the energy at the end of the simulation */
    V energy();

    /** @brief Getter for an obtject storing the statistics of the optimization */
    results_collector& get_iteration_results();
            
    /** @brief Gets how many sweeps have been run */
    int get_last_sweep();

    // Run dmrg_meas (measure all measurements and save them to a HDF5 file specified in parameters)
    // Do not store the measurements in a map.
    void run_measure();

    // Run all measurements and save them as a map
    // Currently, all measurements must be added via parameters first, and then measure() may be executed
    // In the future, we should support single measurements directly from the interface
    void measure();

    /** @brief Getter for the measurements */
    const results_map_type & measurements();

    /** @brief Updates the integrals and re-initialize the model */
    void update_integrals(const integral_map<V> & integrals);

    // Get RDMs
    // TODO: This does not work for 2U1/2U1PG symmetry because "oneptdm" measurement is not recognised by the model!
    // Fix the model to recognise it!
    const meas_with_results_type & onerdm();
    const meas_with_results_type & onespdm();
    const meas_with_results_type & twordm();
    const meas_with_results_type & threerdm();
    const meas_with_results_type & fourrdm();

    // Measure 3 and 4-RDM (for now in 2U1), and save it into the corresponding (SU2U1) result file (which should be set with parms["rfile"])
    void measure_and_save_3rdm();
    void measure_and_save_4rdm();

    /** @brief Getter for the mutual information */
    const meas_with_results_type& mutinf();

    // The same for transition 3-RDM. bra_name provides the name of the bra checkpoint
    // (doesn't matter if it is in SU2U1 or 2U1, as in the former case it will be transformed)
    void measure_and_save_trans3rdm(const std::string & bra_name);

    // Load an MPS from a given checkpoint and measure the overlap with the current MPS
    V overlap(const std::string& aux_mps_name);

    // Dump parameters into text file
    void dump_parameters(const std::string & file);

private:
    DmrgParameters& parms;
    results_map_type measurements_;
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // maquis

#endif
