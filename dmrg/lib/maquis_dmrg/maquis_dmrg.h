/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_H
#define MAQUIS_DMRG_H

#include "maquis_dmrg_detail.h"

using chem::Hamiltonian;

namespace maquis {
    
// Types for measurement results
// meas_with_results_type: one measurement = pair of vectors with labels and results
template <typename ScalarType>
using meas_with_results_type = std::pair<std::vector<std::vector<int> >, std::vector<ScalarType> >;

// All measurements -- a map with measurement names as keys and results as values
template <typename ScalarType>
using results_map_type = std::map<std::string, meas_with_results_type<ScalarType> >;

template <typename ScalarType> // real or complex
class DMRGInterface
{
public:
    typedef maquis::meas_with_results_type<ScalarType> meas_with_results_type;
    typedef maquis::results_map_type<ScalarType> results_map_type;

    /** @brief Class constructor */
    explicit DMRGInterface(DmrgParameters& parms_);

    /** @brief Class destructor */
    ~DMRGInterface();

    /** @brief Run a DMRG optimization */
    void optimize();

    /** @brief Run a DMRG propagation */
    void evolve();

    /** @brief Run a DMRG[IPI] calculation */
    void runInversePowerIteration();

    /** @brief Run a DMRG[FEAST] calculation */
    void runFEAST();

    /** @brief Gets the energy at the end of the simulation */
    ScalarType energy();

    /** @brief Gets the energy at the end of the FEAST simulation */
    ScalarType energyFEAST(int iState);

    /** @brief Gets the overlap of the MPS with a given determinant */
    ScalarType getCICoefficient(std::string determinantString);

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
    void update_integrals(const integral_map<ScalarType> & integrals);

    // Get RDMs
    // TODO: This does not work for 2U1/2U1PG symmetry because "oneptdm" measurement is not recognised by the model!
    // Fix the model to recognise it!
    const meas_with_results_type & onerdm();
    const meas_with_results_type & onespdm();
    const meas_with_results_type & twordm();
    const meas_with_results_type & threerdm();
    const meas_with_results_type & fourrdm();
    const meas_with_results_type & getMeasurement(std::string measName);

    // Measure 3 and 4-RDM (for now in 2U1), and save it into the corresponding (SU2U1) result file (which should be set with parms["rfile"])
    void measure_and_save_3rdm();
    void measure_and_save_4rdm();

    /** @brief Getter for the mutual information */
    const meas_with_results_type& mutinf();

    // The same for transition 3-RDM. bra_name provides the name of the bra checkpoint
    // (doesn't matter if it is in SU2U1 or 2U1, as in the former case it will be transformed)
    void measure_and_save_trans3rdm(const std::string & bra_name);

    // Load an MPS from a given checkpoint and measure the overlap with the current MPS
    ScalarType overlap(const std::string& aux_mps_name);

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
